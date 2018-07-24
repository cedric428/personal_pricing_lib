/*
 * Haocheng Wang, July 09, 2018
 * Last modification:
 * 2018-07-10: Add payoff smoothing for LVSBI
*/

#include <cmath>
#include <algorithm>

#include "utils/exceptions.h"

#include "payoff.h"
#include "fd.h"
#include "spline.h"

using namespace std;

namespace wst {
    namespace analytics {

        BaseBI::~BaseBI() {}

        void BaseBI::add_knockout(double ko, bool is_up, double start_time, double end_time, double rebate) {
            if (time_ >= 0)
                THROW(utils::ValueError, "Cannot add knockouts after the first payout has been added to the grid");
            if (ko <= 0) THROW(utils::ValueError, "Knockout must be positive, got ko=" << ko);

            ko_levels_.push_back(ko);
            ko_is_ups_.push_back(is_up);
            ko_start_times_.push_back(start_time);
            ko_end_times_.push_back(end_time);
            ko_rebates_.push_back(rebate);
        }

        vector<double> get_ko_break_times(const vector<bool> &ko_is_ups,
                                          const vector<double> &ko_start_times, const vector<double> &ko_end_times,
                                          double time) {
            // helper fn to track the knockout break times - we may need to resize the grid at
            // those points.

            bool found_it;
            vector<double> ko_break_times;
            for (size_t i = 0; i < ko_is_ups.size(); i++) {
                vector<double> ko_times;
                if (ko_start_times[i] > 0 and ko_start_times[i] < time) ko_times.push_back(ko_start_times[i]);
                if (ko_end_times[i] > 0 and ko_end_times[i] < time) ko_times.push_back(ko_end_times[i]);
                if (ko_times.size() == 0) continue; // nothing to track

                for (size_t j = 0; j < ko_times.size(); j++) {
                    // add the ko time if it's not already in the list
                    found_it = false;
                    for (size_t k = 0; k < ko_break_times.size(); k++)
                        if (ko_break_times[k] == ko_times[j]) {
                            found_it = true;
                            break;
                        }
                    if (not found_it)
                        ko_break_times.push_back(ko_times[j]);
                };
            }

            // sort the knockout break times in ascending order

            sort(ko_break_times.begin(), ko_break_times.end());
            return ko_break_times;
        }

        LVSBIBase::LVSBIBase(double spot,
                             const vector<double> &volsdd, const vector<double> &volsd,
                             const vector<double> &vols0,
                             const vector<double> &volsu, const vector<double> &volsuu,
                             const vector<double> &rds,
                             const vector<double> &rfs,
                             const vector<double> &break_times,
                             double alpha, double beta, double extrap_fact,
                             int nu, int nt, double nsd)
                : spot_(spot), volsdd_(volsdd), volsd_(volsd), vols0_(vols0), volsu_(volsu), volsuu_(volsuu),
                  rds_(rds), rfs_(rfs), break_times_(break_times), alpha_(alpha), beta_(beta),
                  extrap_fact_(extrap_fact),
                  nu_(nu), nt_(nt), nsd_(nsd) {
            time_ = -1; // notes that we haven't set a payoff yet

            // initialize a subset of grid stuff that doesn't change when time changes or
            // grid resizings happen.

            // figure out the layer size and transition probability frequencies

            double gamma = alpha / 2. / sqrt(beta);
            eps_ = 2 * asinh(gamma);
            double ee = exp(eps_);

            p0p_ = beta * (0.5 - gamma / 2. / sqrt(1 + gamma * gamma));
            p0m_ = beta * (0.5 + gamma / 2. / sqrt(1 + gamma * gamma));
            pp0_ = beta;
            pm0_ = beta;

            // For a grid vol level, use the maximum of all the local vol
            // parameters in the high vol state.

            max_vol_ = 0.;
            for (size_t i = 0; i < volsdd.size(); i++) {
                if (volsdd[i] > max_vol_) max_vol_ = volsdd[i];
                if (volsd[i] > max_vol_) max_vol_ = volsd[i];
                if (vols0[i] > max_vol_) max_vol_ = vols0[i];
                if (volsu[i] > max_vol_) max_vol_ = volsu[i];
                if (volsuu[i] > max_vol_) max_vol_ = volsuu[i];
            }

            max_vol_ *= sqrt(ee); // want high vol state vol

            // calculate the probabilities of being in the different vol layers. We assume
            // we're in a stationary distribution at the start (ie the initial vol layer
            // is uncertain too). We only track the probabilities of being in the up and down
            // states because the probability of being in the mid state is 1-those two.

            prob_up_ = 0.25 - gamma / 4. / sqrt(1 + gamma * gamma);
            prob_dn_ = 0.25 + gamma / 4. / sqrt(1 + gamma * gamma);
        }

        LVSBIBase::~LVSBIBase() {}

        void LVSBIBase::initialize_payoff(const Payoff &p, double time) {
            if (time_ >= 0) THROW(utils::Error, "Cannot initialize grid more than once");

            time_ = time;
            us_.resize(nu_); // initialized in setup_grid
            setup_grid(true);
            vals_.resize(nu_);
            vals_up_.resize(nu_);
            vals_dn_.resize(nu_);
            vals_mixed_.resize(nu_); // populated in go_back
            for (int i = 0; i < nu_; i++) vals_[i] = 0; // initialize with zero values in the grid
            add_payoff(p);
        }

        void LVSBIBase::setup_grid(bool recalcdt) {
            double lf = log(this->forward() / spot_);

            double umin = fmin(lf, 0) - nsd_ * max_vol_ * sqrt(time_);
            double umax = fmax(lf, 0) + nsd_ * max_vol_ * sqrt(time_);

            // if there are any knockouts in force at this time, check whether
            // they should cut off the grid. Note: we might end up (in extreme cases)
            // in a situation where umin>=umax due to adding a knockout that's way
            // in the money. We'll handle that case in the backward induction.

            max_u_is_ko_ = false; // overwritten if otherwise
            min_u_is_ko_ = false; // ditto
            max_u_rebate_ = 0;     // overwritten if there's a rebate
            min_u_rebate_ = 0;     // ditto
            double lk;
            for (size_t i = 0; i < ko_levels_.size(); i++) {
                lk = log(ko_levels_[i] / spot_);

                // check if the knockout is in force at the grid time

                if (ko_start_times_[i] >= time_) continue;
                if (ko_end_times_[i] >= 0 and ko_end_times_[i] < time_) continue;

                // if so, check whether it's inside the current max/min

                if (ko_is_ups_[i] and lk < umax) {
                    umax = lk;
                    max_u_is_ko_ = true;
                    max_u_rebate_ = ko_rebates_[i];
                }
                if (not ko_is_ups_[i] and lk > umin) {
                    umin = lk;
                    min_u_is_ko_ = true;
                    min_u_rebate_ = ko_rebates_[i];
                }
            }

            du_ = (umax - umin) / (nu_ - 1);

            if (recalcdt) dt_ = time_ / nt_;

            // initialize the vector of u values (sized properly already)

            for (int i = 0; i < nu_; i++) us_[i] = umin + i * du_;
        }

        double LVSBIBase::interp(double interp_spot) const {
            double u = log(interp_spot / spot_);

            // respect knockouts

            if (u <= us_[0] and min_u_is_ko_) return min_u_rebate_;
            if (u >= us_[nu_ - 1] and max_u_is_ko_) return max_u_rebate_;

            // otherwise interpolate

            if (alpha_ > 0) {
                NaturalSpline s(us_, vals_mixed_);
                return s.value(u);
            } else {
                NaturalSpline s(us_, vals_);
                return s.value(u);
            }
        }

        void LVSBIBase::add_payoff(const Payoff &p) {
            // add the payoff onto the grid. Respect any knockouts by zeroing the value
            // at the edges if appropriate. However: that works only for knockouts that start
            // before the current time. Convention in that case: if a knockout starts on a particular
            // date and that happens to be an expiration date, we'll *ignore* it, not treat it
            // as a European ko. If you want to treat it that way, add the knockout start time
            // as a teeny time right before the expiration date.

            double val;
            for (int i = 0; i < nu_; i++) {
                if (i == 0 and min_u_is_ko_)
                    vals_[i] = vals_up_[i] = vals_dn_[i] = min_u_rebate_;
                else if (i == nu_ - 1 and max_u_is_ko_)
                    vals_[i] = vals_up_[i] = vals_dn_[i] = max_u_rebate_;
                else {
                    val = p.value(spot_ * exp(us_[i]));
                    vals_[i] += val;
                    vals_up_[i] += val;
                    vals_dn_[i] += val;
                }
            }
        }

        void LVSBIBase::multiply_payoff(const Payoff &p) {
            // multiply grid values by the payoff value

            double v;
            for (int i = 0; i < nu_; i++) {
                v = p.value(spot_ * exp(us_[i]));
                vals_[i] *= v;
                vals_up_[i] *= v;
                vals_dn_[i] *= v;
            }
        }

        void LVSBIBase::max_payoff(const Payoff &p) {
            // max grid values with the payoff value

            double v;
            for (int i = 0; i < nu_; i++) {
                v = p.value(spot_ * exp(us_[i]));
                vals_[i] = fmax(vals_[i], v);
                vals_up_[i] = fmax(vals_up_[i], v);
                vals_dn_[i] = fmax(vals_dn_[i], v);
            }
        }

        double LVSBIBase::forward() const {
            double fwd = spot_;
            double last_time, end_time, rate_time;
            int n = break_times_.size();

            for (int i = 0; i < n; i++) {
                if (i == 0)
                    last_time = 0;
                else
                    last_time = break_times_[i - 1];
                rate_time = break_times_[i];
                if (rate_time < time_)
                    end_time = rate_time;
                else
                    end_time = time_;

                fwd *= exp((rds_[i] - rfs_[i]) * (end_time - last_time));
                if (end_time > time_ - 1e-12) break;
            }

            if (time_ > break_times_[n - 1])
                fwd *= exp((rds_[n - 1] - rfs_[n - 1]) * (time_ - break_times_[n - 1]));

            return fwd;
        }

        LVSBIExp::LVSBIExp(double spot,
                           const vector<double> &volsdd, const vector<double> &volsd,
                           const vector<double> &vols0,
                           const vector<double> &volsu, const vector<double> &volsuu,
                           const vector<double> &rds,
                           const vector<double> &rfs,
                           const vector<double> &break_times,
                           double alpha, double beta, double extrap_fact,
                           int nu, int nt, double nsd)
                : LVSBIBase(spot, volsdd, volsd, vols0, volsu, volsuu, rds, rfs, break_times, alpha, beta, extrap_fact,
                            nu, nt, nsd) {
        }

        LVSBIExp::~LVSBIExp() {}

        void extrap_vals(vector<double> &vals, int nu, bool min_u_is_ko, bool max_u_is_ko, double min_u_rebate_,
                         double max_u_rebate_) {
            double dvdu;
            if (min_u_is_ko)
                vals[0] = min_u_rebate_;
            else {
                dvdu = vals[2] - vals[1];
                vals[0] = vals[1] - dvdu;
            }

            if (max_u_is_ko)
                vals[nu - 1] = max_u_rebate_;
            else {
                dvdu = vals[nu - 2] - vals[nu - 3];
                vals[nu - 1] = vals[nu - 2] + dvdu;
            }
        }

        void LVSBIExp::go_back(double new_time) {
            if (time_ < 0) THROW(utils::Error, "Grid not initialized yet");
            if (new_time > time_)
                THROW(utils::ValueError, "Cannot backward induct to future time=" << new_time << ", time=" << time_);
            if (new_time == time_) return; // nothing to do

            vector<double> xs(5);
            xs[0] = -1.28;
            xs[1] = -0.68;
            xs[2] = 0;
            xs[3] = 0.68;
            xs[4] = 1.28;
            vector<double> vols(5);

            double ee = exp(eps_);
            double stepdt, next_time, voldd = 0, vold = 0, vol0 = 0, volu = 0, voluu = 0, rd = 0, rf = 0, vol, sqrtT, disc;
            int last_break_ind = -1, next_ind, break_ind;

            int np = break_times_.size();
            FlatVolSpline *sp = 0;

            vector<double> grid_vol_sqs(nu_);
            vector<double> new_vals(nu_), new_vals_up(nu_), new_vals_dn(nu_);

            double P0p, P0m, Pp0, Pm0;

            // track the knockout break times - we may need to resize the grid at
            // those points.

            vector<double> ko_break_times(get_ko_break_times(ko_is_ups_, ko_start_times_, ko_end_times_, time_));

            // track the ko break time index that we're working on and the ko break time corresponding
            // to it

            int ko_break_ind = ko_break_times.size() - 1;
            double ko_break_time = -1;
            if (ko_break_ind >= 0) ko_break_time = ko_break_times[ko_break_ind];

            double driftT = log(this->forward() / spot_);

            while (true) {
                // figure out the next time step to move back to. Land on an integer multiple of
                // dt_ unless that takes us_ past new_time, or past the latest knockout break time.

                next_ind = int(time_ / dt_);
                next_time = next_ind * dt_;
                if (fabs(next_time - time_) < 1e-12) next_time -= dt_;
                if (next_time < new_time) next_time = new_time;
                if (next_time < ko_break_time) next_time = ko_break_time;

                stepdt = time_ - next_time; // not always this->dt_

                // figure out the parameters to use for this interval - use the ones at the
                // end of the time interval

                vector<double>::const_iterator it = lower_bound(break_times_.begin(), break_times_.end(), time_);
                if (it == break_times_.end())
                    break_ind = np - 1; // extrapolate flat
                else
                    break_ind = it - break_times_.begin();

                if (break_ind != last_break_ind) {
                    voldd = volsdd_[break_ind];
                    vold = volsd_[break_ind];
                    vol0 = vols0_[break_ind];
                    volu = volsu_[break_ind];
                    voluu = volsuu_[break_ind];

                    if (sp != 0) delete sp; // remove any old spline

                    if (fabs(voldd - vold) < 1e-12 and fabs(voldd - vol0) < 1e-12 and fabs(voldd - volu) < 1e-12 and
                        fabs(voldd - voluu) < 1e-12) {
                        sp = 0; // zero out the pointer as a signal not to use a spline
                    } else {
                        vols[0] = voldd;
                        vols[1] = vold;
                        vols[2] = vol0;
                        vols[3] = volu;
                        vols[4] = voluu;
                        sp = new FlatVolSpline(xs, vols, extrap_fact_);
                    }

                    rd = rds_[break_ind];
                    rf = rfs_[break_ind];
                }

                // get local vols for all grid points - use values at the end of the time interval.
                // If there's no spline defined it's because the vols are flat; just use the vol0 value
                // everywhere.

                if (sp != 0) sqrtT = sqrt(time_);

                for (int i = 0; i < nu_; i++) {
                    if (sp == 0)
                        grid_vol_sqs[i] = vol0 * vol0;
                    else {
                        // x-value for spline = u-value/vol0/sqrtT.
                        vol = sp->value(us_[i] / vol0 / sqrtT);
                        grid_vol_sqs[i] = vol * vol;
                    }
                }

                // do the backward induction

                for (int i = 1; i < nu_ - 1; i++) {
                    new_vals[i] = vals_[i]
                                  + stepdt * (
                            grid_vol_sqs[i] / 2. * (vals_[i + 1] + vals_[i - 1] - 2 * vals_[i]) / du_ / du_
                            + (rd - rf - grid_vol_sqs[i] / 2.) * (vals_[i + 1] - vals_[i - 1]) / 2. / du_);
                    new_vals_up[i] = vals_up_[i]
                                     + stepdt * (
                            grid_vol_sqs[i] * ee / 2. * (vals_up_[i + 1] + vals_up_[i - 1] - 2 * vals_up_[i]) / du_ /
                            du_
                            + (rd - rf - grid_vol_sqs[i] * ee / 2.) * (vals_up_[i + 1] - vals_up_[i - 1]) / 2. / du_);
                    new_vals_dn[i] = vals_dn_[i]
                                     + stepdt * (
                            grid_vol_sqs[i] / ee / 2. * (vals_dn_[i + 1] + vals_dn_[i - 1] - 2 * vals_dn_[i]) / du_ /
                            du_
                            + (rd - rf - grid_vol_sqs[i] / ee / 2.) * (vals_dn_[i + 1] - vals_dn_[i - 1]) / 2. / du_);
                }

                // extrapolate linearly for edge points, unless they're at a knockout, in which
                // case the value is zero.

                extrap_vals(new_vals, nu_, min_u_is_ko_, max_u_is_ko_, min_u_rebate_, max_u_rebate_);
                extrap_vals(new_vals_up, nu_, min_u_is_ko_, max_u_is_ko_, min_u_rebate_, max_u_rebate_);
                extrap_vals(new_vals_dn, nu_, min_u_is_ko_, max_u_is_ko_, min_u_rebate_, max_u_rebate_);

                // mix the layers and apply discounting

                disc = exp(-rd * stepdt);
                P0p = p0p_ * stepdt;
                P0m = p0m_ * stepdt;
                Pp0 = pp0_ * stepdt;
                Pm0 = pm0_ * stepdt;
                for (int i = 0; i < nu_; i++) {
                    if (i == 0 and min_u_is_ko_)
                        vals_[i] = vals_up_[i] = vals_dn_[i] = min_u_rebate_;
                    else if (i == nu_ - 1 and max_u_is_ko_)
                        vals_[i] = vals_up_[i] = vals_dn_[i] = max_u_rebate_;
                    else {
                        vals_[i] = disc * (new_vals[i] * (1 - P0p - P0m) + new_vals_up[i] * P0p + new_vals_dn[i] * P0m);
                        vals_up_[i] = disc * (new_vals_up[i] * (1 - Pp0) + new_vals[i] * Pp0);
                        vals_dn_[i] = disc * (new_vals_dn[i] * (1 - Pm0) + new_vals[i] * Pm0);
                    }
                }

                // move the grid time back and stop if we're at the end

                time_ = next_time;
                if (fabs(time_ - new_time) < 1e-12) break;

                // adjust the accumulated risk neutral drift * time to the new time

                driftT -= (rd - rf) * stepdt;

                // figure out whether to resize the grid. Definitely do it if we're at a knockout
                // break time; otherwise check whether the grid is too narrow based on the forward
                // to this date and the std dev.

                bool resize_grid = false;

                // are we at a knockout break? That requires a resize. Also note the
                // subsequent ko break time to track.

                if (fabs(time_ - ko_break_time) < 1e-12) {
                    ko_break_ind -= 1;
                    if (ko_break_ind == -1)
                        ko_break_time = -1;
                    else
                        ko_break_time = ko_break_times[ko_break_ind];
                    resize_grid = true;
                }

                if (resize_grid) {
                    // construct a spline off the current us_ and grid values for all
                    // three layers. Note these are natural splines with 2nd deriv going to
                    // zero at the edges - not the same as the "flatVolSpline" spline we use
                    // for local vols.

                    NaturalSpline sUp(us_, vals_up_);
                    NaturalSpline sDn(us_, vals_dn_);
                    NaturalSpline s(us_, vals_);

                    // setup the grid again - this properly resizes for new knockouts etc.
                    // Don't reset dt_ though.

                    setup_grid(false);

                    // change the grid values to the interpolated ones

                    for (int i = 0; i < nu_; i++) {
                        vals_[i] = s.value(us_[i]);
                        vals_up_[i] = sUp.value(us_[i]);
                        vals_dn_[i] = sDn.value(us_[i]);
                    }

                    // always respect the knockout values if appropriate

                    if (min_u_is_ko_) vals_[0] = vals_up_[0] = vals_dn_[0] = min_u_rebate_;
                    if (max_u_is_ko_) vals_[nu_ - 1] = vals_up_[nu_ - 1] = vals_dn_[nu_ - 1] = max_u_rebate_;
                }
            }

            // if we need to delete the memory for the spline, do so now

            if (sp != 0) delete sp;

            // generated the mixed layer

            for (int i = 0; i < nu_; i++)
                vals_mixed_[i] = vals_[i] * (1 - prob_up_ - prob_dn_) + vals_up_[i] * prob_up_ + vals_dn_[i] * prob_dn_;
        }

        LVSBICN::LVSBICN(double spot,
                         const vector<double> &volsdd, const vector<double> &volsd,
                         const vector<double> &vols0,
                         const vector<double> &volsu, const vector<double> &volsuu,
                         const vector<double> &rds,
                         const vector<double> &rfs,
                         const vector<double> &break_times,
                         double alpha, double beta, double extrap_fact,
                         int nu, int nt, double nsd,
                         double theta)
                : LVSBIBase(spot, volsdd, volsd, vols0, volsu, volsuu, rds, rfs, break_times, alpha, beta, extrap_fact,
                            nu, nt, nsd),
                  theta_(theta) {
        }

        LVSBICN::~LVSBICN() {}

        void cnbi(vector<double> &new_vals, const vector<double> &vals,
                  const vector<double> &grid_vol_sqs_start, const vector<double> &grid_vol_sqs_end,
                  double vol_sq_fact, double rd, double rf, double dt, int nu, double du, double theta) {
            // if we're in a situation where du_<=0, it's because a knockout's been added that's way in the money
            // and basically certain to knock out. In that case put zeros in for all the values.

            if (du <= 0) {
                for (int i = 0; i < nu; i++)
                    new_vals[i] = 0;
                return;
            }

            // helper function that does the CN backward induction for values in a particular layer
            // (also used for CN forward induction).

            // first calculate the vector that matrix*new values equals - function of the old values

            vector<double> bs(nu);

            double volsq;

            bs[0] = vals[0];
            bs[nu - 1] = vals[nu - 1];
            for (int i = 1; i < nu - 1; i++) {
                volsq = grid_vol_sqs_end[i] * vol_sq_fact;
                bs[i] = (1 - theta) * (volsq / 2. / du / du - (rd - rf - volsq / 2.) / 2. / du) * vals[i - 1]
                        + (1 / dt - (1 - theta) * volsq / du / du) * vals[i]
                        + (1 - theta) * (volsq / 2. / du / du + (rd - rf - volsq / 2.) / 2. / du) * vals[i + 1];
            }

            // do LU-decomposition to solve the problem. The traditional alphas and betas (that
            // define the components of U and L respectively, eg as per Numerical Recipes) are simple
            // to calculate in this tridiagonal case, and are zero if the index is more than one
            // step away from the diagonal. So we don't bother with a matrix. And no need for
            // clever pivoting rules because we know the diagonals are always sensible.

            vector<double> fm(nu), f0(nu), fp(nu);

            // fm is the matrix value at [i,i-1]; f0 is at [i,i], and fp is at [i,i+1].

            fm[0] = 0; // never referenced
            f0[0] = 1;
            fp[0] = 0;
            fm[nu - 1] = 0;
            f0[nu - 1] = 1;
            fp[nu - 1] = 0; // never referenced

            for (int i = 1; i < nu - 1; i++) {
                volsq = grid_vol_sqs_start[i] * vol_sq_fact;
                fm[i] = theta *
                        (-volsq / 2. / du / du + (rd - rf - volsq / 2.) / 2. / du); // at matrix position (i,i-1)
                f0[i] = 1 / dt + theta * volsq / du / du;                         // at matrix position (i,i)
                fp[i] = -theta *
                        (volsq / 2. / du / du + (rd - rf - volsq / 2.) / 2. / du); // at matrix position (i,i+1)
            }

            // betas is beta[i,i]; betasm is beta[i-1,i]; and alphasp is alpha[i+1,i]. All the
            // other components of beta and alpha are zero (except for alpha[i,i] which is always
            // one).

            vector<double> betas(nu), betasm(nu), alphasp(nu);
            for (int i = 0; i < nu; i++) {
                betas[i] = f0[i];
                if (i == 0)
                    betasm[i] = 0; // never referenced
                else {
                    betasm[i] = fp[i - 1];
                    betas[i] -= alphasp[i - 1] * betasm[i];
                }
                if (i == nu - 1)
                    alphasp[i] = 0; // never referenced
                else
                    alphasp[i] = fm[i + 1] / betas[i];
            }

            // solve the linear system in two steps: L*(U*x)=b, so solve
            // L*y=b and get y; then solve U*x=y to get x, what we're really looking for.

            // first put the solution for y into new_vals - this is just temporary so that
            // we don't need to reserve memory for a separate vector.

            new_vals[0] = bs[0];
            for (int i = 1; i < nu; i++)
                new_vals[i] = bs[i] - alphasp[i - 1] * new_vals[i - 1];

            // then replace it with the solution for x - these are the prices in the layer
            // on the time step backward (aside from discounting, which we apply later).

            new_vals[nu - 1] = new_vals[nu - 1] / betas[nu - 1];
            for (int i = nu - 2; i >= 0; i--)
                new_vals[i] = 1 / betas[i] * (new_vals[i] - betasm[i + 1] * new_vals[i + 1]);
        }

        void LVSBICN::go_back(double new_time) {
            if (time_ < 0) THROW(utils::Error, "Cannot backward induct until an initial payout has been added");
            if (new_time > time_)
                THROW(utils::Error, "Cannot backward induct to future time=" << new_time << ", time=" << time_);
            if (new_time == time_) return; // nothing to do

            vector<double> xs(5);
            xs[0] = -1.28;
            xs[1] = -0.68;
            xs[2] = 0;
            xs[3] = 0.68;
            xs[4] = 1.28;
            vector<double> vols(5);

            double ee = exp(eps_);
            double stepdt, next_time, voldd = 0, vold = 0, vol0 = 0, volu = 0, voluu = 0, rd = 0, rf = 0, vol, sqrtT, disc;

            FlatVolSpline *sp = 0;

            vector<double> grid_vol_sqs_start(nu_), grid_vol_sqs_end(nu_), last_grid_vol_sqs(nu_);
            vector<double> new_vals(nu_), new_vals_up(nu_), new_vals_dn(nu_);

            double P0p, P0m, Pp0, Pm0;
            bool first_time = true;

            // track the knockout break times - we may need to resize the grid at
            // those points.

            vector<double> ko_break_times(get_ko_break_times(ko_is_ups_, ko_start_times_, ko_end_times_, time_));

            // track the ko break time index that we're working on and the ko break time corresponding
            // to it

            int ko_break_ind = ko_break_times.size() - 1;
            double ko_break_time = -1;
            if (ko_break_ind >= 0) ko_break_time = ko_break_times[ko_break_ind];

            double umax, umin;
            bool old_min_u_is_ko, old_max_u_is_ko;

            double driftT = log(this->forward() / spot_), driftT_next;

            // remember which parameter break time we're at (note these are not the ko break times -
            // they're breaks in the piecewise-constant parameters: vols and rates).

            int last_break_ind = -1, break_ind;
            double next_break_time;

            vector<double>::const_iterator it = lower_bound(break_times_.begin(), break_times_.end(), time_);
            if (it == break_times_.end())
                break_ind = break_times_.size() - 1; // extrapolate flat
            else
                break_ind = it - break_times_.begin();

            if (break_ind == 0)
                next_break_time = -1;
            else
                next_break_time = break_times_[break_ind - 1];

            int next_ind;
            vector<double> spots(nu_); // used in grid resizings
            double grid_spot;

            while (true) {
                // figure out the next time step to move back to. Land on an integer multiple of
                // dt_ unless that takes us_ past new_time or the latest knockout break time, or
                // the next parameter break time.

                next_ind = int(time_ / dt_);
                next_time = next_ind * dt_;
                if (fabs(next_time - time_) < 1e-12) next_time -= dt_;
                if (next_time < new_time) next_time = new_time;
                if (next_time < ko_break_time) next_time = ko_break_time;
                if (next_time < next_break_time and time_ > next_break_time) next_time = next_break_time;

                stepdt = time_ - next_time; // not always this->dt_
                if (stepdt == 0)
                    THROW(utils::Error, "Ended up with a zero time step, time=" << time_ << "; will never proceed!");

                // figure out the parameters to use for this interval - use the ones at the
                // end of the time interval

                if (break_ind != 0 and time_ <= next_break_time)
                    break_ind -= 1;

                if (break_ind != last_break_ind) {
                    if (break_ind == 0)
                        next_break_time = -1;
                    else
                        next_break_time = break_times_[break_ind - 1];

                    voldd = volsdd_[break_ind];
                    vold = volsd_[break_ind];
                    vol0 = vols0_[break_ind];
                    volu = volsu_[break_ind];
                    voluu = volsuu_[break_ind];

                    if (sp != 0) delete sp; // remove any old spline

                    if (fabs(voldd - vold) < 1e-12 and fabs(voldd - vol0) < 1e-12 and fabs(voldd - volu) < 1e-12 and
                        fabs(voldd - voluu) < 1e-12)
                        sp = 0; // zero out the pointer as a signal not to use a spline; local vol is constant at vol0
                    else {
                        vols[0] = voldd;
                        vols[1] = vold;
                        vols[2] = vol0;
                        vols[3] = volu;
                        vols[4] = voluu;
                        sp = new FlatVolSpline(xs, vols, extrap_fact_);
                    }

                    rd = rds_[break_ind];
                    rf = rfs_[break_ind];

                    last_break_ind = break_ind;
                }

                // get local vols for all grid points, for the interval start and end times.

                if (first_time) {
                    if (sp != 0) sqrtT = sqrt(time_);
                    for (int i = 0; i < nu_; i++) {
                        if (sp == 0)
                            grid_vol_sqs_end[i] = vol0 * vol0;
                        else {
                            // x-value for spline = log(spot/fwd)/vol0/sqrtT. Note that u=log(spot/spot0),
                            // not spot/fwd, so we need to correct it by the drift to this time.
                            vol = sp->value((us_[i] - driftT) / vol0 / sqrtT);
                            grid_vol_sqs_end[i] = vol * vol;
                        }
                    }
                } else
                    grid_vol_sqs_end = last_grid_vol_sqs;

                driftT_next = driftT - (rd - rf) * stepdt;
                if (sp != 0) sqrtT = sqrt(next_time);
                for (int i = 0; i < nu_; i++) {
                    if (sp == 0 or next_time <= 0)
                        grid_vol_sqs_start[i] = vol0 * vol0;
                    else {
                        vol = sp->value((us_[i] - driftT_next) / vol0 / sqrtT);
                        grid_vol_sqs_start[i] = vol * vol;
                    }
                }

                last_grid_vol_sqs = grid_vol_sqs_start; // remember for the next step
                first_time = false; // so it re-uses this calculation in the next step

                // do the backward induction - tridiagonal linear system to solve. Do
                // this for each layer separately. Don't bother doing the high and low layers
                // if the vol of vol is zero - in this case it's pure local vol.

                cnbi(new_vals, vals_, grid_vol_sqs_start, grid_vol_sqs_end, 1., rd, rf, stepdt, nu_, du_, theta_);
                if (alpha_ > 0) {
                    cnbi(new_vals_up, vals_up_, grid_vol_sqs_start, grid_vol_sqs_end, ee, rd, rf, stepdt, nu_, du_,
                         theta_);
                    cnbi(new_vals_dn, vals_dn_, grid_vol_sqs_start, grid_vol_sqs_end, 1. / ee, rd, rf, stepdt, nu_, du_,
                         theta_);
                }

                // extrapolate linearly for edge points

                extrap_vals(new_vals, nu_, min_u_is_ko_, max_u_is_ko_, min_u_rebate_, max_u_rebate_);
                if (alpha_ > 0) {
                    extrap_vals(new_vals_up, nu_, min_u_is_ko_, max_u_is_ko_, min_u_rebate_, max_u_rebate_);
                    extrap_vals(new_vals_dn, nu_, min_u_is_ko_, max_u_is_ko_, min_u_rebate_, max_u_rebate_);
                }

                // mix the layers and apply discounting

                disc = exp(-rd * stepdt);
                P0p = p0p_ * stepdt;
                P0m = p0m_ * stepdt;
                Pp0 = pp0_ * stepdt;
                Pm0 = pm0_ * stepdt;
                for (int i = 0; i < nu_; i++) {
                    if (alpha_ > 0) {
                        if (i == 0 and min_u_is_ko_)
                            vals_[i] = vals_up_[i] = vals_dn_[i] = min_u_rebate_;
                        else if (i == nu_ - 1 and max_u_is_ko_)
                            vals_[i] = vals_up_[i] = vals_dn_[i] = max_u_rebate_;
                        else {
                            vals_[i] = disc *
                                       (new_vals[i] * (1 - P0p - P0m) + new_vals_up[i] * P0p + new_vals_dn[i] * P0m);
                            vals_up_[i] = disc * (new_vals_up[i] * (1 - Pp0) + new_vals[i] * Pp0);
                            vals_dn_[i] = disc * (new_vals_dn[i] * (1 - Pm0) + new_vals[i] * Pm0);
                        }
                    } else if (i == 0 and min_u_is_ko_)
                        vals_[i] = min_u_rebate_;
                    else if (i == nu_ - 1 and max_u_is_ko_)
                        vals_[i] = max_u_rebate_;
                    else
                        vals_[i] = disc * new_vals[i]; // just one layer
                }

                // move the grid time back and stop if we're at the end

                time_ = next_time;
                if (fabs(time_ - new_time) < 1e-12) break;

                // adjust the accumulated risk neutral drift * time to the new time

                driftT = driftT_next;

                // figure out whether to resize the grid. Definitely do it if we're at a knockout
                // break time; otherwise check whether the grid is too narrow based on the forward
                // to this date and the std dev.

                bool resize_grid = false;

                // are we at a knockout break? That requires a resize. Also note the
                // subsequent ko break time to track.

                if (fabs(time_ - ko_break_time) < 1e-12) {
                    ko_break_ind -= 1;
                    if (ko_break_ind == -1)
                        ko_break_time = -1;
                    else
                        ko_break_time = ko_break_times[ko_break_ind];
                    resize_grid = true;
                }

                if (resize_grid) {
                    // construct a spline off the current us_ and grid values for all
                    // three layers. Note these are natural splines with 2nd deriv going to
                    // zero at the edges - not the same as the "flatVolSpline" spline we use
                    // for local vols. We extrapolate linear in spot, though, not in log(spot).

                    for (int i = 0; i < nu_; i++)
                        spots[i] = spot_ * exp(us_[i]);

                    NaturalSpline sUp(spots, vals_up_);
                    NaturalSpline sDn(spots, vals_dn_);
                    NaturalSpline s(spots, vals_);

                    // remember some grid properties so we can properly handle knockouts

                    old_min_u_is_ko = min_u_is_ko_;
                    old_max_u_is_ko = max_u_is_ko_;
                    umin = us_[0];
                    umax = us_[nu_ - 1];

                    // setup the grid again - this properly resizes for new knockouts etc.
                    // Don't reset dt_ though.

                    setup_grid(false);

                    // change the grid values to the interpolated values. If we're stepping
                    // back before a knockout level, explicitly set values to zero to avoid
                    // extrapolation noise.

                    for (int i = 0; i < nu_; i++) {
                        if ((old_min_u_is_ko and us_[i] <= umin) or (old_max_u_is_ko and us_[i] >= umax))
                            vals_[i] = vals_up_[i] = vals_dn_[i] = 0;
                        else {
                            grid_spot = spot_ * exp(us_[i]);
                            vals_[i] = s.value(grid_spot);
                            if (alpha_ > 0) {
                                vals_up_[i] = sUp.value(grid_spot);
                                vals_dn_[i] = sDn.value(grid_spot);
                            }
                        }
                    }

                    // always respect the knockout values if appropriate

                    if (min_u_is_ko_) vals_[0] = vals_up_[0] = vals_dn_[0] = 0;
                    if (max_u_is_ko_) vals_[nu_ - 1] = vals_up_[nu_ - 1] = vals_dn_[nu_ - 1] = 0;
                }
            }

            // if we need to delete the memory for the spline, do so now

            if (sp != 0) delete sp;

            // generated the mixed layer if it's not pure local vol

            if (alpha_ > 0) {
                for (int i = 0; i < nu_; i++)
                    vals_mixed_[i] = vals_[i] * (1 - prob_up_ - prob_dn_) + vals_up_[i] * prob_up_ +
                                     vals_dn_[i] * prob_dn_;
            }
        }

        BaseFI::~BaseFI() {};

        LVSFICN::LVSFICN(double spot,
                         const vector<double> &volsdd, const vector<double> &volsd,
                         const vector<double> &vols0,
                         const vector<double> &volsu, const vector<double> &volsuu,
                         const vector<double> &rds,
                         const vector<double> &rfs,
                         const vector<double> &break_times,
                         double alpha, double beta, double extrap_fact,
                         int nu, int nt, double nsd, double time_max,
                         double theta)
                : spot_(spot), volsdd_(volsdd), volsd_(volsd), vols0_(vols0), volsu_(volsu), volsuu_(volsuu),
                  rds_(rds), rfs_(rfs), break_times_(break_times), alpha_(alpha), beta_(beta),
                  extrap_fact_(extrap_fact),
                  nu_(nu), nt_(nt), nsd_(nsd), time_max_(time_max), theta_(theta) {
            // figure out the layer size and transition probability frequencies

            double gamma = alpha / 2. / sqrt(beta);
            eps_ = 2 * asinh(gamma);
            double ee = exp(eps_);

            p0p_ = beta * (0.5 - gamma / 2. / sqrt(1 + gamma * gamma));
            p0m_ = beta * (0.5 + gamma / 2. / sqrt(1 + gamma * gamma));
            pp0_ = beta;
            pm0_ = beta;

            // initialize the grid. For a vol level, use the maximum of all the local vol
            // parameters in the high vol state.

            double max_vol = 0.;
            for (size_t i = 0; i < volsdd.size(); i++) {
                if (volsdd[i] > max_vol) max_vol = volsdd[i];
                if (volsd[i] > max_vol) max_vol = volsd[i];
                if (vols0[i] > max_vol) max_vol = vols0[i];
                if (volsu[i] > max_vol) max_vol = volsu[i];
                if (volsuu[i] > max_vol) max_vol = volsuu[i];
            }

            max_vol *= sqrt(ee); // want high vol state vol

            double umin = -nsd * max_vol * sqrt(time_max);
            double umax = nsd * max_vol * sqrt(time_max);
            du_ = (umax - umin) / (nu - 1);

            dt_ = time_max / nt;

            us_.resize(nu);
            vals_.resize(nu);
            vals_up_.resize(nu);
            vals_dn_.resize(nu);
            vals_mixed_.resize(nu); // values added later after forward induction
            double u, val;
            for (int i = 0; i < nu; i++) {
                u = umin + i * du_;
                us_[i] = u;
                val = 1 - exp(u);
                if (val < 0) val = 0;
                vals_[i] = vals_up_[i] = vals_dn_[i] = val;
            }

            // note the current grid time to expiration is zero

            time_ = 0;

            // add a placer for the total discount factor - it'll update as we step forward

            tot_disc_ = 1.;

            // calculate the probabilities of being in the different vol layers. We assume
            // we're in a stationary distribution at the start (ie the initial vol layer
            // is uncertain too). We only track the probabilities of being in the up and down
            // states because the probability of being in the mid state is 1-those two.

            prob_up_ = 0.25 - gamma / 4. / sqrt(1 + gamma * gamma);
            prob_dn_ = 0.25 + gamma / 4. / sqrt(1 + gamma * gamma);
        }

        LVSFICN::~LVSFICN() {}

        double LVSFICN::forward() const {
            double fwd = spot_;
            double last_time, end_time, rate_time;
            int const n = break_times_.size();

            for (int i = 0; i < n; i++) {
                if (i == 0)
                    last_time = 0;
                else
                    last_time = break_times_[i - 1];
                rate_time = break_times_[i];
                if (rate_time < time_)
                    end_time = rate_time;
                else
                    end_time = time_;

                fwd *= exp((rds_[i] - rfs_[i]) * (end_time - last_time));
                if (end_time > time_ - 1e-12) break;
            }

            if (time_ > break_times_[n - 1])
                fwd *= exp((rds_[n - 1] - rfs_[n - 1]) * (time_ - break_times_[n - 1]));

            return fwd;
        }

        void LVSFICN::go_forward(double new_time) {
            if (new_time < time_)
                THROW(utils::ValueError, "Cannot forward induct to past time=" << new_time << ", time=" << time_);
            if (new_time > time_max_)
                THROW(utils::ValueError, "Cannot forward induct to " << new_time << " which is past the maximum time="
                                                                     << time_max_);
            if (new_time == time_) return; // nothing to do

            vector<double> xs(5);
            xs[0] = -1.28;
            xs[1] = -0.68;
            xs[2] = 0;
            xs[3] = 0.68;
            xs[4] = 1.28;
            vector<double> vols(5);

            double ee = exp(eps_);
            double stepdt, next_time, voldd = 0, vold = 0, vol0 = 0, volu = 0, voluu = 0, rd = 0, vol, sqrtT, disc;

            FlatVolSpline *sp = 0;

            vector<double> grid_vol_sqs_start(nu_), grid_vol_sqs_end(nu_), last_grid_vol_sqs(nu_);
            vector<double> new_vals(nu_), new_vals_up(nu_), new_vals_dn(nu_);

            double P0p, P0m, Pp0, Pm0;
            bool first_time = true;

            // figure out where we are in the piecewise constant parameters

            int np = break_times_.size();
            int last_break_ind = -1, break_ind;
            double break_time;
            vector<double>::const_iterator it = lower_bound(break_times_.begin(), break_times_.end(), time_);
            if (it == break_times_.end())
                break_ind = np - 1; // extrapolate flat
            else
                break_ind = it - break_times_.begin();
            break_time = break_times_[break_ind];

            int next_ind;

            while (true) {
                // figure out the next time step to move back to. Land on an integer multiple of
                // dt_ unless that takes us_ past new_time or the next break time.

                next_ind = int(time_ / dt_) + 1;
                next_time = next_ind * dt_;
                if (next_time <= time_) next_time = time_ + dt_;
                if (next_time > new_time) next_time = new_time;
                if (next_time > break_time and break_time > time_) next_time = break_time;

                stepdt = next_time - time_; // not always this->dt_
                if (stepdt == 0)
                    THROW(utils::ValueError, "Step dt_ equal to zero at time=" << time_
                                                                               << "; calculation will never proceed!");

                // figure out the parameters to use for this interval - use the ones at the
                // end of the time interval

                if (break_ind != np - 1 and next_time >= break_times_[break_ind])
                    break_ind++;

                if (break_ind != last_break_ind) {
                    break_time = break_times_[break_ind];

                    voldd = volsdd_[break_ind];
                    vold = volsd_[break_ind];
                    vol0 = vols0_[break_ind];
                    volu = volsu_[break_ind];
                    voluu = volsuu_[break_ind];

                    if (sp != 0) delete sp; // remove any old spline

                    if (fabs(voldd - vold) < 1e-12 and fabs(voldd - vol0) < 1e-12 and fabs(voldd - volu) < 1e-12 and
                        fabs(voldd - voluu) < 1e-12)
                        sp = 0; // zero out the pointer as a signal not to use a spline
                    else {
                        vols[0] = voldd;
                        vols[1] = vold;
                        vols[2] = vol0;
                        vols[3] = volu;
                        vols[4] = voluu;
                        sp = new FlatVolSpline(xs, vols, extrap_fact_);
                    }

                    rd = rds_[break_ind];

                    last_break_ind = break_ind;
                }

                // get local vols for all grid points, for the interval start and end times.

                if (first_time) {
                    if (sp != 0) sqrtT = sqrt(time_);
                    for (int i = 0; i < nu_; i++) {
                        if (sp == 0 or time_ <= 0)
                            grid_vol_sqs_start[i] = vol0 * vol0;
                        else {
                            // x-value for spline = u/vol0/sqrtT.
                            vol = sp->value(us_[i] / vol0 / sqrtT);
                            grid_vol_sqs_start[i] = vol * vol;
                        }
                    }
                } else
                    grid_vol_sqs_start = last_grid_vol_sqs;

                if (sp != 0) sqrtT = sqrt(next_time);
                for (int i = 0; i < nu_; i++) {
                    if (sp == 0)
                        grid_vol_sqs_end[i] = vol0 * vol0;
                    else {
                        vol = sp->value(us_[i] / vol0 / sqrtT);
                        grid_vol_sqs_end[i] = vol * vol;
                    }
                }

                last_grid_vol_sqs = grid_vol_sqs_end; // remember for the next step
                first_time = false; // so it re-uses this calculation in the next step

                // do the backward induction - tridiagonal linear system to solve. Do
                // this for each layer separately.

                cnbi(new_vals, vals_, grid_vol_sqs_end, grid_vol_sqs_start, 1., 0, 0, stepdt, nu_, du_, theta_);
                cnbi(new_vals_up, vals_up_, grid_vol_sqs_end, grid_vol_sqs_start, ee, 0, 0, stepdt, nu_, du_, theta_);
                cnbi(new_vals_dn, vals_dn_, grid_vol_sqs_end, grid_vol_sqs_start, 1. / ee, 0, 0, stepdt, nu_, du_,
                     theta_);

                // extrapolate linearly for edge points

                extrap_vals(new_vals, nu_, false, false, 0, 0);
                extrap_vals(new_vals_up, nu_, false, false, 0, 0);
                extrap_vals(new_vals_dn, nu_, false, false, 0, 0);

                // mix the layers

                P0p = p0p_ * stepdt;
                P0m = p0m_ * stepdt;
                Pp0 = pp0_ * stepdt;
                Pm0 = pm0_ * stepdt;
                for (int i = 0; i < nu_; i++) {
                    vals_[i] = new_vals[i] * (1 - P0p - P0m) + new_vals_up[i] * P0p + new_vals_dn[i] * P0m;
                    vals_up_[i] = new_vals_up[i] * (1 - Pp0) + new_vals[i] * Pp0;
                    vals_dn_[i] = new_vals_dn[i] * (1 - Pm0) + new_vals[i] * Pm0;
                }

                // keep track of the total discount

                disc = exp(-rd * stepdt);
                tot_disc_ *= disc;

                // move the grid time forward and stop if we're at the end

                time_ = next_time;
                if (fabs(time_ - new_time) < 1e-12) break;
            }

            // if we need to delete the memory for the spline, do so now

            if (sp != 0) delete sp;

            // construct the mixed values here - this is what we interpolate over to get the
            // actual call option price, averaging over values in the different layers by the
            // probability of ending in the layer.

            for (int i = 0; i < nu_; i++)
                vals_mixed_[i] = vals_[i] * (1 - prob_up_ - prob_dn_) + vals_up_[i] * prob_up_ + vals_dn_[i] * prob_dn_;
        }

        double LVSFICN::interp(double strike) const {
            // if we're past the edges of the grid, return the intrinsic value for
            // the call option (since in this case we know the payoff - no time value).

            double fwd = this->forward();
            double u = log(strike / fwd);

            if (u < us_[0]) return (fwd - strike) * tot_disc_;
            if (u > us_[nu_ - 1]) return 0;

            // construct the spline over the call price vs strike values;
            // though strike is in log(strike/fwd) units. Weight by probabilities
            // of ending in different states.

            NaturalSpline s(us_, vals_mixed_);

            // interpolate the adjusted call price, which is the real call
            // price divided by the discount factor and the forward. Convert
            // it back to a real call price when returning it.

            return s.value(u) * tot_disc_ * fwd;
        }

    } // namespace analytics
} // namespace wst
