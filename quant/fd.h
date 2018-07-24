/*
 * Haocheng Wang, July 09, 2018
 * Last modification:
 * 2018-07-10: Add payoff smoothing for LVSBI
 *
 * Finite difference grid classes, for backward and forward induction to price derivatives
 * under various models.
*/


#ifndef _IN_FD_H
#define _IN_FD_H

#include <vector>
#include "payoff.h"

namespace wst {
    namespace analytics {

        class BaseBI {
            // base class for backward inductors

        public:
            BaseBI() : time_(-1) {}; // time<0 a flag that the grid isn't initialized
            virtual ~BaseBI();

            // go_back moves the grid back from the current grid time to the selected newTime.

            virtual void go_back(double new_time)=0;

            // interp interpolates a value off the grid.

            virtual double interp(double interp_spot) const =0;

            // initialize_payoff puts on the first payoff, at the final grid time.

            virtual void initialize_payoff(const Payoff &p, double time)=0;

            // add_payoff adds a new payoff onto whatever values are currently on the grid

            virtual void add_payoff(const Payoff &p)=0;

            // multiply_payoff scales the current grid values by the value of the payoff.

            virtual void multiply_payoff(const Payoff &p)=0;

            // max_payoff sets all grid values to max(payoff value,grid value)

            virtual void max_payoff(const Payoff &p)=0;

            // add_knockout adds a new knockout level to the grid - must be added before the
            // grid is initialized with a payoff. start_time<0 means the start time is the current
            // pricing time. end_time<0 means the end time is always farther in the future than
            // any payoff. rebate!=0 means that on knock a cash rebate is paid.

            virtual void
            add_knockout(double ko, bool is_up, double start_time = -1, double end_time = -1, double rebate = 0);

        protected:
            // time represents the current time of the grid. <0 is a flag that the grid has not
            // been initialized yet - we assume this is done in the constructors for the grid
            // classes.

            double time_;

            // knockout-related data

            std::vector<double> ko_levels_;
            std::vector<bool> ko_is_ups_;
            std::vector<double> ko_start_times_;
            std::vector<double> ko_end_times_;
            std::vector<double> ko_rebates_;
        };

        class LVSBIBase : public BaseBI {
            // LV-state model base class. The LV-state model is a local vol/stochastic vol mixture
            // model that represents the stochastic vol as three separate local vol states. The model
            // spec is
            //     dS(t)/S(t) = mu dt + sigma(S,t) sqrt(v(t)) dz(t)
            //     v(t) = exp(Y(t) eps)
            // Y(t) is a variable that can be +1, -1, or 0, and represents the three stoch vol
            // states. sigma(S,t) is a local vol function, which is parameterized as a cubic
            // spline in x=ln(S(t)/F(0,t))/vol0/sqrt(t), where F(0,t) is the initial forward to
            // t and vol0 is a marked parameter of the cubic spline. The cubic spline is defined
            // by five marked vol parameters: voldd, vold, vol0, volu, and voluu. Those correspond
            // to five specific x values: -1.28, -0.68, 0, +0.68, and +1.28. Those correspond
            // roughly to ATM, 25d, and 10d levels.
            //   The transition probabilities between the three states are determined by the
            // stoch vol parameters alpha and beta, meant to be similar to the Heston vol of
            // vol and mean reversion parameters: the layer transition probabilities in time dt
            // are set to match the mean and variance of Heston vol changes. The layer separation
            // variable eps is defined by those as well.
            //   The risk-neutral drift and local vol parameters are assumed to be piecewise-
            // constant in time, with common break times.
            //   Note: we assume that the initial distribution of Y is its stationary distribution;
            // it does not start with Y=0. This is like assuming the Heston v is drawn from its
            // stationary distribution at t=0 instead of starting with a known value.

        public:
            // initializer for grid takes:
            //   spot: current spot price of asset
            //   volsdd, volsd, vols0, volsu, volsuu: piecewise-constant local vols for reference
            //                                        Xs (-1.28,-0.68,0,+0.68,+1.28)
            //   rds, rfs: piecewise-constant instantaneous interest rates for denominated and asset
            //   break_times: break times for piecewise-constant vols & rates
            //   alpha: vol of vol parameter - same units as Heston
            //   beta: vol mean reversion parameter - same units as Heston
            //   extrap_fact: extrapolation factor for custom local vol cubic spline - vols flatten
            //                out this many std devs from either side of the marked points
            //   nu: number of log(spot) points
            //   nt: number of time points (to time at first payoff initialization)
            //   nsd: number of std devs to include in the grid

            LVSBIBase(double spot,
                      const std::vector<double> &volsdd, const std::vector<double> &volsd,
                      const std::vector<double> &vols0,
                      const std::vector<double> &volsu, const std::vector<double> &volsuu,
                      const std::vector<double> &rds,
                      const std::vector<double> &rfs,
                      const std::vector<double> &break_times,
                      double alpha, double beta, double extrap_fact,
                      int nu, int nt, double nsd);

            virtual ~LVSBIBase();

            virtual double interp(double interp_spot) const;

            virtual void initialize_payoff(const Payoff &p, double time);

            virtual void add_payoff(const Payoff &p);

            virtual void multiply_payoff(const Payoff &p);

            virtual void max_payoff(const Payoff &p);

            double forward() const;

        protected:
            double spot_;
            std::vector<double> volsdd_, volsd_, vols0_, volsu_, volsuu_;
            std::vector<double> rds_, rfs_;
            std::vector<double> break_times_;
            double alpha_, beta_, extrap_fact_;
            int nu_, nt_;
            double nsd_;
            double du_, dt_;

            std::vector<double> us_;
            std::vector<double> vals_, vals_up_, vals_dn_, vals_mixed_;

            double eps_;
            double p0p_, p0m_, pp0_, pm0_;

            bool min_u_is_ko_, max_u_is_ko_;
            double min_u_rebate_, max_u_rebate_;

            double max_vol_; // used for grid sizing

            double prob_up_, prob_dn_;

            void setup_grid(bool recalcdt);
        };

        class LVSBIExp : public LVSBIBase {
            // LV-state model backward inductor - explicit backward induction

        public:
            LVSBIExp(double spot,
                     const std::vector<double> &volsdd, const std::vector<double> &volsd,
                     const std::vector<double> &vols0,
                     const std::vector<double> &volsu, const std::vector<double> &volsuu,
                     const std::vector<double> &rds,
                     const std::vector<double> &rfs,
                     const std::vector<double> &break_times,
                     double alpha, double beta, double extrap_fact,
                     int nu, int nt, double nsd);

            virtual ~LVSBIExp();

            virtual void go_back(double new_time);
        };

        class LVSBICN : public LVSBIBase {
            // LV-state model backward inductor - Crank-Nicholson backward induction.
            // theta controls the induction type: theta=0 is explicit; theta=1 is fully
            // implicit; and theta=0.5 is traditional CN.

        public:
            LVSBICN(double spot,
                    const std::vector<double> &volsdd, const std::vector<double> &volsd,
                    const std::vector<double> &vols0,
                    const std::vector<double> &volsu, const std::vector<double> &volsuu,
                    const std::vector<double> &rds,
                    const std::vector<double> &rfs,
                    const std::vector<double> &break_times,
                    double alpha, double beta, double extrap_fact,
                    int nu, int nt, double nsd,
                    double theta = 0.5);

            virtual ~LVSBICN();

            virtual void go_back(double new_time);

        protected:
            double theta_;
        };

        class BaseFI {
            // base class for forward inductors: calculates call option price as a function
            // of strike and time to expiration.

        public:
            virtual ~BaseFI();

            virtual void go_forward(double new_time)=0;

            virtual double interp(double strike) const =0;
        };

        class LVSFICN : public BaseFI {
            // LV-state model, as with LVSBIBase.
            // Crank-Nicholson forward inductor for call option price as a function of strike
            // and time to expiration.

        public:
            LVSFICN(double spot,
                    const std::vector<double> &volsdd, const std::vector<double> &volsd,
                    const std::vector<double> &vols0,
                    const std::vector<double> &volsu, const std::vector<double> &volsuu,
                    const std::vector<double> &rds,
                    const std::vector<double> &rfs,
                    const std::vector<double> &break_times,
                    double alpha, double beta, double extrap_fact,
                    int nu, int nt, double nsd, double time_max,
                    double theta = 0.5);

            virtual ~LVSFICN();

            virtual void go_forward(double new_time);

            virtual double interp(double strike) const;

            double forward() const;

        protected:
            double spot_;
            std::vector<double> volsdd_, volsd_, vols0_, volsu_, volsuu_;
            std::vector<double> rds_, rfs_;
            std::vector<double> break_times_;
            double alpha_, beta_, extrap_fact_;
            int nu_, nt_;
            double nsd_;
            double time_max_, time_;
            double du_, dt_;
            double theta_;
            double tot_disc_;

            std::vector<double> us_;
            std::vector<double> vals_, vals_up_, vals_dn_, vals_mixed_;

            double prob_up_, prob_dn_;

            double eps_;
            double p0p_, p0m_, pp0_, pm0_;
        };

    } // namespace analytics
} // namespace wst

#endif
