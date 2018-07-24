/*
 * Haocheng Wang, July 09, 2018
 * Last modification:
 * 2018-07-10: Add option for smoothing the digital payoff
*/

#include <cmath>

#include "payoff.h"

namespace wst {
    namespace analytics {

        Payoff::~Payoff() {}

        ZeroPayoff::~ZeroPayoff() {}

        UnitPayoff::~UnitPayoff() {}

        OptionPayoff::OptionPayoff(bool is_call, double strike, double quant)
                : is_call_(is_call), strike_(strike), quant_(quant) {
        }

        OptionPayoff::~OptionPayoff() {}

        double OptionPayoff::value(double spot) const {
            double imp_val = (spot - strike_);
            if (!is_call_) imp_val *= -1;
            if (imp_val > 0)
                return imp_val * quant_;
            else
                return 0;
        }

        DigitalPayoff::DigitalPayoff(bool is_call, double strike, double quant, bool smoothing, double epsilon)
                : is_call_(is_call), strike_(strike), quant_(quant), smoothing_(smoothing), epsilon_(epsilon) {
        }

        DigitalPayoff::~DigitalPayoff() {}

        double DigitalPayoff::value(double spot) const {
            if (smoothing_) {
                if (is_call_) {
                    return 0.5 * std::tanh((spot - strike_) / epsilon_) + 0.5;
                } else {
                    return 0.5 * std::tanh((strike_ - spot) / epsilon_) +
                           0.5;  // this smoothing method is sensitive to epsilon_
                }
            } else {
                if ((is_call_ and spot >= strike_) or (not is_call_ and spot < strike_))
                    return quant_;
                else
                    return 0;
            }


        }

    } // namespace analytics
} // namespace wst
