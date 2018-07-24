/*
 * Haocheng Wang, July 09, 2018
 * Last modification:
 * 2018-07-10: Add option for smoothing the digital payoff
 * Payoff functions - used in various pricing algorithms that take generic boundary
 * conditions.
*/

#ifndef _IN_PAYOFF_H
#define _IN_PAYOFF_H

namespace wst {
    namespace analytics {

        class Payoff {
            // Base payoff class - supports returning a value for a given asset spot price

        public:
            virtual double value(double spot) const =0;

            virtual ~Payoff();
        };

        class ZeroPayoff : public Payoff {
            // Payoff always returning zero

        public:
            double value(double spot) const override { return 0; };

            ~ZeroPayoff() override;
        };

        class UnitPayoff : public Payoff {
            // Payoff that returns a constant value

        public:
            explicit UnitPayoff(double quant = 1.) : quant_(quant) {};

            ~UnitPayoff() override;

            double value(double spot) const override { return quant_; };

        private:
            double quant_;
        };

        class OptionPayoff : public Payoff {
            // Payoff that returns a call or put payoff

        public:
            OptionPayoff(bool is_call, double strike, double quant = 1.);

            ~OptionPayoff() override;

            double value(double spot) const override;

        private:
            bool is_call_;
            double strike_;
            double quant_;
        };

        class DigitalPayoff : public Payoff {
            // Payoff that returns a digital payoff: a quantity (defaulting to 1) or
            // zero depending on whether the spot price is above or below the strike.
            // Value at the strike is quant for a call, 0 for a put (so that digital
            // call + digital put always equals quant, even at the strike).

        public:
            DigitalPayoff(bool is_call, double strike, double quant = 1., bool smoothing = false,
                          double epsilon = 0.02);

            ~DigitalPayoff() override;

            double value(double spot) const override;

        private:
            bool is_call_;
            double strike_;
            double quant_;
            bool smoothing_;
            double epsilon_;
        };

    } // namespace analytics
} // namespace wst

#endif
