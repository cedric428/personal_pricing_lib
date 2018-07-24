/*
 * Haocheng Wang, July 09, 2018
 * Cubic spline classes, different versions for different boundary conditions.
*/

#ifndef _IN_SPLINE_H
#define _IN_SPLINE_H

#include <vector>

namespace wst
{
namespace analytics
{

class Spline
{
    // base class for all splines
public:
    Spline() {};
    virtual ~Spline();
    
    virtual double value(double x) const=0;
};

class NaturalSpline : public Spline
{
    // Cubic spline with natural boundary conditions: ie assume that second
    // deriv goes to zero at edges.
    
public:
    NaturalSpline(const std::vector<double>& xs,const std::vector<double>& ys);
    virtual ~NaturalSpline();
    
    virtual double value(double x) const;
    
private:
    std::vector<double> xs_;
    std::vector<double> ys_;
    std::vector<double> ypps_; // second derivs
    int n_;
};

class FlatVolSpline : public Spline
{
    // Cubic spline for volatility interpolation. 
    // Extrapolates flat after extrap_fact std devs on either side of the edges. 
    // Assumes that x-values are log(spot/fwd)/vol0/sqrt(T),
    // where vol0 = middle vol level, for x=0 value. So x=0 needs to be in xs.

public:
    FlatVolSpline(const std::vector<double>& xs,const std::vector<double>& vols,double extrap_fact);
    virtual ~FlatVolSpline();
    
    virtual double value(double x) const;

private:
    std::vector<double> xs_;
    std::vector<double> vols_;
    std::vector<double> cubicParams_;
    int n_;
    
    double xl_, xr_;
};

} // namespace analytics
} // namespace wst

#endif