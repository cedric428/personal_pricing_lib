/*
 * Haocheng Wang, July 09, 2018
*/

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include "spline.h"
#include "utils/exceptions.h"

using namespace std;

namespace wst
{
namespace analytics
{

Spline::~Spline() {}

NaturalSpline::NaturalSpline(const vector<double>& xs,const vector<double>& ys)
  : xs_(xs), ys_(ys)
{
    // initialize
    
    if(xs.size()!=ys.size()) THROW(utils::ValueError, "Number of xs and ys must be equal");
    if(xs.size()<3) THROW(utils::ValueError, "Number of xs must be at least 3");
    
    n_=xs.size();
    ypps_.resize(n_);
    
    // calculate the ypps - the second derivatives of the functions at each point
    
    vector<double> us(n_-1);
    ypps_[0] = us[0] = 0;
    
    int i, k;
    double p,sig;
    
    for (i=1; i<n_-1; i++)
    {
        sig      = (xs[i]-xs[i-1])/(xs[i+1]-xs[i-1]);
        p        = sig*ypps_[i-1]+2.;
        ypps_[i] = (sig-1.)/p;
        us[i]    = (ys[i+1]-ys[i])/(xs[i+1]-xs[i]) - (ys[i]-ys[i-1])/(xs[i]-xs[i-1]);
        us[i]    = (6.*us[i]/(xs[i+1]-xs[i-1])-sig*us[i-1])/p;
    }
    
    ypps_[n_-1] = 0;
    for (k=n_-2;k>=0;k--)
        ypps_[k] = ypps_[k]*ypps_[k+1]+us[k];
}

NaturalSpline::~NaturalSpline() {}

double NaturalSpline::value(double x) const
{
    // figure out where x lies in the set of xs
    
    vector<double>::const_iterator it=lower_bound(xs_.begin(),xs_.end(),x);
    
    // if we're doing extrapolation, use linear extrapolation (since we assume y''=0 at each edge).
    // Use the spline parameters to figure out y' at the edge.
    
    if (it==xs_.begin())
    {
        double yp=(ys_[1]-ys_[0])/(xs_[1]-xs_[0])-ypps_[1]/6.*(xs_[1]-xs_[0]);
        return ys_[0] + yp*(x-xs_[0]);
    }
    
    if (it==xs_.end())
    {
        double yp=(ys_[n_-1]-ys_[n_-2])/(xs_[n_-1]-xs_[n_-2])+ypps_[n_-2]/6.*(xs_[n_-1]-xs_[n_-2]);
        return ys_[n_-1] + yp*(x-xs_[n_-1]);
    }
    
    // otherwise do the cubic spline interpolation
    
    int ind=it-xs_.begin(); // index of the next element
    double const y1=ys_[ind-1];
    double const y2=ys_[ind];
    double const x1=xs_[ind-1];
    double const x2=xs_[ind];
    double const ypp1=ypps_[ind-1];
    double const ypp2=ypps_[ind];
    
    double const A=(x2-x)/(x2-x1);
    double const B=1-A;
    double const C=1./6*(A*A*A-A)*(x2-x1)*(x2-x1);
    double const D=1./6*(B*B*B-B)*(x2-x1)*(x2-x1);
    
    return A*y1+B*y2+C*ypp1+D*ypp2;
}

void linear_system_solve(vector< vector<double> >& A, vector<double>& b)
{
    // use Gauss-Jordan elimination just because it's simple, and the systems are relatively small here anyways.
    // Solution ends up in b; A ends up with its inverse.
    
    int icol, irow, n;
    n = b.size();
    
    double big, dum, pivinv, temp;
    
    vector<int> indxc(n), indxr(n), ipiv(n);
    
    for (int j=0;j<n;j++) ipiv[j] = 0;
    for (int i=0;i<n;i++)
    {
        big = 0;
        for (int j=0;j<n;j++)
            if (ipiv[j]!=1)
            {
                for (int k=0;k<n;k++)
                {
                    if (ipiv[k]==0)
                    {
                        if (fabs(A[j][k])>=big)
                        {
                            big = fabs(A[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        
        ++(ipiv[icol]);
        
        if (irow!=icol)
        {
            for (int l=0;l<n;l++)
            {
                temp = A[irow][l];
                A[irow][l] = A[icol][l];
                A[icol][l] = temp;
            }
            temp    = b[irow];
            b[irow] = b[icol];
            b[icol] = temp;
        }
        indxr[i] = irow;
        indxc[i] = icol;
        
        if (A[icol][icol]==0) THROW(utils::ValueError, "Singular matrix");
        pivinv = 1./A[icol][icol];
        
        A[icol][icol] = 1.;
        for (int l=0;l<n;l++) A[icol][l] *= pivinv;
        b[icol] *= pivinv;
        for (int ll=0;ll<n;ll++)
            if (ll!=icol)
            {
                dum = A[ll][icol];
                A[ll][icol] = 0;
                for (int l=0;l<n;l++) A[ll][l] -= A[icol][l]*dum;
                b[ll] -= b[icol]*dum;
            }
    }
    
    for (int l=n-1;l>=0;l--)
        if (indxr[l]!=indxc[l])
            for (int k=0;k<n;k++)
            {
                temp = A[k][indxr[l]];
                A[k][indxr[l]] = A[k][indxc[l]];
                A[k][indxc[l]] = temp;
            }
}

FlatVolSpline::FlatVolSpline(const vector<double>& xs,const vector<double>& vols,double extrap_fact)
  : xs_(xs), vols_(vols)
{
    // figure out which x-value corresponds to x=0 and get the vol - that determines scale
    // for extrapolating.
    
    n_ = xs.size();
    
    bool found_zero=false;
    int zero_ind;
    for (int i=0;i<n_;i++)
    {
        if (fabs(xs[i])<1e-12)
        {
            found_zero = true;
            zero_ind   = i;
            break;
        }
    }
    if (!found_zero) THROW(utils::ValueError, "Could not find 0 in the list of xs");
    
    double vol0=vols[zero_ind];
    
    xl_ = xs[0] - vols[0]/vol0*extrap_fact;
    xr_ = xs[n_-1] + vols[n_-1]/vol0*extrap_fact;
    
    // set up the matrixes we'll use for solving for the cubic parameters
    
    int m=4*n_+4;
    vector< vector<double> > A(m);
    vector<double> b(m);
    for (int i=0;i<m;i++)
    {
        A[i].resize(m);
        for (int j=0;j<m;j++) A[i][j] = 0.;
        b[i] = 0;
    }
    
    // handle the edge conditions first
    
    A[m-4][1] = 1;
    A[m-4][2] = 2*xl_;
    A[m-4][3] = 3*xl_*xl_;
    
    A[m-3][2] = 2;
    A[m-3][3] = 6*xl_;
    
    A[m-2][m-3] = 1;
    A[m-2][m-2] = 2*xr_;
    A[m-2][m-1] = 3*xr_*xr_;
    
    A[m-1][m-2] = 3;
    A[m-1][m-1] = 6*xr_;
    
    // then the internals
    
    double x;
    for (int i=0;i<n_;i++)
    {
        x = xs[i];
        
        // values are equal at the xs[i] point for region i-1 and i
        
        A[3*i][4*i]   = 1;
        A[3*i][4*i+4] = -1;
        A[3*i][4*i+1] = x;
        A[3*i][4*i+5] = -x;
        A[3*i][4*i+2] = x*x;
        A[3*i][4*i+6] = -x*x;
        A[3*i][4*i+3] = x*x*x;
        A[3*i][4*i+7] = -x*x*x;
        
        // slopes are equal there too
        
        A[3*i+1][4*i+1] = 1;
        A[3*i+1][4*i+5] = -1;
        A[3*i+1][4*i+2] = 2*x;
        A[3*i+1][4*i+6] = -2*x;
        A[3*i+1][4*i+3] = 3*x*x;
        A[3*i+1][4*i+7] = -3*x*x;
        
        // and 2nd derivs too
        
        A[3*i+2][4*i+2] = 2;
        A[3*i+2][4*i+6] = -2;
        A[3*i+2][4*i+3] = 6*x;
        A[3*i+2][4*i+7] = -6*x;
    }
    
    // then the actual values of the functions at each point
    
    for (int i=0;i<n_;i++)
    {
        x = xs[i];
        A[3*n_+i][4*i] = 1;
        A[3*n_+i][4*i+1] = x;
        A[3*n_+i][4*i+2] = x*x;
        A[3*n_+i][4*i+3] = x*x*x;
        
        b[3*n_+i] = vols_[i];
    }
    
    // then solve the linear system - this gives us a, b, c, and d such that
    // vol = a + b*x + c*x^2 + d*x^3 in each region.
    
    linear_system_solve(A,b);
    cubicParams_ = b;
}

FlatVolSpline::~FlatVolSpline() {}

double FlatVolSpline::value(double x) const
{
    if (x<xl_) x = xl_;
    if (x>xr_) x = xr_;
    
    // get the position in the input array
    
    vector<double>::const_iterator it=lower_bound(xs_.begin(),xs_.end(),x);
    int ind=it-xs_.begin();
    
    double const a=cubicParams_[4*ind];
    double const b=cubicParams_[4*ind+1];
    double const c=cubicParams_[4*ind+2];
    double const d=cubicParams_[4*ind+3];
    
    return a+b*x+c*x*x+d*x*x*x;
}

} // namespace analytics
} // namespace wst
