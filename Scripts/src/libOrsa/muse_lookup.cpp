// Adapted from https://github.com/vxl/vxl/tree/c6c899aaf9cbf7ccabfca1ba35b2248fb611ffbc
// Algorithm from MUSE: Robust Surface Fitting using Unbiased Scale Estimates, James V. Miller and Charles V. Stewart
// By Clément Riu
// 2021


#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>

#include "muse_lookup.hpp"

bool operator< ( MuseKeyType const& left_t, MuseKeyType const& right_t )
{
  return left_t.n_ < right_t.n_
    || ( left_t.n_ == right_t.n_ && left_t.k_ < right_t.k_ );
}

double
MuseTable::expected_kth( unsigned int k, unsigned int n )
{
  assert( 0<k && k<= n );
  MuseKeyType key(k,n);
  MuseTableEntry& entry = _table[key];
  if ( ! entry.initialized_ )
    calculate_all( k, n, entry );
  return entry.expected_;
}

double
MuseTable::standard_dev_kth( unsigned int k, unsigned int n )
{
  assert( 0<k && k<=n );
  MuseKeyType key(k,n);
  MuseTableEntry& entry = _table[key];
  if ( ! entry.initialized_ )
    calculate_all( k, n, entry );
  return entry.standard_dev_;
}

double
MuseTable::muset_divisor( unsigned int k, unsigned int n )
{
  assert( 0<k && k<= n );
  MuseKeyType key(k,n);
  MuseTableEntry& entry = _table[key];
  if ( ! entry.initialized_ )
    calculate_all( k, n, entry );
  return entry.muse_t_divisor_;
}


double
MuseTable::muset_sq_divisor( unsigned int k, unsigned int n )
{
  assert( 0<k && k<= n );
  MuseKeyType key(k,n);
  MuseTableEntry& entry = _table[key];
  if ( ! entry.initialized_ )
    calculate_all( k, n, entry );
  return entry.muse_t_sq_divisor_;
}

void
MuseTable::calculate_all( unsigned int k, unsigned int n,
                                MuseTableEntry & entry )
{
  entry.initialized_ = true;
  entry.expected_ = calculate_expected( k, n );
  entry.standard_dev_ = calculate_standard_dev( k, n, entry.expected_ );
  entry.muse_t_divisor_ = calculate_divisor( k, n, entry.expected_ );
  entry.muse_t_sq_divisor_ = calculate_sq_divisor( k, n, entry.expected_ );
}

double ERFCC( double x )
{
    double t,z,ans;
    z=std::fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*std::exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
           t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
                 t*(-0.82215223+t*0.17087277)))))))));
    return  x >= 0.0 ? ans : 2.0-ans;
}


double GaussianCDFInv( double p )
{
   double dp , dx , dt , ddq , dq ;
   int newt ;

   dp = (p <= 0.5) ? p : 1.0-p;   // make between 0 and 0.5

   // if original value is invalid, return +infinity or -infinity
   // changed from original code to reflect the fact that the
   // the inverse of P(x) not Q(x) is desired.
   //
   // Original line: used for inverse of Q(x)
   // if ( dp <= 0.0 ){ dx = 13.0 ;  return (p <= 0.5) ? (dx) : (-dx); }

   // replaced with this if construct for the inverse of P(x)
   if (p <= 0.0)
     return -std::numeric_limits<double>::infinity();
   else if (p >= 1.0)
     return std::numeric_limits<double>::infinity();


   //  Step 1:  use 26.2.23 from Abramowitz and Stegun

   dt = std::sqrt( -2.0 * std::log(dp) ) ;
   dx = dt
     - ((.010328e+0*dt + .802853e+0)*dt + 2.515517e+0)
     /(((.001308e+0*dt + .189269e+0)*dt + 1.432788e+0)*dt + 1.e+0) ;

   //  Step 2:  do 3 Newton steps to improve this

   for ( newt=0 ; newt < 3 ; newt++ ){
     dq  = 0.5e+0 * ERFCC( dx / 1.414213562373095e+0 ) - dp ;
     ddq = std::exp( -0.5e+0 * dx * dx ) / 2.506628274631000e+0 ;
     dx  = dx + dq / ddq ;
   }

   // original line when computing the inverse of Q(x)
   // return (p <= 0.5) ? dx : (-dx)   // return with correct sign
   //
   // Note that P(-x) = Q(x), so whatever x was calculated for Q(x) = p,
   // we simply need to return the negative of x to get P(xp) = p.
   //
   return (p <= 0.5) ? (-dx) : dx; // return with correct sign
}

double
MuseTable::calculate_expected( unsigned int k, unsigned int n ) const
{
  return GaussianCDFInv(0.5*(1.0+((double)k / (double)(n+1))));
}


double
MuseTable::calculate_standard_dev( unsigned int k, unsigned int n,
                                         double expected_kth ) const
{
  double pk, qk, Qk, pQk, Qk_prime, Qk_dprime, Qk_tprime, vrk;

  pk = (double) k / (double) (n+1); // might want alpha beta correction
  qk = 1.0 - pk;

  // calculate_ the inverse cdf (might want to an alpha-beta correction)
  Qk = expected_kth;  // ak(k, N);   // inverse cdf of absolute residuals

  // density of absolute residual evaluated at Qk
  pQk = std::exp( -0.5 * Qk*Qk) * std::sqrt(2.0 / M_PI);

  // first derivative of Qk
  Qk_prime = 1.0/pQk;

  /*
  //  Low order approximation
  vrk = (pk*qk/(double)(n+2)) * Qk_prime*Qk_prime;
  */

  // second derivative of Qk
  Qk_dprime = Qk/(pQk*pQk);

  // third derivative of Qk
  Qk_tprime = ( 1.0 + 2.0 * Qk*Qk ) / (pQk*pQk*pQk);

  //  Higher order approximation
  vrk = (pk*qk/(double)(n+2)) * Qk_prime*Qk_prime
        + (pk*qk/((double)((n+2)*(n+2)))) * ( 2.0*(qk - pk)*Qk_prime*Qk_dprime
        + pk*qk*(Qk_prime*Qk_tprime + 0.5*Qk_dprime*Qk_dprime));

  return std::sqrt(vrk);
}


double
MuseTable::calculate_divisor( unsigned int k, unsigned int n,
                                    double expected_kth ) const
{
  return (n+1)*std::sqrt(2/M_PI)*(1.0-std::exp(-(expected_kth * expected_kth)/2.0));
}

double
MuseTable::calculate_sq_divisor( unsigned int k, unsigned int n,
                                       double expected_kth ) const
{
  return k - (n+1) * expected_kth * std::sqrt(2/M_PI)
    * std::exp(-(expected_kth*expected_kth)/2.0);
}
