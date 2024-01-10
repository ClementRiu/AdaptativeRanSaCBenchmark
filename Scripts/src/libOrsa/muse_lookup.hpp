// Adapted from https://github.com/vxl/vxl/tree/c6c899aaf9cbf7ccabfca1ba35b2248fb611ffbc
// Algorithm from MUSE: Robust Surface Fitting using Unbiased Scale Estimates, James V. Miller and Charles V. Stewart
// By Clément Riu
// 2021

#ifndef MUSE_LOOKUP_H
#define MUSE_LOOKUP_H

#include <iostream>
#include <map>

//: Look-up table for the MUSET objective function.
//  Look-up table for the MUSET objective function, derived in James
//  V. Miller's 1997 PhD dissertation at Rensselaer.  An earlier
//  version of appeared in CVPR 1996.  The class computes and stores
//  statistics on the order statistics of Gaussian random variates.
//  Actually, these are for the order statistics of the absolute
//  values of Gaussian random variates.  See rrel_muset_obj for more
//  details.

class MuseKeyType
{
 public:
  MuseKeyType( unsigned int k, unsigned int n ) : k_(k), n_(n) {}
  unsigned int k_;
  unsigned int n_;
};

bool operator< ( MuseKeyType const& left_t, MuseKeyType const& right_t );


class MuseTableEntry
{
 public:
   MuseTableEntry() = default;
   bool initialized_{false};
   double expected_;
   double standard_dev_;
   double muse_t_divisor_;
   double muse_t_sq_divisor_;
};

typedef std::map< MuseKeyType, MuseTableEntry > MuseMapType;


class MuseTable
{
 public:
  //: Constructor.
  //  \a table_size is the size of table (= max number of residuals
  //  pre-computed).
  MuseTable( unsigned int /* max_n_stored */ ) {}

  MuseTable( ) = default;

  //: Destructor
  ~MuseTable() = default;

  //: Expected value of the kth ordered residual from n samples.
  //  The value is retrieved from the lookup table when possible.
  double expected_kth( unsigned int k, unsigned int n );

  //: Standard deviation of the kth ordered residual from n samples.
  //  The value is retrieved from the lookup table when possible.
  double standard_dev_kth( unsigned int k, unsigned int n );

  //: The divisor for trimmed statistics.
  //  The value is retrieved from the lookup table when possible.
  double muset_divisor( unsigned int k, unsigned int n );


  //: The divisor for trimmed square statistics.
  //  The value is retrieved from the lookup table when possible.
  double muset_sq_divisor( unsigned int k, unsigned int n );

 private:
  void calculate_all( unsigned int k, unsigned int n, MuseTableEntry & entry );

  //: Expected value of the kth ordered residual from n samples.
  //  The value is computed "from scratch".
  double calculate_expected( unsigned int k, unsigned int n ) const;

  //: Standard deviation of the kth ordered residual from n samples.
  //  The value is computed "from scratch".
  double calculate_standard_dev( unsigned int k, unsigned int n, double expected_kth ) const;

  //: The divisor for trimmed statistics.
  //  The value is computed "from scratch".
  double calculate_divisor( unsigned int k, unsigned int n, double expected_kth ) const;

  //: The divisor for trimmed squared statistics.
  //  The value is computed "from scratch".
  double calculate_sq_divisor( unsigned int k, unsigned int n, double expected_kth ) const;

 private:
  MuseMapType _table;
};

#endif // MuseTable_h_
