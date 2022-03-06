# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;
# include "brownian_motion.hpp"

# include "universal/number/posit/posit.hpp"
# include "universal/number/cfloat/cfloat.hpp"

using namespace sw::universal;



template<typename Real>
int main_brownian (string type);

int main()
{
  string type;

	constexpr bool hasSubnormals   = true;
	constexpr bool hasSupernormals = false;
	constexpr bool isSaturating    = false;

  type = "double";
  main_brownian<double> (type);

  type = "float";
  main_brownian<float> (type);

  type = "float16";
  main_brownian<cfloat<16, 5, uint8_t, hasSubnormals, hasSupernormals, isSaturating>> (type);

  type = "bfloat";
  main_brownian<cfloat<16, 7, uint8_t, hasSubnormals, hasSupernormals, isSaturating>> (type);

  type = "posit_32";
  main_brownian<posit<32,2>> (type);

  type = "posit_28";
  main_brownian<posit<28,2>> (type);

  type = "posit_24";
  main_brownian<posit<24,2>> (type);

  type = "posit_20";
  main_brownian<posit<20,2>> (type);

  type = "posit_16";
  main_brownian<posit<16,2>> (type);

  type = "posit_14";
  main_brownian<posit<14,2>> (type);

  type = "posit_12";
  main_brownian<posit<12,2>> (type);

  type = "posit_10";
  main_brownian<posit<10,2>> (type);

  type = "posit_8";
  main_brownian<posit<8,2>> (type);


  return 0;
}


//****************************************************************************80

template<typename Real>
int main_brownian (string type)

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for BROWNIAN_MOTION_SIMULATION_TEST.
//
//  Discussion:
//
//    BROWNIAN_MOTION_SIMULATION_TEST tests the BROWNIAN_MOTION_SIMULATION library.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
{
  double d;
  double *dsq;
  string header;
  int k;
  int m;
  int n;
  int seed;
  double t;
  Real *x;

  timestamp ( );
  cout << "\n";
  cout << "BROWNIAN_MOTION_SIMULATION_TEST\n";
  cout << "  C++ version\n";
  cout << "  Test the BROWNIAN_MOTION_SIMULATION library.\n";
//
//  Compute the path of a particle undergoing Brownian motion.
//
  for ( m = 1; m <= 2; m++ )
  {
    n = 1001;
    d = 10.0;
    t = 1.0;
    seed = 123456789;
    x = brownian_motion_simulation<Real> ( m, n, d, t, seed );
    if ( m == 1 )
    {
      header = "motion_1d_" + type;
    }
    else if ( m == 2 )
    {
      header = "motion_2d_" + type;
    }
    brownian_motion_display ( m, n, x, header );
    delete [] x;
  }
// //
// //  Estimate the average displacement of the particle from the origin
// //  as a function of time.
// //
//   for ( m = 1; m <= 3; m++ )
//   {
//     k = 40;
//     n = 1001;
//     d = 10.0;
//     t = 1.0;
//     seed = 123456789;

//     dsq = brownian_displacement_simulation ( k, n, m, d, t, seed );
//     if ( m == 1 )
//     {
//       header = "displacement_1d";
//     }
//     else if ( m == 2 )
//     {
//       header = "displacement_2d";
//     }
//     else if ( m == 3 )
//     {
//       header = "displacement_3d";
//     }
//     brownian_displacement_display ( k, n, m, d, t, dsq, header );
//     delete [] dsq;
//   }
/*
  Terminate.
*/
  cout << "\n";
  cout << "BROWNIAN_MOTION_SIMULATION_TEST\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}