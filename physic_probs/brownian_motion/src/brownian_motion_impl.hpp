// https://people.math.sc.edu/Burkardt/cpp_src/brownian_motion_simulation/brownian_motion_simulation.html
# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>

using namespace std;

# include "brownian_motion.hpp"

//****************************************************************************80

// void brownian_displacement_display ( int k, int n, int m, double d, double t, 
//   double dsq[], string header )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    BROWNIAN_DISPLACEMENT_DISPLAY displays average Brownian motion displacement.
// //
// //  Discussion:
// //
// //    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
// //    displacement calculation.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    18 August 2018
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int K, the number of repetitions.
// //
// //    Input, int N, the number of time steps.  
// //
// //    Input, int M, the spatial dimension.
// //
// //    Input, double D, the diffusion coefficient.
// //
// //    Input, double T, the total time.
// //
// //    Input, double DSQ[K*N], the displacements over time for 
// //    each repetition.
// //
// //    Input, string HEADER, an identifier for the output files.
// //
// {
//   string command_filename;
//   ofstream command;
//   string data_filename;
//   ofstream data;
//   double dsq_ave;
//   double dsq_ideal;
//   int i;
//   int *ii;
//   int j;
//   int seed;
//   double ti;

//   seed = 123456789;
// //
// //  Choose 5 paths at random.
// //
//   ii = i4vec_uniform_new ( 5, 0, k - 1, seed );
// //
// //  Create the data file.
// //
//   data_filename = header + "_data.txt";

//   data.open ( data_filename.c_str ( ) );

//   for ( j = 0; j < n; j++ )
//   {
//     ti = ( double ) ( j ) * t / ( double ) ( n - 1 );
//     dsq_ave = 0.0;
//     for ( i = 0; i < k; i++ )
//     {
//       dsq_ave = dsq_ave + dsq[i+j*k];
//     }
//     dsq_ave = dsq_ave / ( double ) k;
//     dsq_ideal = 2.0 * m * d * ti;
//     data << "  " << ti
//          << "  " << dsq[ii[0]+j*k]
//          << "  " << dsq[ii[1]+j*k]
//          << "  " << dsq[ii[2]+j*k]
//          << "  " << dsq[ii[3]+j*k]
//          << "  " << dsq[ii[3]+j*k]
//          << "  " << dsq_ave
//          << "  " << dsq_ideal << "\n";
//   }

//   data.close ( );

//   cout << "\n";
//   cout << "  BROWNIAN_DISPLACEMENT data stored in \"" << data_filename << "\".\n";
// //
// //  Create the command file.
// //
//   command_filename = header + "_commands.txt";

//   command.open ( command_filename.c_str ( ) );

//   command << "# " << command_filename << "\n";
//   command << "#\n";
//   command << "# Usage:\n";
//   command << "#  gnuplot < " << command_filename << "\n";
//   command << "#\n";
//   command << "set term png\n";
//   command << "set output '" << header << ".png'\n";
//   command << "set xlabel 'T'\n";
//   command << "set ylabel 'D^2'\n";
//   command << "set title 'Squared displacement (Red), Predicted (Black), Samples (Blue)'\n";
//   command << "set grid\n";
//   command << "set style data lines\n";
//   command << "plot '" << data_filename << "' using 1:2 title 'sample 1' linecolor rgb 'blue', \\\n";
//   command << "     '" << data_filename << "' using 1:3 title 'sample 2' linecolor rgb 'blue', \\\n";
//   command << "     '" << data_filename << "' using 1:4 title 'sample 3' linecolor rgb 'blue', \\\n";
//   command << "     '" << data_filename << "' using 1:5 title 'sample 4' linecolor rgb 'blue', \\\n";
//   command << "     '" << data_filename << "' using 1:6 title 'sample 5' linecolor rgb 'blue', \\\n";
//   command << "     '" << data_filename << "' using 1:7 title 'Averaged' lw 3 linecolor rgb 'red', \\\n";
//   command << "     '" << data_filename << "' using 1:8 title 'Ideal' lw 3 linecolor rgb 'black'\n";
//   command << "quit\n";

//   command.close ( );

//   cout << "  BROWNIAN_DISPLACEMENT plot commands stored in \"" << command_filename << "\".\n";

//   delete [] ii;

//   return;
// }
// //****************************************************************************80

// double *brownian_displacement_simulation ( int k, int n, int m, double d, 
//   double t, int &seed )

// //****************************************************************************80
// //
// //  Purpose:
// //
// //    BROWNIAN_DISPLACEMENT_SIMULATION simulates Brownian displacement.
// //
// //  Discussion:
// //
// //    This function computes the square of the distance of the Brownian
// //    particle from the starting point, repeating this calculation 
// //    several times.
// //
// //  Licensing:
// //
// //    This code is distributed under the GNU LGPL license.
// //
// //  Modified:
// //
// //    02 October 2012
// //
// //  Author:
// //
// //    John Burkardt
// //
// //  Parameters:
// //
// //    Input, int K, the number of repetitions.
// //
// //    Input, int N, the number of time steps to take, plus 1.
// //
// //    Input, int M, the spatial dimension.
// //
// //    Input, double D, the diffusion coefficient.
// //    Computationally, this is simply a scale factor between time and space.
// //
// //    Input, double T, the total time.
// //
// //    Input/output, int &SEED, a seed for the random
// //    number generator.
// //
// //    Output, double BROWNIAN_DISPLACEMENT_SIMULATION[K*N], the displacements 
// //    over time for each repetition. 
// //
// {
//   double *dsq;
//   int i;
//   int i2;
//   int j;
//   double temp;
//   double *x;

//   dsq = new double[k*n];

//   for ( i = 0; i < k; i++ )
//   {
//     x = brownian_motion_simulation ( m, n, d, t, seed );

//     for ( j = 0; j < n; j++ )
//     {
//       temp = 0.0;
//       for ( i2 = 0; i2 < m; i2++ )
//       {
//         temp = temp + pow ( x[i2+j*m], 2 );
//       }
//       dsq[i+j*k] = temp;
//     }
//     delete [] x;
//   }

//   return dsq;
// }
// //****************************************************************************80

template<typename Real>
void brownian_motion_display ( int m, int n, Real x[], string header )

//****************************************************************************80
//
//  Purpose:
//
//    BROWNIAN_MOTION_DISPLAY displays successive Brownian motion positions.
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
//  Parameters:
//
//    Input, int M, the spatial dimension.
//    M should be 1 or 2.
//
//    Input, int N, the number of time steps. 
//
//    Input, double X[M*N], the particle positions.
//
//    Input, string HEADER, an identifier for the output files.
//
{
  string command_filename;
  ofstream command;
  string data_filename;
  ofstream data;
  int i;
  double t;

  if ( m != 1 && m != 2 )
  {
    cerr << "\n";
    cerr << "BROWNIAN_MOTION_DISPLAY - Fatal error!\n";
    cerr << "  This routine can only handle M = 1 or 2.\n";
    exit ( 1 );
  }
//
//  Create the data file.
//
  data_filename = header + "_data.txt";

  data.open ( data_filename.c_str ( ) );

  if ( m == 1 )
  {
    for ( i = 0; i < n; i++ )
    {
      t = ( double ) ( i - 1 ) / ( double ) ( n - 1 );
      data << "  " << t
           << "  " << x[0+i*m] << "\n";
    }
  }
  else if ( m == 2 )
  {
    for ( i = 0; i < n; i++ )
    {
      data << "  " << x[0+i*m]
           << "  " << x[1+i*m] << "\n";
    }
  }

  data.close ( );

  cout << "\n";
  cout << "  BROWNIAN_MOTION data stored in \"" << data_filename << "\n";
//
//  Create the command file.
//
  command_filename = header + "_commands.txt";

  command.open ( command_filename.c_str ( ) );

  command << "# " << command_filename << "\n";
  command << "#\n";
  command << "# Usage:\n";
  command << "#  gnuplot < " << command_filename << "\n";
  command << "#\n";
  command << "set term png\n";
  command << "set output '" << header << ".png'\n";
  command << "set xlabel 'X'\n";
  command << "set ylabel 'T'\n";
  command << "set title 'Brownian motion in 1D'\n";
  command << "set grid\n";
  command << "set style data lines\n";
  command << "plot '" << data_filename << "' using 1:2\n";
  command << "quit\n";

  command.close ( );

  cout << "  BROWNIAN_MOTION plot commands stored in \"" << command_filename << "\n";

  // return ;
}
//****************************************************************************80

template<typename Real>
Real *brownian_motion_simulation ( int m, int n, double d, double t, 
  int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    BROWNIAN_MOTION_SIMULATION simulates Brownian motion.
//
//  Discussion:
//
//    Thanks to Feifei Xu for pointing out a missing factor of 2 in the
//    stepsize calculation, 08 March 2016.
//
//    Thanks to Joerg Peter Pfannmoeller for pointing out a missing factor
//    of M in the stepsize calculation, 23 April 2018.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 April 2018
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int N, the number of time steps to take, plus 1. 
//
//    Input, double D, the diffusion coefficient.  
//
//    Input, double T, the total time.
//
//    Input/output, int &SEED, a seed for the random
//    number generator.
//
//    Output, double BROWNIAN_MOTION_SIMULATION[M*N], the initial position at 
//    time 0.0, and the N-1 successive locations of the particle.
//
{
  double dt;
  Real *dx;
  double *dxx;
  int i;
  int j;
  Real norm_dx;
  double s;
  Real *x;

  x = new Real[m*n];
  dx = new Real[m];
//
//  Set the time step.
//
  dt = t / ( double ) ( n - 1 );
//
//  Start at the origin.
//
  for ( i = 0; i < m; i++ )
  {
    x[i+0*m] = (Real)(0.0);
  }
//
//  Take N - 1 steps.
//
  for ( j = 1; j < n; j++ )
  {
//
//  S is the stepsize.
//
    s = sqrt ( 2 * m * d * dt ) * r8_normal_01 ( seed );
//
//  Direction DX is random, unit norm.
//
    dxx = r8vec_normal_01_new ( m, seed );

    // Cast to Real dtype
    for (i = 0; i < m; i++){
      dx[i] = (Real)(dxx[i]);
    }
    norm_dx = (Real)(0.0);
    for ( i = 0; i < m; i++ )
    {
      // norm_dx = norm_dx + pow ( dx[i], 2 );
      norm_dx = norm_dx + (dx[i] * dx[i]);
    }
    norm_dx = sqrt ( norm_dx );
    for ( i = 0; i < m; i++ )
    {
      dx[i] = (Real)(s) * dx[i] / norm_dx;
    }
//
//  Add the step to the current position.
//
    for ( i = 0; i < m; i++ )
    {
      x[i+j*m] = x[i+(j-1)*m] + dx[i];
    }
  }
  return x;
}
//****************************************************************************80

int *i4vec_uniform_new ( int n, int a, int b, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    I4VEC_UNIFORM_NEW returns a scaled pseudorandom I4VEC.
//
//  Discussion:
//
//    An I4VEC is a vector of I4's.
//
//    The pseudorandom numbers should be uniformly distributed
//    between A and B.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 May 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the dimension of the vector.
//
//    Input, int A, B, the limits of the interval.
//
//    Input/output, int &SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int IVEC_UNIFORM_NEW[N], a vector of random values between A and B.
//
{
  int c;
  int i;
  int i4_huge = 2147483647;
  int k;
  float r;
  int value;
  int *x;
  
  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "I4VEC_UNIFORM_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }
//
//  Guarantee A <= B.
//
  if ( b < a )
  {
    c = a;
    a = b;
    b = c;
  }

  x = new int[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r = ( float ) ( seed ) * 4.656612875E-10;
//
//  Scale R to lie between A-0.5 and B+0.5.
//
    r = ( 1.0 - r ) * ( ( float ) a - 0.5 ) 
      +         r   * ( ( float ) b + 0.5 );
//
//  Use rounding to convert R to an integer between A and B.
//
    value = round ( r );
//
//  Guarantee A <= VALUE <= B.
//
    if ( value < a )
    {
      value = a;
    }
    if ( b < value )
    {
      value = b;
    }

    x[i] = value;
  }

  return x;
}
//****************************************************************************80

double r8_normal_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_NORMAL_01 samples the standard normal probability distribution.
//
//  Discussion:
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    The Box-Muller method is used, which is efficient, but
//    generates two values at a time.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input/output, int SEED, a seed for the random number generator.
//
//    Output, double R8_NORMAL_01, a normally distributed random value.
//
{
  double pi = 3.141592653589793;
  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
//
//  If we've used an even number of values so far, generate two more, return one,
//  and save one.
//
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = r8_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }

    r2 = r8_uniform_01 ( seed );

    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {

    x = y;

  }

  used = used + 1;

  return x;
}
//****************************************************************************80

double r8_uniform_01 ( int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//
{
  int i4_huge = 2147483647;
  int k;
  double r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed / 127773;

  seed = 16807 * ( seed - k * 127773 ) - k * 2836;

  if ( seed < 0 )
  {
    seed = seed + i4_huge;
  }
  r = ( double ) ( seed ) * 4.656612875E-10;

  return r;
}
//****************************************************************************80

double *r8vec_normal_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_NORMAL_01_NEW returns a unit pseudonormal R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//    The standard normal probability distribution function (PDF) has
//    mean 0 and standard deviation 1.
//
//    This routine can generate a vector of values on one call.  It
//    has the feature that it should provide the same results
//    in the same order no matter how we break up the task.
//
//    Before calling this routine, the user may call RANDOM_SEED
//    in order to set the seed of the random number generator.
//
//    The Box-Muller method is used, which is efficient, but
//    generates an even number of values each time.  On any call
//    to this routine, an even number of new values are generated.
//    Depending on the situation, one value may be left over.
//    In that case, it is saved for the next call.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of values desired.  If N is negative,
//    then the code will flush its internal memory; in particular,
//    if there is a saved value to be used on the next call, it is
//    instead discarded.  This is useful if the user has reset the
//    random number seed, for instance.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_NORMAL_01_NEW[N], a sample of the standard normal PDF.
//
//  Local parameters:
//
//    Local, int MADE, records the number of values that have
//    been computed.  On input with negative N, this value overwrites
//    the return value of N, so the user can get an accounting of
//    how much work has been done.
//
//    Local, double R[N+1], is used to store some uniform random values.
//    Its dimension is N+1, but really it is only needed to be the
//    smallest even number greater than or equal to N.
//
//    Local, int SAVED, is 0 or 1 depending on whether there is a
//    single saved value left over from the previous call.
//
//    Local, int X_LO, X_HI, records the range of entries of
//    X that we need to compute.  This starts off as 1:N, but is adjusted
//    if we have a saved value that can be immediately stored in X(1),
//    and so on.
//
//    Local, double Y, the value saved from the previous call, if
//    SAVED is 1.
//
{
  int i;
  int m;
  static int made = 0;
  double pi = 3.141592653589793;
  double *r;
  static int saved = 0;
  double *x;
  int x_hi;
  int x_lo;
  static double y = 0.0;
//
//  I'd like to allow the user to reset the internal data.
//  But this won't work properly if we have a saved value Y.
//  I'm making a crock option that allows the user to signal
//  explicitly that any internal memory should be flushed,
//  by passing in a negative value for N.
//
  if ( n < 0 )
  {
    made = 0;
    saved = 0;
    y = 0.0;
    return NULL;
  }
  else if ( n == 0 )
  {
    return NULL;
  }

  x = new double[n];
//
//  Record the range of X we need to fill in.
//
  x_lo = 1;
  x_hi = n;
//
//  Use up the old value, if we have it.
//
  if ( saved == 1 )
  {
    x[0] = y;
    saved = 0;
    x_lo = 2;
  }
//
//  Maybe we don't need any more values.
//
  if ( x_hi - x_lo + 1 == 0 )
  {
  }
//
//  If we need just one new value, do that here to avoid null arrays.
//
  else if ( x_hi - x_lo + 1 == 1 )
  {
    r = r8vec_uniform_01_new ( 2, seed );

    x[x_hi-1] = sqrt ( -2.0 * log ( r[0] ) ) * cos ( 2.0 * pi * r[1] );
    y =         sqrt ( -2.0 * log ( r[0] ) ) * sin ( 2.0 * pi * r[1] );

    saved = 1;

    made = made + 2;

    delete [] r;
  }
//
//  If we require an even number of values, that's easy.
//
  else if ( ( x_hi - x_lo + 1 ) % 2 == 0 )
  {
    m = ( x_hi - x_lo + 1 ) / 2;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-2; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }
    made = made + x_hi - x_lo + 1;

    delete [] r;
  }
//
//  If we require an odd number of values, we generate an even number,
//  and handle the last pair specially, storing one in X(N), and
//  saving the other for later.
//
  else
  {
    x_hi = x_hi - 1;

    m = ( x_hi - x_lo + 1 ) / 2 + 1;

    r = r8vec_uniform_01_new ( 2*m, seed );

    for ( i = 0; i <= 2*m-4; i = i + 2 )
    {
      x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
      x[x_lo+i  ] = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );
    }

    i = 2*m - 2;

    x[x_lo+i-1] = sqrt ( -2.0 * log ( r[i] ) ) * cos ( 2.0 * pi * r[i+1] );
    y           = sqrt ( -2.0 * log ( r[i] ) ) * sin ( 2.0 * pi * r[i+1] );

    saved = 1;

    made = made + x_hi - x_lo + 2;

    delete [] r;
  }

  return x;
}
//****************************************************************************80

double *r8vec_uniform_01_new ( int n, int &seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01_NEW returns a new unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01_NEW[N], the vector of pseudorandom values.
//
{
  int i;
  int i4_huge = 2147483647;
  int k;
  double *r;

  if ( seed == 0 )
  {
    cerr << "\n";
    cerr << "R8VEC_UNIFORM_01_NEW - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  r = new double[n];

  for ( i = 0; i < n; i++ )
  {
    k = seed / 127773;

    seed = 16807 * ( seed - k * 127773 ) - k * 2836;

    if ( seed < 0 )
    {
      seed = seed + i4_huge;
    }

    r[i] = ( double ) ( seed ) * 4.656612875E-10;
  }

  return r;
}
//****************************************************************************80

void timestamp ( )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}