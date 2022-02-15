# include <cmath>
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <omp.h>

#include <universal/number/posit/posit.hpp>

using namespace std;
using namespace sw::universal;

int main ( int argc, char *argv[] );
template<typename Real>
int main_md ( int argc, char *argv[], string dtype );
template<typename Real>
void compute ( int np, int nd, Real pos[], Real vel[], 
  Real mass, Real f[], Real *pot, Real *kin );
template<typename Real>
Real dist ( int nd, Real r1[], Real r2[], Real dr[] );
template<typename Real>
void initialize ( int np, int nd, Real box[], int *seed, Real pos[], 
  Real vel[], Real acc[] );
template<typename Real>
Real r8_uniform_01 ( int *seed );
template<typename Real>
void update ( int np, int nd, Real pos[], Real vel[], Real f[], 
  Real acc[], Real mass, Real dt );
void timestamp ( );

//****************************************************************************80

int main ( int argc, char *argv[] )
{
  // Double - FP64
  if(main_md< double >(argc, argv, "double") != 0)
    return 1;
  // Float - FP32
  if(main_md< float >(argc, argv, "float") != 0)
    return 1;
  // Posit<32,2>
  if(main_md< posit<32,2> >(argc, argv, "posit<32,2>") != 0)
    return 1;
  // Posit<16,2>
  if(main_md< posit<16,2> >(argc, argv, "posit<16,2>") != 0)
    return 1;

  return 0;
}

template<typename Real>
int main_md ( int argc, char *argv[], string dtype )

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for MD_OPENMP.
//
//  Discussion:
//
//    MD implements a simple molecular dynamics simulation.
//
//    The program uses Open MP directives to allow parallel computation.
//
//    The velocity Verlet time integration scheme is used. 
//
//    The particles interact with a central pair potential.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    Original FORTRAN90 version by Bill Magro.
//    C++ version by John Burkardt.
//
{
  using prec = long double;

  Real *acc;
  Real *box;
  Real dt = 0.0001;
  prec e0;
  Real *force;
  int i;
  Real kinetic;
  Real mass = 1.0;
  int nd = 3;     // default
  int np = 1000;  // default
  Real *pos;
  Real potential;
  int seed = 123456789;
  int step;
  int step_num = 400; // default
  int step_print;
  int step_print_index;
  int step_print_num;
  Real *vel;
  double wtime;

  timestamp ( );

//
//  Get the spatial dimension.
//
  if ( 1 < argc )
  {
    nd = atoi ( argv[1] );
  }
//
//  Get the number of particles.
//
  if ( 2 < argc )
  {
    np = atoi ( argv[2] );
  }
//
//  Get the number of time steps.
//
  if ( 3 < argc )
  {
    step_num = atoi ( argv[3] );
  }
//
//  Get the time step.
//
  if ( 4 < argc )
  {
    dt = atof ( argv[4] );
  }
//
//  Report.
//
  cout << "\n";
  cout << "MD_OPENMP\n";
  cout << "  C++/OpenMP version\n";
  cout << "\n";
  cout << "  A molecular dynamics program.\n";
  cout << "\n";
  cout << "  Data type: " << dtype << "\n";
  cout << "  ND, the spatial dimension, is " << nd << "\n";
  cout << "  NP, the number of particles in the simulation is " << np << "\n";
  cout << "  STEP_NUM, the number of time steps, is " << step_num << "\n";
  cout << "  DT, the size of each time step, is " << dt << "\n";;

  cout << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs ( ) << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads ( ) << "\n";
//
//  Allocate memory.
//
  acc = new Real[nd*np];
  box = new Real[nd];
  force = new Real[nd*np];
  pos = new Real[nd*np];
  vel = new Real[nd*np];
//
//  Set the dimensions of the box.
//
  for ( i = 0; i < nd; i++ )
  {
    box[i] = 10.0;
  }

  cout << "\n";
  cout << "  Initializing positions, velocities, and accelerations.\n" << flush;
//
//  Set initial positions, velocities, and accelerations.
//
  initialize ( np, nd, box, &seed, pos, vel, acc );
//
//  Compute the forces and energies.
//
  cout << "\n";
  cout << "  Computing initial forces and energies.\n";

  compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );

  e0 = (prec)(potential) + (prec)(kinetic);
//
//  This is the main time stepping loop:
//    Compute forces and energies,
//    Update positions, velocities, accelerations.
//
  cout << "\n";
  cout << "  At each step, we report the potential and kinetic energies.\n";
  cout << "  The sum of these energies should be a constant.\n";
  cout << "  As an accuracy check, we also print the relative error\n";
  cout << "  in the total energy.\n";
  cout << "\n";
  cout << "      Step      Potential       Kinetic        (P+K-E0)/E0\n";
  cout << "                Energy P        Energy K       Relative Energy Error\n";
  cout << "\n";

  step_print = 0;
  step_print_index = 0;
  step_print_num = 10;
  
  step = 0;
  cout << "  " << setw(8) << step
       << "  " << setw(14) << potential
       << "  " << setw(14) << kinetic
       << "  " << setw(14) << ( (prec)(potential) + (prec)(kinetic) - e0 ) / e0 << "\n";
  step_print_index = step_print_index + 1;
  step_print = ( step_print_index * step_num ) / step_print_num;

  wtime = omp_get_wtime ( );

  for ( step = 1; step <= step_num; step++ )
  {
    compute ( np, nd, pos, vel, mass, force, &potential, &kinetic );

    if ( step == step_print )
    {
      cout << "  " << setw(8) << step
           << "  " << setw(14) << potential
           << "  " << setw(14) << kinetic
           << "  " << setw(14) << ( (prec)(potential) + (prec)(kinetic) - e0 ) / e0 << "\n";
      step_print_index = step_print_index + 1;
      step_print = ( step_print_index * step_num ) / step_print_num;
    }
    update ( np, nd, pos, vel, force, acc, mass, dt );
  }

  wtime = omp_get_wtime ( ) - wtime;
  cout << "\n";
  cout << "  Elapsed cpu time for main computation:\n";
  cout << "  " << wtime << " seconds.\n";

  delete [] acc;
  delete [] box;
  delete [] force;
  delete [] pos;
  delete [] vel;
//
//  Terminate.
//
  cout << "\n";
  cout << "MD_OPENMP\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

template<typename Real>
void compute ( int np, int nd, Real pos[], Real vel[], 
  Real mass, Real f[], Real *pot, Real *kin )

//****************************************************************************80
//
//  Purpose:
//
//    COMPUTE computes the forces and energies.
//
//  Discussion:
//
//    The computation of forces and energies is fully parallel.
//
//    The potential function V(X) is a harmonic well which smoothly
//    saturates to a maximum value at PI/2:
//
//      v(x) = ( sin ( min ( x, PI2 ) ) )**2
//
//    The derivative of the potential is:
//
//      dv(x) = 2.0 * sin ( min ( x, PI2 ) ) * cos ( min ( x, PI2 ) )
//            = sin ( 2.0 * min ( x, PI2 ) )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    Original FORTRAN90 version by Bill Magro.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, Real POS[ND*NP], the position of each particle.
//
//    Input, Real VEL[ND*NP], the velocity of each particle.
//
//    Input, Real MASS, the mass of each particle.
//
//    Output, Real F[ND*NP], the forces.
//
//    Output, Real *POT, the total potential energy.
//
//    Output, Real *KIN, the total kinetic energy.
//
{
  Real d;
  Real d2;
  int i;
  int j;
  int k;
  Real ke;
  Real pe;
  Real PI2 = 3.141592653589793 / 2.0;
  Real rij[3];

  pe = 0.0;
  ke = 0.0;
  Real pe_v[np*np];
  Real ke_v[nd*np];

# pragma omp parallel \
  shared ( f, nd, np, pos, vel, pe_v, ke_v ) \
  private ( i, j, k, rij, d, d2 )
  

// Reduction is not supported with other types
// Store values in aux arrays and accumulate later
// # pragma omp for reduction ( + : pe, ke )
# pragma omp for
  for ( k = 0; k < np; k++ )
  {
//
//  Compute the potential energy and forces.
//
    for ( i = 0; i < nd; i++ )
    {
      f[i+k*nd] = 0.0;
    }

    for ( j = 0; j < np; j++ )
    {
      if ( k != j )
      {
        d = dist ( nd, pos+k*nd, pos+j*nd, rij );
//
//  Attribute half of the potential energy to particle J.
//
        if ( d < PI2 )
        {
          d2 = d;
        }
        else
        {
          d2 = PI2;
        }
// #pragma omp critical
//         {
//           pe = pe + 0.5 * pow ( sin ( d2 ), 2 );
//         }
        pe_v[j+k*np] = 0.5 * ( sin ( d2 )*sin ( d2 ) );

        for ( i = 0; i < nd; i++ )
        {
          f[i+k*nd] = f[i+k*nd] - rij[i] * sin ( 2.0 * d2 ) / d;
        }
      }
    }
//
//  Compute the kinetic energy.
//
// #pragma omp critical
    for ( i = 0; i < nd; i++ )
    {
      // ke = ke + vel[i+k*nd] * vel[i+k*nd];
      ke_v[i+k*nd] = vel[i+k*nd] * vel[i+k*nd];
    }
  }

  // Accumulate results
  for (k = 0; k < np; ++k)
  {
    for ( j = 0; j < np; ++j )
      pe = pe + pe_v[j+k*np];

    for ( i = 0; i < nd; ++i )
      ke = ke + ke_v[i+k*nd];
  }

  ke = ke * 0.5 * mass;
  
  *pot = pe;
  *kin = ke;

  return;
}
//****************************************************************************80

template<typename Real>
Real dist ( int nd, Real r1[], Real r2[], Real dr[] )

//****************************************************************************80
//
//  Purpose:
//
//    DIST computes the displacement (and its norm) between two particles.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    Original FORTRAN90 version by Bill Magro.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, Real R1[ND], R2[ND], the positions of the particles.
//
//    Output, Real DR[ND], the displacement vector.
//
//    Output, Real D, the Euclidean norm of the displacement.
//
{
  Real d;
  int i;

  d = 0.0;
  for ( i = 0; i < nd; i++ )
  {
    dr[i] = r1[i] - r2[i];
    d = d + dr[i] * dr[i];
  }
  d = sqrt ( d );

  return d;
}
//****************************************************************************80

template<typename Real>
void initialize ( int np, int nd, Real box[], int *seed, Real pos[], 
  Real vel[], Real acc[] )

//****************************************************************************80
//
//  Purpose:
//
//    INITIALIZE initializes the positions, velocities, and accelerations.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    Original FORTRAN90 version by Bill Magro.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input, Real BOX[ND], specifies the maximum position
//    of particles in each dimension.
//
//    Input, int *SEED, a seed for the random number generator.
//
//    Output, Real POS[ND*NP], the position of each particle.
//
//    Output, Real VEL[ND*NP], the velocity of each particle.
//
//    Output, Real ACC[ND*NP], the acceleration of each particle.
//
{
  int i;
  int j;
//
//  Give the particles random positions within the box.
//
  for ( i = 0; i < nd; i++ )
  {
    for ( j = 0; j < np; j++ )
    {
      pos[i+j*nd] = box[i] * r8_uniform_01<Real> ( seed );
    }
  }
//
//  Set the velocities.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      vel[i+j*nd] = 0.0;
    }
  }
//
//  Set the accelerations.
//
  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      acc[i+j*nd] = 0.0;
    }
  }

  return;
}
//****************************************************************************80

template<typename Real>
Real r8_uniform_01 ( int *seed )

//****************************************************************************80
//
//  Purpose:
//
//    R8_UNIFORM_01 is a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      unif = seed / ( 2**31 - 1 )
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
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//  Parameters:
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, Real R8_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  Real r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }

  r = ( Real ) (( *seed ) * 4.656612875E-10);

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
//    24 September 2003
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
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************80

template<typename Real>
void update ( int np, int nd, Real pos[], Real vel[], Real f[], 
  Real acc[], Real mass, Real dt )

//****************************************************************************80
//
//  Purpose:
//
//    UPDATE updates positions, velocities and accelerations.
//
//  Discussion:
//
//    The time integration is fully parallel.
//
//    A velocity Verlet algorithm is used for the updating.
//
//    x(t+dt) = x(t) + v(t) * dt + 0.5 * a(t) * dt**2
//    v(t+dt) = v(t) + 0.5 * ( a(t) + a(t+dt) ) * dt
//    a(t+dt) = f(t) / m
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 November 2007
//
//  Author:
//
//    Original FORTRAN90 version by Bill Magro.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int NP, the number of particles.
//
//    Input, int ND, the number of spatial dimensions.
//
//    Input/output, Real POS[ND*NP], the position of each particle.
//
//    Input/output, Real VEL[ND*NP], the velocity of each particle.
//
//    Input, Real F[ND*NP], the force on each particle.
//
//    Input/output, Real ACC[ND*NP], the acceleration of each particle.
//
//    Input, Real MASS, the mass of each particle.
//
//    Input, Real DT, the time step.
//
{
  int i;
  int j;
  Real rmass;

  rmass = 1.0 / mass;

# pragma omp parallel \
  shared ( acc, dt, nd, np, pos, rmass, vel ) \
  private ( i, j )

# pragma omp for

  for ( j = 0; j < np; j++ )
  {
    for ( i = 0; i < nd; i++ )
    {
      pos[i+j*nd] = pos[i+j*nd] + vel[i+j*nd] * dt + 0.5 * acc[i+j*nd] * dt * dt;
      vel[i+j*nd] = vel[i+j*nd] + 0.5 * dt * ( f[i+j*nd] * rmass + acc[i+j*nd] );
      acc[i+j*nd] = f[i+j*nd] * rmass;
    }
  }

  return;
}