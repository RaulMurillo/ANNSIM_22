#include <math.h>
// #include <stdio.h>      /* printf, scanf, puts, NULL */
# include <iostream>
#include <random>
#include "fft.hpp"
#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

using namespace sw::universal;

int main()
{
    using Real = posit<8,2>;
    /* initialize random seed: */
    std::random_device rd{};
    std::mt19937 gen{rd()}; gen.seed(1);
    std::normal_distribution<> d{0,1};

    // https://www.educba.com/matlab-fft/
    int Ls = 1024;  // Signal length
    // int Fs = 2000;  // Sampling frequency
    int Fs = 1000;  // Sampling frequency
    double Ts = 1.0/Fs;   // Sampling period
    double* tv = new double[Ls];   // Time vector
    double* f = new double[Ls];  // Input signal
    for (int i = 0; i< Ls; ++i)
    {
        tv[i] = i*Ts;
        // Add Random Noise to Signal { 0.6*sin(2*PI*50*tv[i]) + sin(2*PI*120*tv[i]) }
        // f[i] = 0.6*sin(2*PI*50*tv[i]) + sin(2*PI*120*tv[i]) + 3*d(gen);
        f[i] = 0.7*sin(2*PI*50*tv[i]) + sin(2*PI*120*tv[i]) + 2*d(gen);
    }
    std::cout << tv[0];
    for (int i = 1; i< Ls; ++i)
    {
        std::cout << "," << tv[i];
    }
    std::cout << std::endl;
    std::cout << f[0];
    for (int i = 1; i< Ls; ++i)
    {
        std::cout << "," << f[i];
    }
    std::cout << std::endl;

    complex<Real> *x = new complex<Real>[Ls];
    complex<Real> *DFT = new complex<Real>[Ls];
    for (int i = 0; i< Ls; ++i)
    {
        x[i].Re = f[i];
        x[i].Im = 0.0;
    }
    

    // Calling fft() function for signal ‘f’
    rad2FFT(Ls, x, DFT);

    std::cout << DFT[0].Re;
    for (int i = 1; i< Ls; ++i)
    {
        std::cout << "," << DFT[i].Re;
    }
    std::cout << std::endl;
    std::cout << DFT[0].Im;
    for (int i = 1; i< Ls; ++i)
    {
        std::cout << "," << DFT[i].Im;
    }
    std::cout << std::endl;



    return 0;
}