#ifndef GAUSS_HPP
#define GAUSS_HPP


#include <iostream>
#include <cmath>
#include <vector>

void Gauss(int N, double a, double b, std::vector<double> &x, std::vector<double> &w, const double EPS){
    N = N-1;
    int N1 = N+1;
    int N2 = N+2;

    std::vector<double> xu(N1, 0.0);
    for (int ii = 0; ii < N1; ++ii){
        xu[ii] = -1.0 + 2.0 * ii / N;
    }

    // Initial guess
    std::vector<double> y(N1, 0.0);
    for (int ii = 0; ii < N1; ++ii){
        y[ii] = cos( (2*ii + 1) * M_PI / (2*N+2) )+(0.27 / N1) * sin( M_PI * xu[ii] * N / N2 );
    }
    

    // Legendre-Gauss Vandermonde Matrix
    std::vector<double> L(N1 * N2, 0.0);

    // Derivative of LGM
    std::vector<double> Lp(N1 * N2, 0.0);
    std::vector<double> Lpp(N1, 0.0);

    // Compute the zeros of the N+1 Legendre Polynomial
    // using the recursion relation and the Newton-Raphson method

    std::vector<double> y0(N1, 0.0);

    // iterate until new points are uniformly within epsilon of old points
    while (true){
        for (int ii = 0; ii < N1; ++ii){
            L[N2 * ii] = 1.0;
            L[N2 * ii + 1] = y[ii];
            Lp[N2 * ii] = 0.0;
            Lp[N2 * ii + 1] = 1.0;
        }

        for (int kk = 1; kk < N1; ++kk){
            for( int ii = 0; ii < N1; ++ii){
                L[N2*ii + kk + 1] = ( (2 * (kk + 1) - 1) * y[ii] * L[N2 * ii + kk] - kk * L[N2*ii + kk - 1] )/ (kk + 1);
            }
        }

        for (int ii = 0; ii < N1; ++ii){
            Lpp[ii] = N2 * ( L[N2 * ii + N1 - 1] - y[ii] * L[N2 * ii + N2 - 1] ) / (1 - y[ii] * y[ii]);
        }


        for (int ii = 0; ii < N1; ++ii){
            y0[ii] = y[ii];
        }

        for (int ii = 0; ii < N1; ++ii){
            y[ii] = y0[ii] - L[N2 * ii + N2 - 1] / Lpp[ii];
        }

        double diff = 0.0;
        for (int ii = 0; ii < N1; ++ii){
            diff = std::max(diff, std::abs(y[ii] - y0[ii]));
        }

        if (diff < EPS) break;

    }

    // linear map from [-1, 1] to [a,b]
    for(int ii = 0; ii < N1; ++ii){
        x[ii] = (a * (1 - y[ii]) + b * (1 + y[ii])) / 2;
        w[ii] = ( (b - a) / ((1.0 - y[ii] * y[ii]) * Lpp[ii] * Lpp[ii]) ) * ( double(N2) / double(N1) ) * ( double(N2) / double(N1) ) ;
    }

}

#endif