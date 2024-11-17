#include <iostream>
#include <cmath>
#include <vector>


constexpr double EPS = 1e-14;

void Gauss(int N, double a, double b, std::vector<double> &x, std::vector<double> &w){
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
    std::cout << "y: \n";
    for (int ii = 0; ii < N1; ++ii){
        std::cout << y[ii] << "\t";
    }
    std::cout << std::endl;
    

    // Legendre-Gauss Vandermonde Matrix
    std::vector<double> L(N1 * N2, 0.0);

    // Derivative of LGM
    std::vector<double> Lp(N1 * N2, 0.0);
    std::vector<double> Lpp(N1, 0.0);

    // Compute the zeros of the N+1 Legendre Polynomial
    // using the recursion relation and the Newton-Raphson method

    std::vector<double> y0(N1, 0.0);

    int tag = 0;

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
        std::cout << "L: \n";
        for (int ii = 0; ii < N1; ++ii){
            for (int jj = 0; jj < N2; ++jj) std::cout << L[N2 * ii + jj] << "\t";
        }
        std::cout << std::endl;

        for (int ii = 0; ii < N1; ++ii){
            Lpp[ii] = N2 * ( L[N2 * ii + N1 - 1] - y[ii] * L[N2 * ii + N2 - 1] ) / (1 - y[ii] * y[ii]);
        }

        std::cout << "Lpp: \n";
        for (int ii = 0; ii < N1; ++ii){
            std::cout << Lpp[ii] << "\t";
        }
        std::cout << std::endl;

        for (int ii = 0; ii < N1; ++ii){
            y0[ii] = y[ii];
        }

        std::cout << "y0: \n";
        for (int ii = 0; ii < N1; ++ii){
            std::cout << y0[ii] << "\t";
        }
        std::cout << std::endl;

        for (int ii = 0; ii < N1; ++ii){
            y[ii] = y0[ii] - L[N2 * ii + N2 - 1] / Lpp[ii];
        }

        std::cout << "y: \n";
        for (int ii = 0; ii < N1; ++ii){
            std::cout << y[ii] << "\t";
        }
        std::cout << std::endl;

        double diff = 0.0;
        for (int ii = 0; ii < N1; ++ii){
            diff = std::max(diff, std::abs(y[ii] - y0[ii]));
        }

        std::cout << "diff: \n";
        std::cout << diff << std::endl;

        if (diff < EPS) break;
        // tag += 1;
        // if (tag > 50) break;

    }

    std::cout << "w: \n";
    for (int ii = 0; ii < N1; ++ii){
        std::cout << w[ii] << "\t";
    }
    std::cout << "tag: " << tag << std::endl;
    std::cout << std::endl;
    std::cout << "w[0] : " << ( (b - a) / ((1.0 - y[0] * y[0]) * Lpp[0] * Lpp[0]) ) * (N2 / N1) * (N2 / N1) << std::endl;
    // linear map from [-1, 1] to [a,b]
    for(int ii = 0; ii < N1; ++ii){
        x[ii] = (a * (1 - y[ii]) + b * (1 + y[ii])) / 2;
        w[ii] = ( (b - a) / ((1.0 - y[ii] * y[ii]) * Lpp[ii] * Lpp[ii]) ) * ( double(N2) / double(N1) ) * ( double(N2) / double(N1) ) ;
        std::cout << "w: " << w[ii] << std::endl;
    }

}

void Gauss2D(int N1, int N2, std::vector<double> &xi, std::vector<double> &eta, std::vector<double> &w) {
    // preallocation
    int N = N1 * N2;
    xi.resize(N, 0.0);
    eta.resize(N, 0.0);
    w.resize(N, 0);

    // generate 1D rule
    std::vector<double> x1(N1, 0.0);
    std::vector<double> w1(N1, 0.0);
    std::vector<double> x2(N2, 0.0);
    std::vector<double> w2(N2, 0.0);
    Gauss(N1, -1.0, 1.0, x1, w1);
    Gauss(N2, -1.0, 1.0, x2, w2);

    for (int ii = 0; ii < N1; ++ii){
        for (int jj = 0; jj < N2; ++jj){
            xi[jj * N1 + ii] = x1[ii];
            eta[jj * N1 + ii] = x2[jj];
            w[jj * N1 + ii] = w1[ii] * w2[jj];
        }
    }

}

int main(){
    int N = 3;
    std::vector<double> x(N, 0.0);
    std::vector<double> w(N, 0.0);
    double a = -1.0, b = 1.0;
    Gauss(N, a, b, x, w);
    std::cout << "x: \n";
    for (int ii = 0; ii < N; ++ii) std::cout << x[ii] << "\n";
    std::cout << "w: \n";
    for (int ii = 0; ii < N; ++ii) std::cout << w[ii] << "\n";

    int N1 = 3, N2 = 3;
    std::vector<double> xi(N1*N2, 0.0);
    std::vector<double> eta(N1*N2, 0.0);
    std::vector<double> w_2d(N1*N2, 0.0);
    Gauss2D(N1, N2, xi, eta, w_2d);
    std::cout << "xi: \n";
    for (int ii = 0; ii < N1; ++ii){
        for (int jj = 0; jj < N2; ++jj) std::cout << xi[ii * N2 + jj] << "\n";
    } 
    std::cout << "eta: \n";
    for (int ii = 0; ii < N1; ++ii){
        for (int jj = 0; jj < N2; ++jj) std::cout << eta[ii * N2 + jj] << "\n";
    } 
    std::cout << "w_2d: \n";
    for (int ii = 0; ii < N1; ++ii){
        for (int jj = 0; jj < N2; ++jj) std::cout << w_2d[ii * N2 + jj] << "\n";
    } 


}