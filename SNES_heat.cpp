#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"

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

// % =========================================================================
// % This is the shape function routine for one-dimensional finite element code
// %
// % In this function, we give the Lagrange type basis function over the
// % computational domain [-1,1]. Interior nodes are uniformly distributed in
// % the reference domain.
// %
// % degree: the degree of our interpolation basis function space. In our code
// %         we give the basis function up to degree six.
// % i     : the number of the basis function. i takes value from 1 to degree+1.
// % der   : if der == 0, return the value of the basis function;
// %         if der == 1, return the 1st derivative of the basis function;
// %         if der == 2, return the 2nd derivative of the basis funciton.
// % x     : the point we which to perform the evaluation.
// %
// % Output: the value of the basis function or the 1st and 2nd derivatives of 
// %         the basis function.
// % -------------------------------------------------------------------------
// % By Ju Liu, 2009 Dec. 24th.
// % =========================================================================

void PolyBasis(int degree , int i , int der , double x, double &poly){
    // linear basis function
    if (degree == 1){
        if (i == 1){
            if (der == 0) poly = 0.5 * (1.0 - x);
            else if (der == 1) poly = -0.5;
            else if (der == 2) poly = 0.0;
        }
        else if (i == 2){
            if (der == 0) poly = 0.5 * (1.0 + x);
            else if (der == 1) poly = 0.5;
            else if (der == 2) poly = 0.0;
        }
    }
    // quadratic basis function
    else if (degree == 2){
        if (i == 1){
            if (der == 0) poly = 0.5 * x * (x - 1.0);
            else if (der == 1) poly = x - 0.5;
            else if (der == 2) poly = 1.0;
        }
        else if (i == 2){
            if (der == 0) poly = 1.0 - x*x;
            else if (der == 1) poly = -2.0 * x;
            else if (der == 2) poly = -2.0;
        }
        else if (i == 3){
            if (der == 0) poly = 0.5 * x * (x + 1.0);
            else if (der == 1) poly = x + 0.5;
            else if (der == 2) poly = 1.0;
        }
    }
    // cubic basis function
    else if (degree == 3){
        if (i == 1){
            if (der == 0) poly = -9.0 *( x - (1.0/3.0) ) * (x + (1.0/3.0) ) * (x - 1.0) / 16.0;
            else if (der == 1) poly = -9.0*(2.0*x*(x-1.0)+x*x-(1.0/9.0))/16.0;
            else if (der == 2) poly = -27.0/8.0*x+9.0/8.0;
        }
        else if (i == 2){
            if (der == 0) poly = 27.0 *(x*x-1.0)*(x-(1.0/3.0))/16.0;
            else if (der == 1) poly = 27.0 * (2.0*x*(x-(1.0/3.0))+x*x-1.0)/16.0;
            else if (der == 2) poly = 81.0/8.0*x-9.0/8.0;
        }
        else if (i == 3){
            if (der == 0) poly = -27.0 * (x*x-1.0)*(x+(1.0/3.0))/16.0;
            else if (der == 1) poly =  -27.0 * (2.0*x*(x+(1.0/3.0))+x*x-1.0)/16.0;
            else if (der == 2) poly = -81.0/8.0*x-9.0/8.0;
        }
        else if (i == 4){
            if (der == 0) poly = 9.0*(x+1.0)*(x*x-(1.0/9.0))/16.0;
            else if (der == 1) poly = 9.0*(x*x-(1.0/9.0)+2.0*x*(x+1.0))/16.0;
            else if (der == 2) poly = 27.0/8.0*x+9.0/8.0;
        }
    }
    // quartic basis function
    else if (degree == 4){
        if (i == 1){
            if (der == 0) poly = 2.0*x*(x*x-(1.0/4.0))*(x-1.0)/3.0;
            else if (der == 1) poly = 2.0*((x*x-(1.0/4.0))*(x-1.0)+2.0*x*x*(x-1.0)+x*(x*x-(1.0/4.0)))/3.0;
            else if (der == 2) poly = 4.0*x*(x-1.0)+4.0*x*x-1.0/3.0;
        }
        else if (i == 2){
            if (der == 0) poly = -8.0*x*(x*x-1.0)*(x-0.5)/3.0;
            else if (der == 1) poly = -8.0*((x*x-1.0)*(x-0.5)+x*x*(2.0*x-1.0)+x*(x*x-1.0))/3.0;
            else if (der == 2) poly = -16.0/3.0*x*(x-1.0/2.0)-16.0*x*x+16.0/3.0-16.0/3.0*x*(2.0*x-1.0);
        }
        else if (i == 3){
            if (der == 0) poly = 4.0*(x*x-1.0)*(x*x-0.25);
            else if (der == 1) poly = 4.0*(2.0*x*(x*x-0.25)+2.0*x*(x*x-1.0));
            else if (der == 2) poly = 48.0*x*x-10.0;
        }
        else if (i == 4){
            if (der == 0) poly = -8.0*x*(x*x-1.0)*(x+0.5)/3.0;
            else if (der == 1) poly = -8.0*((x*x-1.0)*(x+0.5)+x*x*(2.0*x+1.0)+x*(x*x-1.0))/3.0;
            else if (der == 2) poly = -16.0/3.0*x*(x+1.0/2.0)-16.0*x*x+16.0/3.0-16.0/3.0*x*(2.0*x+1.0);
        }
        else if (i == 5){
            if (der == 0) poly = 2.0*x*(x*x-0.25)*(x+1.0)/3.0;
            else if (der == 1) poly = 2.0*((x*x-0.25)*(x+1.0)+2.0*x*x*(x+1.0)+x*(x*x-0.25))/3.0;
            else if (der == 2) poly = 4.0*x*(x+1.0)+4.0*x*x-1.0/3.0;
        }
    }
    // quintic basis function
    else if (degree == 5){
        if (i == 1){
            if (der == 0) poly = -625.0*(x*x-(9.0/25.0))*(x*x-(1.0/25.0))*(x-1.0)/768.0;
            else if (der == 1) poly = -3125.0/768.0*x*x*x*x+625.0/192.0*x*x*x+125.0/128.0*x*x-125.0/192.0*x-3.0/256.0;
            else if (der == 2) poly = -3125.0/192.0*x*x*x+625.0/64.0*x*x+125.0/64.0*x-125.0/192.0;
        }
        else if (i == 2){
            if (der == 0) poly = 3125.0/768.0*(x+1.0)*(x+1.0/5.0)*(x-1.0/5.0)*(x-3.0/5.0)*(x-1.0);
            else if (der == 1) poly = 15625.0/768.0*x*x*x*x-625.0/64.0*x*x*x-1625.0/128.0*x*x+325.0/64.0*x+125.0/768.0;
            else if (der == 2) poly = 15625.0/192.0*x*x*x-1875.0/64.0*x*x-1625.0/64.0*x+325.0/64.0;
        }
        else if (i == 3){
            if (der == 0) poly = -3125.0/384.0*(x+1.0)*(x+3.0/5.0)*(x-1.0/5.0)*(x-3.0/5.0)*(x-1.0);
            else if (der == 1) poly = -15625.0/384.0*x*x*x*x+625.0/96.0*x*x*x+2125.0/64.0*x*x-425.0/96.0*x-375.0/128.0;
            else if (der == 2) poly = -15625.0/96.0*x*x*x+625.0/32.0*x*x+2125.0/32.0*x-425.0/96.0;
        }
        else if (i == 4){
            if (der == 0) poly = 3125.0/384.0*(x+1.0)*(x+3.0/5.0)*(x+1.0/5.0)*(x-3.0/5.0)*(x-1.0);
            else if (der == 1) poly = 15625.0/384.0*x*x*x*x+625.0/96.0*x*x*x-2125.0/64.0*x*x-425.0/96.0*x+375.0/128.0;
            else if (der == 2) poly = 15625.0/96.0*x*x*x+625.0/32.0*x*x-2125.0/32.0*x-425.0/96.0;
        }
        else if (i == 5){
            if (der == 0) poly = -3125.0/768.0*(x+1.0)*(x+3.0/5.0)*(x+1.0/5.0)*(x-1.0/5.0)*(x-1.0);
            else if (der == 1) poly = -15625.0/768.0*x*x*x*x-625.0/64.0*x*x*x+1625.0/128.0*x*x+325.0/64.0*x-125.0/768.0;
            else if (der == 2) poly = -15625.0/192.0*x*x*x-1875.0/64.0*x*x+1625.0/64.0*x+325.0/64.0;
        }
        else if (i == 6){
            if (der == 0) poly = 625.0/768.0*(x+1.0)*(x+3.0/5.0)*(x+1.0/5.0)*(x-1.0/5.0)*(x-3.0/5.0);
            else if (der == 1) poly = 3125.0/768.0*x*x*x*x-125.0/128.0*x*x+3.0/256.0+625.0/192.0*x*x*x-125.0/192.0*x;
            else if (der == 2) poly = 3125.0/192.0*x*x*x-125.0/64.0*x+625.0/64.0*x*x-125.0/192.0;
        }
    }
    // 
    else if (degree == 6){
        if (i == 1){
            if (der == 0) poly = 81.0/80.0*(x+2.0/3.0)*(x+1.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0) ;
            else if (der == 1) poly = 1.0/80.0*(6.0*x-1.0)*(81.0*x*x*x*x-54.0*x*x*x-39.0*x*x+16.0*x+4.0);
            else if (der == 2) poly = 243.0/40.0*x*x*x*x-81.0/20.0*x*x*x-117.0/40.0*x*x+6.0/5.0*x+3.0/10.0+(3.0/40.0*x-1.0/80.0)*(324.0*x*x*x-162.0*x*x-78.0*x+16.0);
        }
        else if (i == 2){
            if (der == 0) poly = -243.0/40.0*(x+1.0)*(x+1.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
            else if (der == 1) poly = -27.0/20.0*x-27.0/2.0*x*x+81.0/4.0*x*x*x*x-729.0/20.0*x*x*x*x*x+27.0*x*x*x+9.0/20.0;
            else if (der == 2) poly = -27.0/20.0-27.0*x+81.0*x*x*x-729.0/4.0*x*x*x*x+81.0*x*x;
        }
        else if (i == 3){
            if (der == 0) poly = 243.0/16.0*(x+1.0)*(x+2.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
            else if (der == 1) poly = 27.0/2.0*x+351.0/16.0*x*x-405.0/16.0*x*x*x*x-9.0/4.0+729.0/8.0*x*x*x*x*x-351.0/4.0*x*x*x;
            else if (der == 2) poly = 27.0/2.0+351.0/8.0*x-405.0/4.0*x*x*x+3645.0/8.0*x*x*x*x-1053.0/4.0*x*x;
        }
        else if (i == 4){
            if (der == 0) poly = -81.0/4.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
            else if (der == 1) poly = -49.0/2.0*x-243.0/2.0*x*x*x*x*x+126.0*x*x*x;
            else if (der == 2) poly = -49.0/2.0-1215.0/2.0*x*x*x*x+378.0*x*x;
        }
        else if (i == 5){
            if (der == 0) poly = 243.0/16.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*x*(x-2.0/3.0)*(x-1.0);
            else if (der == 1) poly = 27.0/2.0*x-351.0/16.0*x*x+405.0/16.0*x*x*x*x+729.0/8.0*x*x*x*x*x+9.0/4.0-351.0/4.0*x*x*x;
            else if (der == 2) poly = 27.0/2.0-351.0/8.0*x+405.0/4.0*x*x*x+3645.0/8.0*x*x*x*x-1053.0/4.0*x*x;
        }
        else if (i == 6){
            if (der == 0) poly = -243.0/40.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*x*(x-1.0/3.0)*(x-1.0);
            else if (der == 1) poly = -27.0/20.0*x+27.0/2.0*x*x-81.0/4.0*x*x*x*x-729.0/20.0*x*x*x*x*x+27.0*x*x*x-9.0/20.0;
            else if (der == 2) poly = -27.0/20.0+27.0*x-81.0*x*x*x-729.0/4.0*x*x*x*x+81.0*x*x;
        }
        else if (i == 7){
            if (der == 0) poly = 81.0/80.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0);
            else if (der == 1) poly = 1.0/80.0*(6.0*x+1.0)*(81.0*x*x*x*x+54.0*x*x*x-39.0*x*x-16.0*x+4.0);
            else if (der == 2) poly = 243.0/40.0*x*x*x*x+81.0/20.0*x*x*x-117.0/40.0*x*x-6.0/5.0*x+3.0/10.0+(3.0/40.0*x+1.0/80.0)*(324.0*x*x*x+162.0*x*x-78.0*x-16.0);
        }
    }
}

// exact solution 
double exact (double x){
    return sin(x);
}

double exact_x(double x){
    return cos(x);
}

// nonlinear kappa
double fun_kappa(double x){
    return 1.0  + x * x;
}

double fun_dkappa(double x){
    return 2*x;
}

double f_body(double x){
    return -2.0*cos(x)*cos(x)*sin(x) + sin(x)*(sin(x)*sin(x) + 1.0);
}

double h(double x){
    return -1.0 * fun_kappa( exact(0.0) );
}

double g(double x){
    return exact(1.0);
}

static char help[] = "Solves a nonlinear system with SNES.\n\n";

extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);

typedef struct {
    // domain
    double omega_l;
    double omega_r;

    std::vector<double> qp;
    std::vector<double> wq;
    std::vector<double> uh;
    std::vector<double> x_coor;

    std::vector<int> IEN;
    std::vector<int> ID; 

    int pp;                 
    int nElem;              
    int nqp;                 
    int n_np;   
    int n_en;         
    int n_eq; 

} AppCtx;
void AssemblyResidual(const PetscScalar * &xx, AppCtx &user);
void AssemblyJacobian(const PetscScalar * &xx, AppCtx &user);

int main(int argc, char **args) {

    int pp = 1;                 // interpolation degree
    int nElem = 10;              // number of elements
    int nqp = 2;                 // quadrature rule  
    int n_np = nElem * pp + 1;   // number of nodal points
    int n_en = pp + 1;           // number of element nodes

    // Setup the stiffness matrix and load vector 
    // number of equations equals the number of nodes minus the number of Dirichlet nodes
    int n_eq = n_np - 1; 

    AppCtx    user; /* user-defined work context */

    user.pp = pp;
    user.nElem = nElem;
    user.nqp = nqp;
    user.n_np = n_np;
    user.n_en = n_en;
    user.n_eq = n_eq;

    user.omega_l = 0.0;
    user.omega_r = 1.0;

    std::cout << "pp " << pp << std::endl;
    std::cout << "nElem " << nElem << std::endl;
    std::cout << "nqp " << nqp << std::endl;
    std::cout << "n_np " << n_np << std::endl;
    std::cout << "n_en " << n_en << std::endl;
    std::cout << "n_eq " << n_eq << std::endl;

    std::vector<double> qp(nqp, 0.0);
    std::vector<double> wq(nqp, 0.0);
    user.qp = qp;
    user.wq = wq;
    
    Gauss(nqp, -1.0, 1.0, user.qp, user.wq);

    std::cout << "qp: \n";
    for (int ii = 0; ii < nqp; ++ii){
        std::cout << user.qp[ii] << "\t";
    }
    std::cout << std::endl;

    std::cout << "wq: \n";
    for (int ii = 0; ii < nqp; ++ii){
        std::cout << user.wq[ii] << "\t";
    }
    std::cout << std::endl;

    std::vector<int> IEN(nElem * n_en, 0.0);
    std::vector<double> x_coor(n_np, 0.0);
    std::vector<int> ID(n_np, 0.0);
    user.IEN = IEN;
    user.x_coor = x_coor;
    user.ID = ID;

    for (int ee = 0; ee < nElem; ++ee){
        for (int aa = 0; aa < n_en; ++aa){
            user.IEN[ee + aa * nElem] = ee * pp + aa + 1;
        }
    }

    std::cout << "IEN: \n";
    for (int ii = 0; ii < nElem * n_en; ++ii){
        std::cout << user.IEN[ii] << "\t";
    }
    std::cout << std::endl;

    // mesh is assumed to have uniform size hh
    double hh = (user.omega_r - user.omega_l) / nElem;
    hh = hh / pp;

    for(int nx = 0; nx < n_np; ++nx){
        user.x_coor[nx] = nx * hh;
    }

    std::cout << "x_coor: \n";
    for (int ii = 0; ii < n_np; ++ii){
        std::cout << user.x_coor[ii] << "\t";
    }
    std::cout << std::endl;

    // setup ID array based on the boundary condition
    int counter = 1;
    for (int ii = 0; ii < n_np - 1; ++ii){
        user.ID[ii] = counter;
        counter = counter + 1;
    }

    std::cout << "ID: \n";
    for (int ii = 0; ii < n_np; ++ii){
        std::cout << user.ID[ii] << "\t";
    }
    std::cout << std::endl;

    // initial guess
    std::vector<double> uh(n_np, 0.0);
    uh[n_np-1] = g(user.omega_r);
    user.uh = uh;

    SNES           snes;         /* nonlinear solver context */
    KSP            ksp;          /* linear solver context */
    PC             pc;           /* preconditioner context */
    Vec            x,r;          /* solution, residual vectors */
    Mat            J;            /* Jacobian matrix */
    PetscErrorCode ierr;
    PetscMPIInt    size;
    PetscScalar    pfive = .5,*xx;
    PetscReal      rtol = 1.0e-9;  // 相对残差阈值
    PetscReal      atol = 1.0e-50;  // 绝对残差阈值
    PetscReal      stol = 1.0e-20;  // 范数残差
    PetscInt       maxit = 20;  // 最大迭代次数
    PetscInt       maxf = 10000; // 最大函数评估次数

    ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return ierr;
    ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Example is only for sequential runs");

     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

    ierr = SNESSetTolerances(snes, atol, rtol, PETSC_DECIDE, maxit, maxf); CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix and vector data structures; set corresponding routines
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
        Create vectors for solution and nonlinear function
    */
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,n_eq);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecDuplicate(x,&r);CHKERRQ(ierr);

    /*
        Create Jacobian matrix data structure
    */
    ierr = MatCreate(PETSC_COMM_WORLD,&J);CHKERRQ(ierr);
    ierr = MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n_eq,n_eq);CHKERRQ(ierr);
    ierr = MatSetFromOptions(J);CHKERRQ(ierr);
    ierr = MatSetUp(J);CHKERRQ(ierr);

    /*
    Set function evaluation routine and vector.
    */
    ierr = SNESSetFunction(snes,r,FormFunction,&user);CHKERRQ(ierr);

    // std::cout << "r: \n";
    // ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    /*
    Set Jacobian matrix data structure and Jacobian evaluation routine
    */
    // ierr = SNESSetJacobian(snes,J,J,FormJacobian,&user);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Customize nonlinear solver; set runtime options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    Set linear solver defaults for this problem. By extracting the
    KSP and PC contexts from the SNES context, we can then
    directly call any KSP and PC routines to set various options.
    */
    // ierr = SNESGetKSP(snes,&ksp);CHKERRQ(ierr);
    // ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
    // ierr = PCSetType(pc,PCNONE);CHKERRQ(ierr);
    // ierr = KSPSetTolerances(ksp,1.e-4,PETSC_DEFAULT,PETSC_DEFAULT,20);CHKERRQ(ierr);

    /*
    Set SNES/KSP/KSP/PC runtime options, e.g.,
        -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
    These options will override those specified above as long as
    SNESSetFromOptions() is called _after_ any other customization
    routines.
    */
    
    ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Evaluate initial guess; then solve nonlinear system
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    ierr = VecSet(x,0.0);CHKERRQ(ierr);
    
    /*
    Note: The user should initialize the vector, x, with the initial guess
    for the nonlinear solver prior to calling SNESSolve().  In particular,
    to employ an initial guess of zero, the user should explicitly set
    this vector to zero by calling VecSet().
    */
    // std::cout << "x: \n";
    // ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    ierr = SNESSolve(snes,NULL,x);CHKERRQ(ierr);

    Vec f;
    
    ierr = SNESGetFunction(snes,&f,0,0);CHKERRQ(ierr);
    std::cout << "x: \n";
    ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    std::cout << "f: \n";
    ierr = VecView(f,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // std::cout << "r: \n";
    // ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    // std::cout << "J: \n";
    // ierr = MatView(J,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = VecDestroy(&x);CHKERRQ(ierr); ierr = VecDestroy(&r);CHKERRQ(ierr);
    ierr = MatDestroy(&J);CHKERRQ(ierr); ierr = SNESDestroy(&snes);CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;

}

/* ------------------------------------------------------------------- */
/*
   FormFunction1 - Evaluates nonlinear function, F(x).

   Input Parameters:
.  snes - the SNES context
.  x    - input vector
.  ctx  - optional user-defined context

   Output Parameter:
.  f - function vector
 */
PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void *ctx)
{
    AppCtx  *user = (AppCtx *)ctx;
    PetscErrorCode    ierr;
    const PetscScalar *xx;
    PetscScalar       *F;

    /*
    Get pointers to vector data.
        - For default PETSc vectors, VecGetArray() returns a pointer to
            the data array.  Otherwise, the routine is implementation dependent.
        - You MUST call VecRestoreArray() when you no longer need access to
            the array.
    */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecGetArray(f,&F);CHKERRQ(ierr);

    VecZeroEntries(f);

    // AssemblyResidual(xx, *user);

    // Assembly stiffness matrix and load vector
    for(int ii = 0; ii < user->n_eq; ++ii){
        user->uh[ii] = xx[ii];
    }

    for (int ee = 0; ee < user->nElem; ++ee) {
        // Allocate zero element stiffness matrix 
        // std::vector<double> k_ele(user->n_en * user->n_en, 0.0);
        std::vector<double> f_ele(user->n_en, 0.0);

        std::vector<double> x_ele(user->n_en);
        std::vector<double> d_ele(user->n_en);
        std::vector<double> u_ele(user->n_en);
        for (int aa = 0; aa < user->n_en; ++aa) {
            x_ele[aa] = user->x_coor[user->IEN[ee + aa * user->nElem] - 1];
            // if (user->ID[user->IEN[ee + aa * user->nElem] - 1] == user->n_eq) u_ele[aa] = user->uh[user->n_np-1];
            u_ele[aa] = user->uh[user->IEN[ee + aa * user->nElem] - 1];
        }

        // loop over quadrature points
        for (int ll = 0; ll < user->nqp; ++ll) {
            // we need the geometric mapping at each quadrature points
            double x_qua = 0.0;
            double dx_dxi = 0.0;
            double u_qua = 0.0;
            double du_dxi = 0.0;
            for (int aa = 0; aa < user->n_en; ++aa) {
                double Na = 0.0;
                PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
                x_qua += x_ele[aa] * Na;
                u_qua += u_ele[aa] * Na;
                double Na_xi = 0.0;
                PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);
                dx_dxi += x_ele[aa] * Na_xi;
                du_dxi += u_ele[aa] * Na_xi;
            }
            double dxi_dx = 1.0 / dx_dxi;

            double kappa = fun_kappa( u_qua );
            double dkappa = fun_dkappa( u_qua );

            double detJ = dx_dxi;
            // we loop over a and b to assemble the element stiffness matrix and load vector
            for (int aa = 0; aa < user->n_en; ++aa){
                double Na = 0.0;
                PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
                double Na_xi = 0.0;
                PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);

                f_ele[aa] += user->wq[ll] * Na * f_body(x_qua) * detJ;
                f_ele[aa] -= user->wq[ll] * Na_xi * kappa * du_dxi * dxi_dx;

            }
        }
        // global assembly
        // distribute the entries to the global stiffness matrix and global load vector
        for (int aa = 0; aa < user->n_en; ++aa){
            int PP = user->ID[ user->IEN[aa * user->nElem + ee] - 1 ];
            if (PP > 0){
                F[PP - 1] += f_ele[aa];
            }
        }
        if (ee == 0){
            F[ user->ID[ user->IEN[0] - 1 ] - 1 ] += h( user->x_coor[ user->IEN[0] - 1 ]);
        }
        
    } //  end of element loop

    for (int ii = 0; ii < user->n_eq; ++ii){
        F[ii] = -F[ii];
    }
    
    // std::cout << "F: \n";
    // for (int ii = 0; ii < user->n_eq; ++ii){
    //     std::cout << F[ii] << "\n";
    // }
    // std::cout << std::endl;

    // std::cout << "xx: \n";
    // for (int ii = 0; ii < user->n_eq; ++ii){
    //     std::cout << xx[ii] << "\n";
    // }
    // std::cout << std::endl;


    /* Restore vectors */
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
    ierr = VecRestoreArray(f,&F);CHKERRQ(ierr);


    // ierr = VecView(f,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    return 0;
}


/* ------------------------------------------------------------------- */
/*
   FormJacobian1 - Evaluates Jacobian matrix.

   Input Parameters:
.  snes - the SNES context
.  x - input vector
.  dummy - optional user-defined context (not used here)

   Output Parameters:
.  jac - Jacobian matrix
.  B - optionally different preconditioning matrix
.  flag - flag indicating matrix structure
*/
PetscErrorCode FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void *ctx)
{
    AppCtx  *user = (AppCtx *)ctx;
    const PetscScalar *xx;
    PetscScalar       A[100];
    PetscErrorCode    ierr;
    //   PetscInt          idx[2] = {0,1};


    /*
        Get pointer to vector data
    */
    ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);

    MatZeroEntries(B);

    std::vector<double> K(user->n_eq * user->n_eq, 0.0);

    /*
        Compute Jacobian entries and insert into matrix.
        - Since this is such a small problem, we set all entries for
            the matrix at once.
    */
    // AssemblyJacobian(xx, *user);
    // Assembly stiffness matrix and load vector
    for(int ii = 0; ii < user->n_eq; ++ii){
        user->uh[ii] = xx[ii];
    }

    for (int ee = 0; ee < user->nElem; ++ee) {
        // Allocate zero element stiffness matrix 
        std::vector<double> k_ele(user->n_en * user->n_en, 0.0);
        // std::vector<double> f_ele(user->n_en, 0.0);

        std::vector<double> x_ele(user->n_en);
        std::vector<double> d_ele(user->n_en);
        std::vector<double> u_ele(user->n_en);
        for (int aa = 0; aa < user->n_en; ++aa) {
            x_ele[aa] = user->x_coor[user->IEN[ee + aa * user->nElem] - 1];
            // if (user->ID[user->IEN[ee + aa * user->nElem] - 1] == user->n_eq) u_ele[aa] = user->uh[user->n_np-1];
            u_ele[aa] = user->uh[user->IEN[ee + aa * user->nElem] - 1];
        }

        // loop over quadrature points
        for (int ll = 0; ll < user->nqp; ++ll) {
            // we need the geometric mapping at each quadrature points
            double x_qua = 0.0;
            double dx_dxi = 0.0;
            double u_qua = 0.0;
            double du_dxi = 0.0;
            for (int aa = 0; aa < user->n_en; ++aa) {
                double Na = 0.0;
                PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
                x_qua += x_ele[aa] * Na;
                u_qua += u_ele[aa] * Na;
                double Na_xi = 0.0;
                PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);
                dx_dxi += x_ele[aa] * Na_xi;
                du_dxi += u_ele[aa] * Na_xi;
            }
            double dxi_dx = 1.0 / dx_dxi;

            double kappa = fun_kappa( u_qua );
            double dkappa = fun_dkappa( u_qua );

            double detJ = dx_dxi;
            // we loop over a and b to assemble the element stiffness matrix and load vector
            for (int aa = 0; aa < user->n_en; ++aa){
                double Na = 0.0;
                PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
                double Na_xi = 0.0;
                PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);

                for (int bb = 0; bb < user->n_en; ++bb){
                    double Nb = 0.0;
                    PolyBasis(user->pp, bb + 1, 0, user->qp[ll], Nb);
                    double Nb_xi = 0.0;
                    PolyBasis(user->pp, bb + 1, 1, user->qp[ll], Nb_xi);
                    k_ele[aa * user->n_en + bb] += user->wq[ll] * Na_xi * kappa * Nb_xi * dxi_dx;
                    k_ele[aa * user->n_en + bb] += user->wq[ll] * Na_xi * dkappa * Nb * du_dxi * dxi_dx;
                }
            }
        }
        // global assembly
        // distribute the entries to the global stiffness matrix and global load vector
        for (int aa = 0; aa < user->n_en; ++aa){
            int PP = user->ID[ user->IEN[aa * user->nElem + ee] - 1 ];
            if (PP > 0){
                for (int bb = 0; bb < user->n_en; ++bb){
                    int QQ = user->ID[ user->IEN[bb * user->nElem + ee] - 1 ];
                    if (QQ > 0){
                        K[(PP - 1) * user->n_eq + QQ - 1] += k_ele[aa * user->n_en + bb];
                    }
                }
            }
        }
        
    } //  end of element loop

    // std::cout << "K: ";
    // for (int ii = 0; ii < user->n_eq * user->n_eq; ++ii){
    //     if (ii % user->n_eq == 0) std::cout << std::endl;
    //     std::cout << K[ii] << "\t";
    // }
    // std::cout << std::endl;

    for (PetscInt ii = 0; ii < user->n_eq; ++ii){
        for (PetscInt jj = 0; jj < user->n_eq; ++jj){
            ierr = MatSetValues(B, 1, &ii, 1, &jj, &K[ii * user->n_eq + jj], INSERT_VALUES); CHKERRQ(ierr);
        }
    }
        

    /*
        Restore vector
    */
    ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);

    /*
        Assemble matrix
    */
    ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    // ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    if (jac != B) {
        ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
        ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    }
    return 0;
}


void AssemblyResidual(const PetscScalar * &xx, AppCtx &user){

    std::vector<double> F(user.n_eq, 0.0);
    // Assembly stiffness matrix and load vector
    for(int ii = 0; ii < user.n_eq; ++ii){
        user.uh[ii] += xx[ii];
    }

    for (int ee = 0; ee < user.nElem; ++ee) {
        // Allocate zero element stiffness matrix 
        // std::vector<double> k_ele(user.n_en * user.n_en, 0.0);
        std::vector<double> f_ele(user.n_en, 0.0);

        std::vector<double> x_ele(user.n_en);
        std::vector<double> d_ele(user.n_en);
        std::vector<double> u_ele(user.n_en);
        for (int aa = 0; aa < user.n_en; ++aa) {
            x_ele[aa] = user.x_coor[user.IEN[ee + aa * user.nElem] - 1];
            u_ele[aa] = user.uh[user.IEN[ee + aa * user.nElem] - 1];
        }

        // loop over quadrature points
        for (int ll = 0; ll < user.nqp; ++ll) {
            // we need the geometric mapping at each quadrature points
            double x_qua = 0.0;
            double dx_dxi = 0.0;
            double u_qua = 0.0;
            double du_dxi = 0.0;
            for (int aa = 0; aa < user.n_en; ++aa) {
                double Na = 0.0;
                PolyBasis(user.pp, aa + 1, 0, user.qp[ll], Na);
                x_qua += x_ele[aa] * Na;
                u_qua += u_ele[aa] * Na;
                double Na_xi = 0.0;
                PolyBasis(user.pp, aa + 1, 1, user.qp[ll], Na_xi);
                dx_dxi += x_ele[aa] * Na_xi;
                du_dxi += u_ele[aa] * Na_xi;
            }
            double dxi_dx = 1.0 / dx_dxi;

            double kappa = fun_kappa( u_qua );
            double dkappa = fun_dkappa( u_qua );

            double detJ = dx_dxi;
            // we loop over a and b to assemble the element stiffness matrix and load vector
            for (int aa = 0; aa < user.n_en; ++aa){
                double Na = 0.0;
                PolyBasis(user.pp, aa + 1, 0, user.qp[ll], Na);
                double Na_xi = 0.0;
                PolyBasis(user.pp, aa + 1, 1, user.qp[ll], Na_xi);

                f_ele[aa] += user.wq[ll] * Na * f_body(x_qua) * detJ;
                f_ele[aa] -= user.wq[ll] * Na_xi * kappa * du_dxi * dxi_dx;
            }
        }
        // global assembly
        // distribute the entries to the global stiffness matrix and global load vector
        for (int aa = 0; aa < user.n_en; ++aa){
            int PP = user.ID[ user.IEN[aa * user.nElem + ee] - 1 ];
            if (PP > 0){
                F[PP - 1] += f_ele[aa];
            }
        }
        // Modify the load vector by the Natural BC
        // Note: for multi-dimensional cases, one needs to perform line or
        // surface integration for the natural BC.
        if (ee == 0){
            F[ user.ID[ user.IEN[0] - 1 ] - 1 ] += h( user.x_coor[ user.IEN[0] - 1 ]);
        }
        
    } //  end of element loop

    return ;

}


void AssemblyJacobian(const PetscScalar * &xx, AppCtx &user){

    std::vector<double> K(user.n_eq * user.n_eq, 0.0);
    // Assembly stiffness matrix and load vector
    for(int ii = 0; ii < user.n_eq; ++ii){
        user.uh[ii] += xx[ii];
    }

    for (int ee = 0; ee < user.nElem; ++ee) {
        // Allocate zero element stiffness matrix 
        std::vector<double> k_ele(user.n_en * user.n_en, 0.0);

        std::vector<double> x_ele(user.n_en);
        std::vector<double> d_ele(user.n_en);
        std::vector<double> u_ele(user.n_en);
        for (int aa = 0; aa < user.n_en; ++aa) {
            x_ele[aa] = user.x_coor[user.IEN[ee + aa * user.nElem] - 1];
            u_ele[aa] = user.uh[user.IEN[ee + aa * user.nElem] - 1];
        }

        // loop over quadrature points
        for (int ll = 0; ll < user.nqp; ++ll) {
            // we need the geometric mapping at each quadrature points
            double x_qua = 0.0;
            double dx_dxi = 0.0;
            double u_qua = 0.0;
            double du_dxi = 0.0;
            for (int aa = 0; aa < user.n_en; ++aa) {
                double Na = 0.0;
                PolyBasis(user.pp, aa + 1, 0, user.qp[ll], Na);
                x_qua += x_ele[aa] * Na;
                u_qua += u_ele[aa] * Na;
                double Na_xi = 0.0;
                PolyBasis(user.pp, aa + 1, 1, user.qp[ll], Na_xi);
                dx_dxi += x_ele[aa] * Na_xi;
                du_dxi += u_ele[aa] * Na_xi;
            }
            double dxi_dx = 1.0 / dx_dxi;

            double kappa = fun_kappa( u_qua );
            double dkappa = fun_dkappa( u_qua );

            double detJ = dx_dxi;
            // we loop over a and b to assemble the element stiffness matrix and load vector
            for (int aa = 0; aa < user.n_en; ++aa){
                double Na = 0.0;
                PolyBasis(user.pp, aa + 1, 0, user.qp[ll], Na);
                double Na_xi = 0.0;
                PolyBasis(user.pp, aa + 1, 1, user.qp[ll], Na_xi);

                for (int bb = 0; bb < user.n_en; ++bb){
                    double Nb = 0.0;
                    PolyBasis(user.pp, bb + 1, 0, user.qp[ll], Nb);
                    double Nb_xi = 0.0;
                    PolyBasis(user.pp, bb + 1, 1, user.qp[ll], Nb_xi);
                    k_ele[aa * user.n_en + bb] += user.wq[ll] * Na_xi * kappa * Nb_xi * dxi_dx;
                    k_ele[aa * user.n_en + bb] += user.wq[ll] * Na_xi * dkappa * Nb * du_dxi * dxi_dx;
                }
            }
        }
        // global assembly
        // distribute the entries to the global stiffness matrix and global load vector
        for (int aa = 0; aa < user.n_en; ++aa){
            int PP = user.ID[ user.IEN[aa * user.nElem + ee] - 1 ];
            if (PP > 0){
                for (int bb = 0; bb < user.n_en; ++bb){
                    int QQ = user.ID[ user.IEN[bb * user.nElem + ee] - 1 ];
                    if (QQ > 0){
                        K[(PP - 1) * user.n_eq + QQ - 1] += k_ele[aa * user.n_en + bb];
                    }
                }
            }
        }
        
    } //  end of element loop


    return ;

}
