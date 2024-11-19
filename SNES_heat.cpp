#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"


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

void Quad(int aa, double xi, double eta, double &val) {
    if (aa == 1) val = 0.25 * (1 - xi) * (1 - eta);
    else if (aa == 2) val =  0.25 * (1 + xi) * (1 - eta);
    else if (aa == 3) val =  0.25 * (1 + xi) * (1 + eta);
    else if (aa == 4) val =  0.25 * (1 - xi) * (1 + eta);
    else
    {
        std::cerr << "Error: value of a should be 1, 2, 3, or 4." << std::endl;
        std::exit(EXIT_FAILURE);  
    }
}

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
            if (der == 0) poly = 0.5 * x * (1.0 - x);
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
            if (der == 0) poly = ;
            else if (der == 1) poly = ;
            else if (der == 2) poly = ;
        }
        else if (i == 2){
            if (der == 0) poly = ;
            else if (der == 1) poly = ;
            else if (der == 2) poly = ;
        }
        else if (i == 3){
            if (der == 0) poly = ;
            else if (der == 1) poly = ;
            else if (der == 2) poly = ;
        }
        else if (i == 4){
            if (der == 0) poly = ;
            else if (der == 1) poly = ;
            else if (der == 2) poly = ;
        }
        else if (i == 5){
            if (der == 0) poly = ;
            else if (der == 1) poly = ;
            else if (der == 2) poly = ;
        }
    }
}




static char help[] = "Solves a nonlinear system with SNES.\n\n";

extern void FormFunction_SNES(void *ptr, Vec U, Vec F, void *ctx);
extern void FormJacobian_SNES(void *ptr, Vec U, Mat A, Mat B, void *ctx);

int main(int argc, char **args) {
    // domain
    double omega_l = 0.0;
    double omega_r = 1.0;

    // material properties and input data


    PetscErrorCode ierr;
    PetscFunctionBeginUser;

    ierr = PetscInitialize(&argc, &args, NULL, help); if (ierr) return ierr;

    // 创建 SNES 对象
    SNES snes;
    ierr = SNESCreate(PETSC_COMM_WORLD, &snes); CHKERRQ(ierr);

    // 设置函数评估和 Jacobian 计算的回调函数
    ierr = SNESSetFunction(snes, NULL, FormFunction_SNES, NULL); CHKERRQ(ierr);
    ierr = SNESSetJacobian(snes, NULL, NULL, FormJacobian_SNES, NULL); CHKERRQ(ierr);

    // 设置初始猜测解
    Vec U;
    ierr = VecCreate(PETSC_COMM_WORLD, &U); CHKERRQ(ierr);
    ierr = VecSetSizes(U, PETSC_DECIDE, n_eq); CHKERRQ(ierr);
    ierr = VecSetFromOptions(U); CHKERRQ(ierr);
    ierr = VecSet(U, 0.0); // 初始猜测为零

    // 配置 SNES
    ierr = SNESSetFromOptions(snes); CHKERRQ(ierr);
    ierr = SNESSetTolerances(snes, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT); CHKERRQ(ierr);

    // 求解非线性问题
    ierr = SNESSolve(snes, U); CHKERRQ(ierr);

    // 查看结果
    ierr = VecView(U, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);

    // 清理
    ierr = VecDestroy(&U); CHKERRQ(ierr);
    ierr = SNESDestroy(&snes); CHKERRQ(ierr);
    ierr = PetscFinalize();
    return ierr;
}

// 定义函数评估的回调函数
void FormFunction_SNES(void *ptr, Vec U, Vec F, void *ctx) {
    // 这里实现具体的函数评估逻辑
    // 根据 U 计算 F
}

// 定义 Jacobian 计算的回调函数
void FormJacobian_SNES(void *ptr, Vec U, Mat A, Mat B, void *ctx) {
    // 这里实现具体的 Jacobian 计算逻辑
    // 根据 U 计算 A 或 B
}