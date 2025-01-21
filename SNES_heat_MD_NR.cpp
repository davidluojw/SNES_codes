#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"
#include "UserPointers.hpp"
#include "Gauss.hpp"
#include "Res_Jac.hpp"

static char help[] = "Solves a nonlinear system with SNES.\n\n";

int main(int argc, char **args) {

    const double EPS = 1e-14; 

    double omega_l = 0.0;
    double omega_r = 10.0;

    int pp = 2;                 // interpolation degree
    int nElem = 10;              // number of elements
    int nqp = 2;                 // quadrature rule  
    int n_np = nElem * pp + 1;   // number of nodal points
    int n_en = pp + 1;           // number of element nodes

    // Setup the stiffness matrix and load vector 
    // number of equations equals the number of nodes minus the number of Dirichlet nodes
    int n_eq = n_np - 1; 

    UserPointers user; /* user-defined work context */

    user.pp = pp;
    user.nElem = nElem;
    user.nqp = nqp;
    user.n_np = n_np;
    user.n_en = n_en;
    user.n_eq = n_eq;

    user.omega_l = omega_l;
    user.omega_r = omega_r;

    user.EPS = EPS;

    std::cout << "omega_l " << omega_l << std::endl;
    std::cout << "omega_r " << omega_r << std::endl;
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
    
    Gauss(nqp, -1.0, 1.0, user.qp, user.wq, EPS);

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
    double hh = (omega_r - omega_l) / nElem;
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
    // PC             pc;           /* preconditioner context */
    SNESLineSearch linesearch;   /* line search context */
    Vec            x,r;          /* solution, residual vectors */
    Mat            J;            /* Jacobian matrix */
    PetscErrorCode ierr;
    PetscMPIInt    size;
    PetscReal      rtol = 1.0e-50;  // 相对残差阈值
    PetscReal      atol = 1.0e-50;  // 绝对残差阈值
    PetscReal      stol = 1.0e-8;  // 范数残差
    PetscInt       maxit = 100;  // 最大迭代次数
    PetscInt       maxf = 16; // 最大函数评估次数

    PetscInitialize(&argc,&args,(char*)0,help);
    MPI_Comm_size(PETSC_COMM_WORLD,&size);
    if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Example is only for sequential runs");

     /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create nonlinear solver context
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    SNESCreate(PETSC_COMM_WORLD,&snes);

    SNESSetTolerances(snes, atol, rtol, stol, maxit, maxf);
    SNESSetType(snes, SNESNEWTONLS);

    /* 设置 Jacobian 的更新策略为 Modified Newton-Raphson */
    /* -1 表示在整个求解过程中只计算一次 Jacobian */
    // SNESSetLagJacobian(snes, -2);

    SNESGetLineSearch(snes, &linesearch);
    SNESLineSearchSetType(linesearch, SNESLINESEARCHNONE);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Create matrix and vector data structures; set corresponding routines
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
        Create vectors for solution and nonlinear function
    */
    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,n_eq);
    VecSetFromOptions(x);
    VecDuplicate(x,&r);

    /*
        Create Jacobian matrix data structure
    */
    MatCreate(PETSC_COMM_WORLD,&J);
    MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n_eq,n_eq);
    MatSetFromOptions(J);
    MatSetUp(J);

    /*
    Set function evaluation routine and vector.
    */
    SNESSetFunction(snes,r,FormFunction,&user);

    // std::cout << "r: \n";
    // ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    /*
    Set Jacobian matrix data structure and Jacobian evaluation routine
    */
    SNESSetJacobian(snes,J,J,FormJacobian,&user);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Customize nonlinear solver; set runtime options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    Set linear solver defaults for this problem. By extracting the
    KSP and PC contexts from the SNES context, we can then
    directly call any KSP and PC routines to set various options.
    */
    PetscReal l_rtol = 1.0e-5;
    PetscReal l_atol = 1.0e-50;
    PetscReal l_dtol = 1.0e20;
    PetscInt l_maxits = 200;
    SNESGetKSP(snes,&ksp);
    KSPSetTolerances(ksp,l_rtol, l_atol, l_dtol, l_maxits);
    KSPSetFromOptions(ksp);

    /*
    Set SNES/KSP/KSP/PC runtime options, e.g.,
        -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
    These options will override those specified above as long as
    SNESSetFromOptions() is called _after_ any other customization
    routines.
    */
    
    SNESSetFromOptions(snes);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Evaluate initial guess; then solve nonlinear system
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    VecSet(x,0.0);
    
    /*
    Note: The user should initialize the vector, x, with the initial guess
    for the nonlinear solver prior to calling SNESSolve().  In particular,
    to employ an initial guess of zero, the user should explicitly set
    this vector to zero by calling VecSet().
    */
    // std::cout << "x: \n";
    // ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    SNESSolve(snes,NULL,x);

    // Print solution uh
    std::cout << "uh \n"; 
    for(int ii = 0; ii < user.n_np; ++ii){
        std::cout << std::setprecision(15) << user.uh[ii] << "\n";
    }
    std::cout << std::endl;


    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Free work space.  All PETSc objects should be destroyed when they
    are no longer needed.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    VecDestroy(&x); VecDestroy(&r);
    MatDestroy(&J); SNESDestroy(&snes);
    ierr = PetscFinalize();
    return ierr;

}

