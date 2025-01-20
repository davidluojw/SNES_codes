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
// PetscErrorCode FormFunction(SNES snes,Vec x,Vec f,void *ctx)
// {
//     AppCtx  *user = (AppCtx *)ctx;
//     PetscErrorCode    ierr;
//     const PetscScalar *xx;
//     PetscScalar       *F;

//     /*
//     Get pointers to vector data.
//         - For default PETSc vectors, VecGetArray() returns a pointer to
//             the data array.  Otherwise, the routine is implementation dependent.
//         - You MUST call VecRestoreArray() when you no longer need access to
//             the array.
//     */
//     ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);
//     ierr = VecGetArray(f,&F);CHKERRQ(ierr);

//     VecZeroEntries(f);

//     // AssemblyResidual(xx, *user);

//     // Assembly stiffness matrix and load vector
//     for(int ii = 0; ii < user->n_eq; ++ii){
//         user->uh[ii] = xx[ii];
//     }

//     for (int ee = 0; ee < user->nElem; ++ee) {
//         // Allocate zero element stiffness matrix 
//         // std::vector<double> k_ele(user->n_en * user->n_en, 0.0);
//         std::vector<double> f_ele(user->n_en, 0.0);

//         std::vector<double> x_ele(user->n_en);
//         std::vector<double> d_ele(user->n_en);
//         std::vector<double> u_ele(user->n_en);
//         for (int aa = 0; aa < user->n_en; ++aa) {
//             x_ele[aa] = user->x_coor[user->IEN[ee + aa * user->nElem] - 1];
//             // if (user->ID[user->IEN[ee + aa * user->nElem] - 1] == user->n_eq) u_ele[aa] = user->uh[user->n_np-1];
//             u_ele[aa] = user->uh[user->IEN[ee + aa * user->nElem] - 1];
//         }

//         // loop over quadrature points
//         for (int ll = 0; ll < user->nqp; ++ll) {
//             // we need the geometric mapping at each quadrature points
//             double x_qua = 0.0;
//             double dx_dxi = 0.0;
//             double u_qua = 0.0;
//             double du_dxi = 0.0;
//             for (int aa = 0; aa < user->n_en; ++aa) {
//                 double Na = 0.0;
//                 PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
//                 x_qua += x_ele[aa] * Na;
//                 u_qua += u_ele[aa] * Na;
//                 double Na_xi = 0.0;
//                 PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);
//                 dx_dxi += x_ele[aa] * Na_xi;
//                 du_dxi += u_ele[aa] * Na_xi;
//             }
//             double dxi_dx = 1.0 / dx_dxi;

//             double kappa = fun_kappa( u_qua );
//             double dkappa = fun_dkappa( u_qua );

//             double detJ = dx_dxi;
//             // we loop over a and b to assemble the element stiffness matrix and load vector
//             for (int aa = 0; aa < user->n_en; ++aa){
//                 double Na = 0.0;
//                 PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
//                 double Na_xi = 0.0;
//                 PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);

//                 f_ele[aa] += user->wq[ll] * Na * f_body(x_qua) * detJ;
//                 f_ele[aa] -= user->wq[ll] * Na_xi * kappa * du_dxi * dxi_dx;

//             }
//         }
//         // global assembly
//         // distribute the entries to the global stiffness matrix and global load vector
//         for (int aa = 0; aa < user->n_en; ++aa){
//             int PP = user->ID[ user->IEN[aa * user->nElem + ee] - 1 ];
//             if (PP > 0){
//                 F[PP - 1] += f_ele[aa];
//             }
//         }
//         if (ee == 0){
//             F[ user->ID[ user->IEN[0] - 1 ] - 1 ] += h( user->x_coor[ user->IEN[0] - 1 ]);
//         }
        
//     } //  end of element loop

//     for (int ii = 0; ii < user->n_eq; ++ii){
//         F[ii] = -F[ii];
//     }
    
//     // std::cout << "F: \n";
//     // for (int ii = 0; ii < user->n_eq; ++ii){
//     //     std::cout << F[ii] << "\n";
//     // }
//     // std::cout << std::endl;

//     // std::cout << "xx: \n";
//     // for (int ii = 0; ii < user->n_eq; ++ii){
//     //     std::cout << xx[ii] << "\n";
//     // }
//     // std::cout << std::endl;


//     /* Restore vectors */
//     ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);
//     ierr = VecRestoreArray(f,&F);CHKERRQ(ierr);


//     // ierr = VecView(f,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

//     return 0;
// }


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
// PetscErrorCode FormJacobian(SNES snes,Vec x,Mat jac,Mat B,void *ctx)
// {
//     AppCtx  *user = (AppCtx *)ctx;
//     const PetscScalar *xx;
//     PetscScalar       A[100];
//     PetscErrorCode    ierr;
//     //   PetscInt          idx[2] = {0,1};


//     /*
//         Get pointer to vector data
//     */
//     ierr = VecGetArrayRead(x,&xx);CHKERRQ(ierr);

//     MatZeroEntries(B);

//     std::vector<double> K(user->n_eq * user->n_eq, 0.0);

//     /*
//         Compute Jacobian entries and insert into matrix.
//         - Since this is such a small problem, we set all entries for
//             the matrix at once.
//     */
//     // AssemblyJacobian(xx, *user);
//     // Assembly stiffness matrix and load vector
//     for(int ii = 0; ii < user->n_eq; ++ii){
//         user->uh[ii] = xx[ii];
//     }

//     for (int ee = 0; ee < user->nElem; ++ee) {
//         // Allocate zero element stiffness matrix 
//         std::vector<double> k_ele(user->n_en * user->n_en, 0.0);
//         // std::vector<double> f_ele(user->n_en, 0.0);

//         std::vector<double> x_ele(user->n_en);
//         std::vector<double> d_ele(user->n_en);
//         std::vector<double> u_ele(user->n_en);
//         for (int aa = 0; aa < user->n_en; ++aa) {
//             x_ele[aa] = user->x_coor[user->IEN[ee + aa * user->nElem] - 1];
//             // if (user->ID[user->IEN[ee + aa * user->nElem] - 1] == user->n_eq) u_ele[aa] = user->uh[user->n_np-1];
//             u_ele[aa] = user->uh[user->IEN[ee + aa * user->nElem] - 1];
//         }

//         // loop over quadrature points
//         for (int ll = 0; ll < user->nqp; ++ll) {
//             // we need the geometric mapping at each quadrature points
//             double x_qua = 0.0;
//             double dx_dxi = 0.0;
//             double u_qua = 0.0;
//             double du_dxi = 0.0;
//             for (int aa = 0; aa < user->n_en; ++aa) {
//                 double Na = 0.0;
//                 PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
//                 x_qua += x_ele[aa] * Na;
//                 u_qua += u_ele[aa] * Na;
//                 double Na_xi = 0.0;
//                 PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);
//                 dx_dxi += x_ele[aa] * Na_xi;
//                 du_dxi += u_ele[aa] * Na_xi;
//             }
//             double dxi_dx = 1.0 / dx_dxi;

//             double kappa = fun_kappa( u_qua );
//             double dkappa = fun_dkappa( u_qua );

//             double detJ = dx_dxi;
//             // we loop over a and b to assemble the element stiffness matrix and load vector
//             for (int aa = 0; aa < user->n_en; ++aa){
//                 double Na = 0.0;
//                 PolyBasis(user->pp, aa + 1, 0, user->qp[ll], Na);
//                 double Na_xi = 0.0;
//                 PolyBasis(user->pp, aa + 1, 1, user->qp[ll], Na_xi);

//                 for (int bb = 0; bb < user->n_en; ++bb){
//                     double Nb = 0.0;
//                     PolyBasis(user->pp, bb + 1, 0, user->qp[ll], Nb);
//                     double Nb_xi = 0.0;
//                     PolyBasis(user->pp, bb + 1, 1, user->qp[ll], Nb_xi);
//                     k_ele[aa * user->n_en + bb] += user->wq[ll] * Na_xi * kappa * Nb_xi * dxi_dx;
//                     k_ele[aa * user->n_en + bb] += user->wq[ll] * Na_xi * dkappa * Nb * du_dxi * dxi_dx;
//                 }
//             }
//         }
//         // global assembly
//         // distribute the entries to the global stiffness matrix and global load vector
//         for (int aa = 0; aa < user->n_en; ++aa){
//             int PP = user->ID[ user->IEN[aa * user->nElem + ee] - 1 ];
//             if (PP > 0){
//                 for (int bb = 0; bb < user->n_en; ++bb){
//                     int QQ = user->ID[ user->IEN[bb * user->nElem + ee] - 1 ];
//                     if (QQ > 0){
//                         K[(PP - 1) * user->n_eq + QQ - 1] += k_ele[aa * user->n_en + bb];
//                     }
//                 }
//             }
//         }
        
//     } //  end of element loop

//     // std::cout << "K: ";
//     // for (int ii = 0; ii < user->n_eq * user->n_eq; ++ii){
//     //     if (ii % user->n_eq == 0) std::cout << std::endl;
//     //     std::cout << K[ii] << "\t";
//     // }
//     // std::cout << std::endl;

//     for (PetscInt ii = 0; ii < user->n_eq; ++ii){
//         for (PetscInt jj = 0; jj < user->n_eq; ++jj){
//             ierr = MatSetValues(B, 1, &ii, 1, &jj, &K[ii * user->n_eq + jj], INSERT_VALUES); CHKERRQ(ierr);
//         }
//     }
        

//     /*
//         Restore vector
//     */
//     ierr = VecRestoreArrayRead(x,&xx);CHKERRQ(ierr);

//     /*
//         Assemble matrix
//     */
//     ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//     ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//     // ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//     if (jac != B) {
//         ierr = MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//         ierr = MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
//     }
//     return 0;
// }



