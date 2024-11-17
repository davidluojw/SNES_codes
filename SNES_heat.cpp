#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"

// using namespace std;

constexpr double EPS = 1e-14;

double f(double xx, double yy) {
    double dist = sqrt( (xx - 0.5) * (xx - 0.5) + (yy - 0.5) *  (yy - 0.5));
    if (dist < 0.5) return 0.5;
    else return 0.0;
}

double g(double xx, double yy, double tt){
    if (tt < 1.0) return tt;
    else return 1.0;
}

double g_dot(double xx, double yy, double tt){
    if (tt < 1.0) return 1.0;
    else return 0.0;
}

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

void Quad_grad(int aa, double xi, double eta, double &dval_dxi, double &dval_eta) {
    if (aa == 1) {
        dval_dxi = -0.25 * (1 - eta);
        dval_eta = -0.25 * (1 - xi);
    }
    if (aa == 2) {
        dval_dxi = 0.25 * (1 - eta);
        dval_eta = -0.25 * (1 + xi);
    }
    if (aa == 3) {
        dval_dxi = 0.25 * (1 + eta);
        dval_eta = 0.25 * (1 + xi);
    }
    if (aa == 4) {
        dval_dxi = -0.25 * (1 + eta);
        dval_eta = 0.25 * (1 - xi);
    }
}

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


// manufactured solution and source term
double exact(double x, double y) {
    return x * (1 - x) * y * (1 - y);
}
double exact_x(double x, double y) {
    return (1 - 2 * x) * y * (1 - y);
}
double exact_y(double x, double y) {
    return x * (1 - x) * (1 - 2 * y);
}


int main() {
    double kappa = 1.0; // isotropic homogeneous heat conductivity
    double rho = 1.0; // density
    double cap = 1.0; // heat capacity

    double alpha = 0.5;

    double T_final = 10.0;
    double dt = 0.1;
    std::vector<double> t(1 + T_final / dt, 0.0);  // time sub-interval
    for (int i = 0; i <= T_final / dt; ++i) {
        t[i] = i * dt;
    }
    double NN = T_final / dt;

    // quadrature rule
    int n_int_xi = 3; // number of quadrature points in xi-direction
    int n_int_eta = 3; // number of quadrature points in eta-direction
    int n_int = n_int_xi * n_int_eta; 
    std::vector<double> xi(n_int, 0.0);
    std::vector<double> eta(n_int, 0.0);
    std::vector<double> weight(n_int, 0.0);
    Gauss2D(n_int_xi, n_int_xi, xi, eta, weight);

    // FEM mesh settings
    int n_en = 4;      // 4-node quadrilateral element
    int n_ed = 1;      // number of element degrees of freedom

    int n_el_x = 7; // number of element in x-direction
    int n_el_y = 6; // number of element in y-direction
    int n_el = n_el_x * n_el_y; // total number of element in 2D domain

    int n_np_x = n_el_x + 1; // number of node points in x-direction
    int n_np_y = n_el_y + 1; // number of node points in y-direction
    int n_np = n_np_x * n_np_y; // total number of node points in 2D domain

    // generate the coordinates of the nodal points
    std::vector<double> x_coor(n_np, 0.0), y_coor(n_np, 0.0);
    double hh_x = 1.0 / n_el_x; // mesh size in the x-direction
    double hh_y = 1.0 / n_el_y; // mesh size in the y-direction
    for (int ny = 0; ny < n_np_y; ++ny) {
        for (int nx = 0; nx < n_np_x; ++nx) {
            int index = ny * n_np_x + nx;
            x_coor[index] = nx * hh_x;
            y_coor[index] = ny * hh_y;
        }
    }

    // setup the IEN array for element with local node numbering as
    // a=4 ------- a=3
    // |           |
    // |           |
    // |           |
    // |           |
    // a=1 ------- a=2
    std::vector<int> IEN(n_el * n_en, 0);
    for (int ey = 0; ey < n_el_y; ++ey) {
        for (int ex = 0; ex < n_el_x; ++ex) {
            int ee = ey * n_el_x + ex;
            IEN[ee + 0*n_el] = ey * n_np_x + ex + 1;
            IEN[ee + 1*n_el] = ey * n_np_x + ex + 2;
            IEN[ee + 2*n_el] = (ey + 1) * n_np_x + ex + 2;
            IEN[ee + 3*n_el] = (ey + 1) * n_np_x + ex + 1;
        }
    }

    // ID array
    std::vector<int> ID(n_np, 0);
    int counter = 1;
    for (int ny = 1; ny < n_np_y - 1; ++ny) {
        for (int nx = 1; nx < n_np_x - 1; ++nx) {
            ID[ny * n_np_x + nx] = counter;
            counter = counter + 1;
        }
    }

    int n_eq = n_np - n_np_x * 2 - n_np_y * 2 + 4;

    std::vector<int> LM(n_el * n_en, 0);
    for (int ii = 0; ii < n_el * n_en; ++ii) {
        LM[ii] = ID[IEN[ii] - 1];
    }

    // Start the assembly procedure
    std::vector<double> K(n_eq * n_eq, 0.0);
    std::vector<double> M(n_eq * n_eq, 0.0);
    std::vector<double> F(n_eq, 0.0);

    for (int ee = 0; ee < n_el; ++ee) {
        std::vector<double> k_ele(n_en * n_en, 0.0);
        std::vector<double> m_ele(n_en * n_en, 0.0);
        std::vector<double> f_ele(n_en, 0.0);

        std::vector<double> x_ele(n_en);
        std::vector<double> y_ele(n_en);
        for (int ii = 0; ii < n_en; ++ii) {
            x_ele[ii] = x_coor[IEN[ee + ii * n_el] - 1];
            y_ele[ii] = y_coor[IEN[ee + ii * n_el] - 1];
        }

        // loop over quadrature points
        for (int ll = 0; ll < n_int; ++ll) {
            // we need the geometric mapping at each quadrature points
            double x_l = 0.0, y_l = 0.0;
            double dx_dxi = 0.0, dx_deta = 0.0;
            double dy_dxi = 0.0, dy_deta = 0.0;
            for (int aa = 0; aa < n_en; ++aa) {
                double Na = 0.0;
                Quad(aa + 1, xi[ll], eta[ll], Na);
                x_l += x_ele[aa] * Na;
                y_l += y_ele[aa] * Na;
                double Na_xi = 0.0, Na_eta = 0.0;
                Quad_grad(aa + 1, xi[ll], eta[ll], Na_xi, Na_eta);
                dx_dxi += x_ele[aa] * Na_xi;
                dx_deta += x_ele[aa] * Na_eta;
                dy_dxi += y_ele[aa] * Na_xi;
                dy_deta += y_ele[aa] * Na_eta;

            }

            double detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;
            // we loop over a and b to assemble the element stiffness matrix and load vector
            for (int aa = 0; aa < n_en; ++aa){
                double Na = 0.0;
                Quad(aa + 1, xi[ll], eta[ll], Na);
                double Na_xi = 0.0, Na_eta = 0.0;
                Quad_grad(aa + 1, xi[ll], eta[ll], Na_xi, Na_eta);

                // See page 147 of Sec. 3.9 of the textbook for the shape function routine details
                double Na_x = (Na_xi * dy_deta - Na_eta * dy_dxi) / detJ;
                double Na_y = (Na_xi * (-dx_deta) + Na_eta * dx_dxi)  / detJ;

                f_ele[aa] += weight[ll] * detJ * f(x_l, y_l) * Na;

                for (int bb = 0; bb < n_en; ++bb){
                    double Nb = 0.0;
                    Quad(bb + 1, xi[ll], eta[ll], Nb);
                    double Nb_xi = 0.0, Nb_eta = 0.0;
                    Quad_grad(bb + 1, xi[ll], eta[ll], Nb_xi, Nb_eta);

                    double Nb_x = (Nb_xi * dy_deta - Nb_eta * dy_dxi) / detJ;
                    double Nb_y = (Nb_xi * (-dx_deta) + Nb_eta * dx_dxi)  / detJ;

                    // See page 421
                    m_ele[aa * n_en + bb] += weight[ll] * detJ * rho * cap * Na * Nb;
                    k_ele[aa * n_en + bb] += weight[ll] * detJ * kappa * ( Na_x * Nb_x + Na_y * Nb_y);
                }
            }
        }
        // global assembly
        for (int aa = 0; aa < n_en; ++aa){
            int PP = LM[aa * n_el + ee];
            if (PP > 0){
                F[PP - 1] += f_ele[aa];
                for (int bb = 0; bb < n_en; ++bb){
                    int QQ = LM[bb * n_el + ee];
                    if (QQ > 0){
                        M[(PP - 1) * n_eq + QQ - 1] += m_ele[aa * n_en + bb];
                        K[(PP - 1) * n_eq + QQ - 1] += k_ele[aa * n_en + bb];
                    }
                    else{
                        F[PP - 1] -= k_ele[aa * n_en + bb] * g(x_ele[bb], y_ele[bb], 0.0) + m_ele[aa * n_en + bb] * g_dot(x_ele[bb], y_ele[bb], 0.0);
                    }
                }
            }
        }
    } //  end of element loop

    // set initial condition
    Vec         dn, b, vn; /* approx solution, RHS, approx velo */
    Mat         Kstiff, Mass;       /* linear system matrix */
    KSP         ksp;     /* linear solver context */
    PC          pc;      /* preconditioner context */
    PetscReal   norm;    /* norm of solution error */
    PetscInt    i, n = n_eq, col[3], its;
    PetscMPIInt size;
    PetscScalar value[3];

    PetscFunctionBeginUser;
    PetscCall(PetscInitialize(&argc, &args, NULL, help));
    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCheck(size == 1, PETSC_COMM_WORLD, PETSC_ERR_WRONG_MPI_SIZE, "This is a uniprocessor example only!");

    PetscCall(PetscOptionsGetInt(NULL, NULL, "-n", &n, NULL));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            Compute the matrix and right-hand-side vector that define
            the linear system, Ax = b.
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    /*
        Create vectors.  Note that we form 1 vector from scratch and
        then duplicate as needed.
    */
    PetscCall(VecCreate(PETSC_COMM_SELF, &dn));
    PetscCall(PetscObjectSetName((PetscObject)dn, "Solution"));
    PetscCall(VecSetSizes(dn, PETSC_DECIDE, n));
    PetscCall(VecSetFromOptions(dn));
    PetscCall(VecDuplicate(dn, &b));
    PetscCall(VecDuplicate(dn, &vn));

    /*
        Create matrix.  When using MatCreate(), the matrix format can
        be specified at runtime.

        Performance tuning note:  For problems of substantial size,
        preallocation of matrix memory is crucial for attaining good
        performance. See the matrix chapter of the users manual for details.
    */
    PetscCall(MatCreate(PETSC_COMM_SELF, &Kstiff));
    PetscCall(MatSetSizes(Kstiff, PETSC_DECIDE, PETSC_DECIDE, n, n));
    PetscCall(MatSetFromOptions(Kstiff));
    PetscCall(MatSetUp(KStiff));
    PetscCall(MatDuplicate(Kstiff, &Mass));

    /*
        Assemble matrix
    */
    value[0] = -1.0;
    value[1] = 2.0;
    value[2] = -1.0;
    for (i = 1; i < n - 1; i++) {
        col[0] = i - 1;
        col[1] = i;
        col[2] = i + 1;
        PetscCall(MatSetValues(A, 1, &i, 3, col, value, INSERT_VALUES));
    }
    i      = n - 1;
    col[0] = n - 2;
    col[1] = n - 1;
    PetscCall(MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES));
    i        = 0;
    col[0]   = 0;
    col[1]   = 1;
    value[0] = 2.0;
    value[1] = -1.0;
    PetscCall(MatSetValues(A, 1, &i, 2, col, value, INSERT_VALUES));
    PetscCall(MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));

    /*
        Set exact solution; then compute right-hand-side vector.
    */
    PetscCall(VecSet(u, 1.0));
    PetscCall(MatMult(A, u, b));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Create the linear solver and set various options
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));

    /*
        Set operators. Here the matrix that defines the linear system
        also serves as the matrix that defines the preconditioner.
    */
    PetscCall(KSPSetOperators(ksp, A, A));

    /*
        Set linear solver defaults for this problem (optional).
        - By extracting the KSP and PC contexts from the KSP context,
        we can then directly call any KSP and PC routines to set
        various options.
        - The following four statements are optional; all of these
        parameters could alternatively be specified at runtime via
        KSPSetFromOptions();
    */
    if (!PCMPIServerActive) { /* cannot directly set KSP/PC options when using the MPI linear solver server */
        PetscCall(KSPGetPC(ksp, &pc));
        PetscCall(PCSetType(pc, PCJACOBI));
        PetscCall(KSPSetTolerances(ksp, 1.e-5, PETSC_CURRENT, PETSC_CURRENT, PETSC_CURRENT));
    }

    /*
        Set runtime options, e.g.,
            -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
        These options will override those specified above as long as
        KSPSetFromOptions() is called _after_ any other customization
        routines.
    */
    PetscCall(KSPSetFromOptions(ksp));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Solve the linear system
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(KSPSolve(ksp, b, x));

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                        Check the solution and clean up
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    PetscCall(VecAXPY(x, -1.0, u));
    PetscCall(VecNorm(x, NORM_2, &norm));
    PetscCall(KSPGetIterationNumber(ksp, &its));
    PetscCall(PetscPrintf(PETSC_COMM_SELF, "Norm of error %g, Iterations %" PetscInt_FMT "\n", (double)norm, its));

    /* check that KSP automatically handles the fact that the new non-zero values in the matrix are propagated to the KSP solver */
    PetscCall(MatShift(A, 2.0));
    PetscCall(KSPSolve(ksp, b, x));

    /*
        Free work space.  All PETSc objects should be destroyed when they
        are no longer needed.
    */
    PetscCall(KSPDestroy(&ksp));

    /* test if prefixes properly propagate to PCMPI objects */
    if (PCMPIServerActive) {
        PetscCall(KSPCreate(PETSC_COMM_SELF, &ksp));
        PetscCall(KSPSetOptionsPrefix(ksp, "prefix_test_"));
        PetscCall(MatSetOptionsPrefix(A, "prefix_test_"));
        PetscCall(KSPSetOperators(ksp, A, A));
        PetscCall(KSPSetFromOptions(ksp));
        PetscCall(KSPSolve(ksp, b, x));
        PetscCall(KSPDestroy(&ksp));
    }

    PetscCall(VecDestroy(&x));
    PetscCall(VecDestroy(&u));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));

    /*
        Always call PetscFinalize() before exiting a program.  This routine
        - finalizes the PETSc libraries as well as MPI
        - provides summary and diagnostic information if certain runtime
            options are chosen (e.g., -log_view).
    */
    PetscCall(PetscFinalize());







    return 0;
}
