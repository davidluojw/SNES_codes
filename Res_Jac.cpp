#include "Res_Jac.hpp"

PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void* ctx)
{
    UserPointers *user = (UserPointers *) ctx;
    const PetscScalar *xx;
    PetscScalar       *F;

    /*
    Get pointers to vector data.s
        - For default PETSc vectors, VecGetArray() returns a pointer to
            the data array.  Otherwise, the routine is implementation dependent.
        - You MUST call VecRestoreArray() when you no longer need access to
            the array.
    */
    VecGetArrayRead(x,&xx);
    VecGetArray(f,&F);

    VecZeroEntries(f);

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

    
    /* Restore vectors */
    VecRestoreArrayRead(x,&xx);
    VecRestoreArray(f,&F);


    // VecView(f,PETSC_VIEWER_STDOUT_WORLD);

    return 0;


}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void* ctx)
{
    UserPointers *user = (UserPointers *)ctx;
    const PetscScalar *xx;
    /*
        Get pointer to vector data
    */
    VecGetArrayRead(x,&xx);

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

            // double detJ = dx_dxi;
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

    for (PetscInt ii = 0; ii < user->n_eq; ++ii){
        for (PetscInt jj = 0; jj < user->n_eq; ++jj){
            MatSetValues(B, 1, &ii, 1, &jj, &K[ii * user->n_eq + jj], INSERT_VALUES);
        }
    }

    // Print solution uh
    std::cout << "uh \n"; 
    for(int ii = 0; ii < user->n_np; ++ii){
        std::cout << std::setprecision(15) << user->uh[ii] << "\n";
    }
    std::cout << std::endl;
        
    /*
        Restore vector
    */
    VecRestoreArrayRead(x,&xx);

    /*
        Assemble matrix
    */
    MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
    // MatView(B,PETSC_VIEWER_STDOUT_WORLD);
    if (jac != B) {
        MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
    }
    return 0;
}