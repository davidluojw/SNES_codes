#ifndef ASSEMBLYK_SECONDORDER_HPP
#define ASSEMBLYK_SECONDORDER_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"
#include "Eigen"

inline void AssemblyK_SecondOrder(Mat &Kout, void* ctx){

    UserPointersLeeFrame *user = (UserPointersLeeFrame *)ctx;
    // Allocate an empty stiffness matrix

    Eigen::MatrixXd K(user->n_eq, user->n_eq);

    // Assembly the stiffness matrix and load vector
    for (int ee = 0; ee < user->n_el; ee++) {
        // Element length
        double L_1 = user->Elems_n[ee].eleLength; // Initial load step element length
        double hh = elemLength(user->Elems_n[ee], user->disp); // Current element length
        double D_L = hh - L_1; // Element axial deformation

        // Rigid rotation increment
        double rbr = elemRigidRot(user->Elems_n[ee], user->disp, user->disp_n);
        // Absolute angle of the element
        double angle = user->Elems_n[ee].eleAngle + rbr;

        // Update the angle
        user->Elems[ee].eleAngle = angle;

        // Compute the rotation matrix
        Eigen::MatrixXd T = RotationMatrix(angle);

        // Deformation rotation, w.r.t initial load step
        Eigen::VectorXd deltadisp = user->disp - user->disp_n;
        double r1 = deltadisp.eval()[user->Elems_n[ee].Node1.dof(2)] - rbr;
        double r2 = deltadisp.eval()[user->Elems_n[ee].Node2.dof(2)] - rbr;


        // Local coordinate, deformation increment, w.r.t initial load step
        Eigen::VectorXd dl(6);
        dl << 0, 0, r1, D_L, 0, r2;

        // Current iterate internal forces
        double P = user->Elems[ee].eleFint(4);
        double M1 = user->Elems[ee].eleFint(3);
        double M2 = user->Elems[ee].eleFint(6);

        // Local coordinate, displacement & angle increment with respect to the
        Eigen::VectorXd x_ele(2);
        x_ele << 0, hh;

        // Allocate zero element stiffness matrix
        Eigen::MatrixXd ke_ele = Eigen::MatrixXd::Zero(user->n_ee, user->n_ee); // Dimension of element stiffness is n_ee x n_ee
        Eigen::MatrixXd kg_ele = Eigen::MatrixXd::Zero(user->n_ee, user->n_ee);

        // Loop over quadrature points
        for (int ll = 0; ll < user->n_int; ll++) {
            double x_l = 0;
            double dx_dxi = 0;
            double dv_dxi = 0;
            double du_dxi = 0;

            for (int aa = 0; aa < user->n_en; aa++) {
                double Na = 0.0;
                double Na_xi = 0.0;
                PolyBasis(user->deg, aa+1, 0, user->xi[ll], Na);
                PolyBasis(user->deg, aa+1, 1, user->xi[ll], Na_xi);
                x_l += x_ele(aa) * Na;
                dx_dxi += x_ele(aa) * Na_xi;
            }

            double dxi_dx = 1.0 / dx_dxi;

            for (int aa = 0; aa < user->n_en; aa++) {
                for (int ii = 0; ii < user->n_ed; ii++) {
                    int pp = user->n_ed * (aa) + ii;

                    for (int bb = 0; bb < user->n_en; bb++) {
                        for (int jj = 0; jj < user->n_ed; jj++) {
                            int qq = user->n_ed * (bb) + jj;

                            if (ii == 0 && jj == 0) {
                                double Na_xi = 0.0;
                                double Nb_xi = 0.0;
                                PolyBasis(user->deg, aa+1, 1, user->xi[ll], Na_xi);
                                PolyBasis(user->deg, bb+1, 1, user->xi[ll], Nb_xi);
                                ke_ele(pp, qq) = ke_ele(pp, qq) + user->E * user->A * user->weight[ll] * Na_xi * Nb_xi * dxi_dx;
                                kg_ele(pp, qq) = kg_ele(pp, qq) + P * user->weight[ll] * Na_xi * Nb_xi * dxi_dx;
                            }

                            if (ii != 0 && jj != 0) {
                                int pp_ind = (user->n_ed - 1) * (aa) + ii ;
                                int qq_ind = (user->n_ed - 1) * (bb) + jj ;
                                double Np_xixi = 0.0;
                                double Nq_xixi = 0.0;
                                HermiteShape(pp_ind, user->xi[ll], 2, hh, Np_xixi);
                                HermiteShape(qq_ind, user->xi[ll], 2, hh, Nq_xixi);
                                ke_ele(pp, qq) = ke_ele(pp, qq) + user->E * user->I * user->weight[ll] * Np_xixi * Nq_xixi * dxi_dx * dxi_dx * dxi_dx;
                                kg_ele(pp, qq) = kg_ele(pp, qq) + P * user->weight[ll] * Np_xixi * Nq_xixi * dxi_dx;
                            }
                        }
                    }
                }
            }
        } // End of quadrature loop

        // Generate the tangent matrix = elastic matrix + geometric matrix
        Eigen::MatrixXd kt_ele = ke_ele + kg_ele;

        // Store elastic stiffness matrix to be used in computation of internal forces
        user->Elems[ee].eleElasticK = ke_ele;
        user->Elems[ee].eleTangentK = kt_ele;

        // Transform the stiffness matrix in terms of the global coordinates
        kt_ele = T * kt_ele * T.transpose();

        // Now we need to put element k and f into global K and F
        for (int aa = 0; aa < user->n_en; aa++) {
            for (int ii = 0; ii < user->n_ed; ii++) {
                int pp = user->n_ed * (aa) + ii;
                int PP = user->LM.coeff(pp, ee);
                if (PP >= 0) {
                    for (int bb = 0; bb < user->n_en; bb++) {
                        for (int jj = 0; jj < user->n_ed; jj++) {
                            int qq = user->n_ed * (bb) + jj;
                            int QQ = user->LM.coeff(qq, ee);
                            if (QQ >= 0) {
                                K.coeffRef(PP, QQ) += kt_ele.coeff(pp, qq);
                            }
                        }
                    }
                }
            }
        }
    }

    for (PetscInt ii = 0; ii < user->n_eq; ++ii){
        for (PetscInt jj = 0; jj < user->n_eq; ++jj){
            MatSetValues(Kout, 1, &ii, 1, &jj, &K.coeffRef(ii, jj), INSERT_VALUES);
        }
    }

    MatAssemblyBegin(Kout, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(Kout, MAT_FINAL_ASSEMBLY);
}

#endif