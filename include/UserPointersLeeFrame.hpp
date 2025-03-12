#ifndef USERPOINTERSLEEFRAME_HPP
#define USERPOINTERSLEEFRAME_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"
#include "Eigen"
#include "LeeFrameTools.hpp"

class UserPointersLeeFrame
{
    public:
        //  model data
        double L;
        double E;
        double I;
        double A;
        double nu;

        double F_applied;
        double M;
        double Q;
        double g1;
        double g2;
        double g3;
        double g4;

        // arc-length
        double Deltalambda_init = 0.01;  // initial increment of lamdba
        int max_steps = 10000;           // maximum load steps
        int max_iter  = 1000;            // maximum iterate steps
        double max_ratio = 1.00;
        double tol       = 1e-05;        // tolerance of equilibrium

        // parameters of the FEM
        int n_sd;              // space dimension
        int n_el;             // number of elements
        int n_en;              // number of element nodes
        int n_ed;              // number of element degrees of freedom (per node)
        int deg;       // polynomial degree
        int n_np; // number of points
        int n_eq;       // number of equations
        int n_ee;           // number of equations of an element

        //  quadrature rule
        int n_int;              // number of quadrature points
        std::vector<double> xi;
        std::vector<double> weight;

        // mesh
        Eigen::VectorXd x_coor;
        Eigen::VectorXd y_coor;
        Eigen::MatrixXd ID;
        Eigen::MatrixXd IEN;
        Eigen::MatrixXd LM;

        // solution
        Eigen::VectorXd disp;
        Eigen::VectorXd disp_n; 
        Vec d;
        Vec d_n;
        Vec Deltad;
        Vec Deltad_n;
        Vec deltad;
        Vec q_bar;
        Vec q_bar_n;
        Vec R;
        Vec N;
        Vec Deltad_bar;
        double lambda;
        double lambda_n;
        double deltalambda;
        double Deltalambda;
        double Deltalambda_n;
        Results results;
        Eigen::VectorXd Force_ext;
        Vec F_ext;
        std::vector<Node> Nodes;
        std::vector<Element> Elems;
        std::vector<Element> Elems_n;
        Eigen::VectorXd deltadisp;


};

#endif
