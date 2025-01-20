#ifndef USERPOINTERS_HPP
#define USERPOINTERS_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"

class UserPointers
{
    public:
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

        double EPS;

};

#endif
