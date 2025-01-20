#ifndef PROBLEMSETTINGS_HPP
#define PROBLEMSETTINGS_HPP

#include <iostream>
#include <cmath>
#include <vector>


// exact solution 
inline double exact (double x){
    return sin(x);
}

inline double exact_x(double x){
    return cos(x);
}

// nonlinear kappa
inline double fun_kappa(double x){
    return 1.0  + x * x;
}

inline double fun_dkappa(double x){
    return 2*x;
}

inline double f_body(double x){
    return -2.0*cos(x)*cos(x)*sin(x) + sin(x)*(sin(x)*sin(x) + 1.0);
}

inline double h(double x){
    return -1.0 * fun_kappa( exact(x) ) * exact_x(x);
}

inline double g(double x){
    return exact(x);
}

#endif