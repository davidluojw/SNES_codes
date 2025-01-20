#ifndef POLYBASIS_HPP
#define POLYBASIS_HPP
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

inline void PolyBasis(int degree , int i , int der , double x, double &poly){
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

#endif