#ifndef HERMITESHAPE_HPP
#define HERMITESHAPE_HPP

#include <iostream>
#include <cstdio>

inline void HermiteShape(int a, double xi, int der, double h, double &result) {
    if (a == 1) {
        if (der == 0) {
            result = 0.25 * (xi - 1) * (xi - 1) * (xi + 2);
        } else if (der == 1) {
            result = 0.75 * (xi * xi - 1);
        } else if (der == 2) {
            result = 1.5 * xi;
        } else if (der == 3) {
            result = 1.5;
        }
    } else if (a == 2) {
        if (der == 0) {
            result = 0.125 * h * (xi + 1) * (xi - 1) * (xi - 1);
        } else if (der == 1) {
            result = 0.125 * h * (xi - 1) * (3 * xi + 1);
        } else if (der == 2) {
            result = 0.25 * h * (3 * xi - 1);
        } else if (der == 3) {
            result = 0.25 * h * 3;
        }
    } else if (a == 3) {
        if (der == 0) {
            result = 0.25 * (xi + 1) * (xi + 1) * (2 - xi);
        } else if (der == 1) {
            result = 0.75 * (1 - xi * xi);
        } else if (der == 2) {
            result = -1.5 * xi;
        } else if (der == 3) {
            result = -1.5;
        }
    } else if (a == 4) {
        if (der == 0) {
            result = 0.125 * h * (xi + 1) * (xi + 1) * (xi - 1);
        } else if (der == 1) {
            result = 0.125 * h * (xi + 1) * (3 * xi - 1);
        } else if (der == 2) {
            result = 0.25 * h * (3 * xi + 1);
        } else if (der == 3) {
            result = 0.25 * h * 3;
        }
    }
}



#endif