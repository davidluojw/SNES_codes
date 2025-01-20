#ifndef RES_JAC_HPP
#define RES_JAC_HPP

#include "UserPointers.hpp"
#include "PolyBasis.hpp"
#include "ProblemSettings.hpp"
#include <iomanip> 

extern PetscErrorCode FormJacobian(SNES,Vec,Mat,Mat,void*);
extern PetscErrorCode FormFunction(SNES,Vec,Vec,void*);


#endif