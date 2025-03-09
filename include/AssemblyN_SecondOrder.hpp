#ifndef ASSEMBLYN_SECONDORDER_HPP
#define ASSEMBLYN_SECONDORDER_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"
#include "Eigen"

inline void AssemblyN_SecondOrder(Vec &Nout, void* ctx){

    UserPointersLeeFrame *user = (UserPointersLeeFrame *)ctx;
    // Allocate an empty stiffness matrix

    PetscScalar       *N;
    VecGetArray(Nout,&N);

    // Eigen::VectorXd N(user->n_eq);

    // 遍历所有单元
    for (int ee = 0; ee < user->n_el; ++ee) {
        
        // 计算单元长度变化
        const double L_1 = user->Elems_n[ee].eleLength;
        const double hh = elemLength(user->Elems_n[ee], user->disp);
        const double D_L = hh - L_1;
        
        // 计算刚体转动增量
        const double rbr = elemRigidRot(user->Elems_n[ee], user->disp, user->disp_n);
        
        // 更新单元绝对角度
        double angle = user->Elems_n[ee].eleAngle + rbr;
        // Update the angle 
        user->Elems[ee].eleAngle = angle;
        
        // 计算旋转矩阵
        Eigen::MatrixXd T = RotationMatrix(angle);
        
        // 计算位移增量
        const Eigen::VectorXd delta_disp = user->disp - user->disp_n;
        
        // 提取节点转动分量
        const double r1 = delta_disp.eval()[user->Elems_n[ee].Node1.dof(2)] - rbr; // 注意C++的0-based索引
        const double r2 = delta_disp.eval()[user->Elems_n[ee].Node2.dof(2)] - rbr;
        
        // 局部坐标系变形增量
        Eigen::VectorXd dl(6);
        dl << 0.0, 0.0, r1, D_L, 0.0, r2;
        
        // 计算内力增量
        Eigen::VectorXd D_fl = user->Elems_n[ee].eleElasticK * dl;
        
        // 更新总内力
        Eigen::VectorXd fl = user->Elems_n[ee].eleFint + D_fl;

        // Update internal force
        user->Elems[ee].eleFint = fl;
        
        // 转换到全局坐标系
        Eigen::VectorXd n_e = T * user->Elems[ee].eleFint;
        
        // 组装到全局载荷向量
        for (int aa = 0; aa < user->n_en; ++aa) {
            for (int ii = 0; ii < user->n_ed; ++ii) {
                const int pp = user->n_ed * aa + ii; // 0-based索引
                const int PP = user->LM.coeff(pp,ee); // 假设LM存储方式为[DOF][单元]
                
                if (PP >= 0) { // 假设有效自由度为非
                    N[PP] += n_e(pp);
                }
            }
        }
    }

    PetscScalar *F_ext_arr;
    VecGetArray(user->F_ext, &F_ext_arr);

    for (int ii = 0; ii < user->n_eq; ++ii){
        N[ii] = N[ii] - user->lambda * F_ext_arr[ii];
    }

    VecRestoreArray(Nout,&N);

}

#endif