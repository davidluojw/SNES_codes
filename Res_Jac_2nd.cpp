#include "Res_Jac2nd.hpp"

// PetscErrorCode FormFunction2nd(SNES snes, Vec x, Vec f, void* ctx)
// {
//     UserPointersLeeFrame *user = (UserPointersLeeFrame *) ctx;
//     const PetscScalar       *xx;
//     PetscScalar       *R;

//     PetscPrintf(PETSC_COMM_WORLD, "Assembly residual: \n");
//     /*
//     Get pointers to vector data.s
//         - For default PETSc vectors, VecGetArray() returns a pointer to
//             the data array.  Otherwise, the routine is implementation dependent.
//         - You MUST call VecRestoreArray() when you no longer need access to
//             the array.
//     */
//     VecGetArrayRead(x,&xx);
//     VecGetArray(f,&R);

//     std::cout << "f: ";
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscScalar value;
//         VecGetValues(f, 1, &i, &value);
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     PetscInt it;
//     SNESGetIterationNumber(snes, &it);
//     PetscPrintf(PETSC_COMM_WORLD, "iteration number: %d\n", it);

//     VecZeroEntries(f);

    
//     PetscPrintf(PETSC_COMM_WORLD, "Vector f (horizontal view): ");
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscScalar value;
//         VecGetValues(f, 1, &i, &value);
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     PetscPrintf(PETSC_COMM_WORLD, "Vector x (horizontal view): ");
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscScalar value;
//         VecGetValues(x, 1, &i, &value);
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     // if (it >= 0){
//     //     Vec sol;
//     //     SNESGetSolution(snes, &sol);
//     //     user->d = sol;
//     //     PetscPrintf(PETSC_COMM_WORLD, "Corrector solution captured\n");
//     //     PetscPrintf(PETSC_COMM_WORLD, "Vector d (horizontal view): ");
//     //     for (PetscInt i = 0; i < user->n_eq; ++i) {
//     //         PetscScalar value;
//     //         VecGetValues(user->d, 1, &i, &value);
//     //         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     //     }
//     //     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     //     //  get lambda
//     //     // PetscReal   lambda; 
//     //     // SNESNewtonALGetLoadParameter(snes, &lambda);
//     //     // user->lambda = lambda;
//     //     // PetscPrintf(PETSC_COMM_WORLD, "lambda: ");
//     //     // PetscPrintf(PETSC_COMM_WORLD, "%f ", user->lambda);
//     //     // PetscPrintf(PETSC_COMM_WORLD, "\n");


//     //     const PetscScalar *d_array;
//     //     PetscInt d_size;
//     //     VecGetSize(user->d, &d_size);
//     //     VecGetArrayRead(user->d, &d_array);
//     //     for (PetscInt ii = 0; ii < user->n_np * user->n_ed - 1; ++ii){
//     //         user->disp[ii + 2] = d_array[ii];
//     //     }
//     //     user->disp[user->n_np * user->n_ed - 1] = d_array[d_size - 1];
//     //     VecRestoreArrayRead(user->d, &d_array);
//     //     std::cout << "disp: " << user->disp.transpose() << std::endl;

//     //     // 计算增量
//     //     // user.deltad = user.d - user.d_n;
//     //     VecCopy(user->d, user->deltad);         // Copy user.d to user.deltad
//     //     VecAXPY(user->deltad, -1.0, user->d_n);  // Subtract user.d_n from user.deltad,  deltad = -1.0 d_n + deltad

//     //     user->deltalambda = user->lambda - user->lambda_n;
//     //     user->deltadisp = user->disp - user->disp_n;

//     //     VecDestroy(&sol);
//     // }
    
    
//     // 遍历所有单元
//     for (int ee = 0; ee < user->n_el; ++ee) {
        
//         // 计算单元长度变化
//         const double L_1 = user->Elems_n[ee].eleLength;
//         const double hh = elemLength(user->Elems_n[ee], user->disp);
//         const double D_L = hh - L_1;
        
//         // 计算刚体转动增量
//         const double rbr = elemRigidRot(user->Elems_n[ee], user->disp, user->disp_n);
        
//         // 更新单元绝对角度
//         double angle = user->Elems_n[ee].eleAngle + rbr;
//         // Update the angle 
//         user->Elems[ee].eleAngle = angle;
        
//         // 计算旋转矩阵
//         Eigen::MatrixXd T = RotationMatrix(angle);
        
//         // 计算位移增量
//         const Eigen::VectorXd delta_disp = user->disp - user->disp_n;
        
//         // 提取节点转动分量
//         const double r1 = delta_disp.eval()[user->Elems_n[ee].Node1.dof(2)] - rbr; // 注意C++的0-based索引
//         const double r2 = delta_disp.eval()[user->Elems_n[ee].Node2.dof(2)] - rbr;
        
//         // 局部坐标系变形增量
//         Eigen::VectorXd dl(6);
//         dl << 0.0, 0.0, r1, D_L, 0.0, r2;
        
//         // 计算内力增量
//         Eigen::VectorXd D_fl = user->Elems_n[ee].eleElasticK * dl;
        
//         // 更新总内力
//         Eigen::VectorXd fl = user->Elems_n[ee].eleFint + D_fl;

//         // Update internal force
//         user->Elems[ee].eleFint = fl;
        
//         // 转换到全局坐标系
//         Eigen::VectorXd n_e = T * user->Elems[ee].eleFint;
        
//         // 组装到全局载荷向量
//         for (int aa = 0; aa < user->n_en; ++aa) {
//             for (int ii = 0; ii < user->n_ed; ++ii) {
//                 const int pp = user->n_ed * aa + ii; // 0-based索引
//                 const int PP = user->LM.coeff(pp,ee); // 假设LM存储方式为[DOF][单元]
                
//                 if (PP >= 0) { // 假设有效自由度为非负
//                     R[PP] += n_e(pp);
//                 }
//             }
//         }
//     }
//     std::cout << "N: ";
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", R[i]);
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     PetscScalar *F_ext_arr;
//     VecGetArray(user->F_ext, &F_ext_arr);

//     std::cout << "F_ext: ";
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscScalar value;
//         VecGetValues(user->F_ext, 1, &i, &value);
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     for (int ii = 0; ii < user->n_eq; ++ii){
//         if (it == 0) continue;
//         else R[ii] = R[ii] - user->lambda * F_ext_arr[ii];
//     }

//     VecRestoreArrayRead(x,&xx);
//     VecRestoreArray(f,&R);

//     PetscReal FunNorm;
//     SNESGetFunctionNorm(snes, &FunNorm);
//     std::cout << "FunNorm: " << FunNorm << std::endl;

//     return 0;


// }

// PetscErrorCode FormJacobian2nd(SNES snes, Vec x, Mat jac, Mat B, void* ctx)
// {
//     UserPointersLeeFrame *user = (UserPointersLeeFrame *)ctx;
//     const PetscScalar *xx;

    
//     /*
//         Get pointer to vector data
//     */
//     VecGetArrayRead(x,&xx);

//     MatZeroEntries(B);

//     PetscPrintf(PETSC_COMM_WORLD, "Assembly Jacobian: \n");

//     AssemblyK_SecondOrder(B, user);

//     VecRestoreArrayRead(x,&xx);

//     /*
//         Assemble matrix
//     */
//     MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);
//     MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);
//     // MatView(B,PETSC_VIEWER_STDOUT_WORLD);
//     if (jac != B) {
//         MatAssemblyBegin(jac,MAT_FINAL_ASSEMBLY);
//         MatAssemblyEnd(jac,MAT_FINAL_ASSEMBLY);
//     }
//     return 0;
// }


// // 弧长约束方程（单独定义）
// PetscErrorCode ComputeTangentLoad(SNES snes, Vec X, Vec F, void *ctx) {
//     UserPointersLeeFrame *user = (UserPointersLeeFrame*)ctx;
//     const PetscScalar *x;
//     PetscScalar *f;

//     PetscFunctionBeginUser;
//     VecGetArrayRead(X, &x);

//     PetscPrintf(PETSC_COMM_WORLD, "ComputeTangentLoad captured\n");
//     PetscInt it;
//     SNESGetIterationNumber(snes, &it);
//     PetscPrintf(PETSC_COMM_WORLD, "iteration number: %d\n", it);

//     //  get lambda
//     PetscReal   lambda; 
//     SNESNewtonALGetLoadParameter(snes, &lambda);
//     if (it == 0) lambda = user->lambda;
//     else user->lambda = lambda;
//     PetscPrintf(PETSC_COMM_WORLD, "lambda: ");
//     PetscPrintf(PETSC_COMM_WORLD, "%f ", user->lambda);
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

//     VecZeroEntries(F);
//     VecCopy(user->F_ext, F);
//     VecScale(F, user->lambda);
//     VecRestoreArrayRead(X, &x);

//     PetscPrintf(PETSC_COMM_WORLD, "ComputeTangentLoad: Vector F (horizontal view): ");
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscScalar value;
//         VecGetValues(F, 1, &i, &value);
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

    


//     // update variables
//     Vec sol;
//     SNESGetSolution(snes, &sol);
//     user->d = sol;
//     PetscPrintf(PETSC_COMM_WORLD, "Update solution captured\n");
//     PetscPrintf(PETSC_COMM_WORLD, "Vector d (horizontal view): ");
//     for (PetscInt i = 0; i < user->n_eq; ++i) {
//         PetscScalar value;
//         VecGetValues(user->d, 1, &i, &value);
//         PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
//     }
//     PetscPrintf(PETSC_COMM_WORLD, "\n");

    


//     const PetscScalar *d_array;
//     PetscInt d_size;
//     VecGetSize(user->d, &d_size);
//     VecGetArrayRead(user->d, &d_array);
//     for (PetscInt ii = 0; ii < user->n_np * user->n_ed - 1; ++ii){
//         user->disp[ii + 2] = d_array[ii];
//     }
//     user->disp[user->n_np * user->n_ed - 1] = d_array[d_size - 1];
//     VecRestoreArrayRead(user->d, &d_array);
//     std::cout << "disp: " << user->disp.transpose() << std::endl;

//     // 计算增量
//     // user.deltad = user.d - user.d_n;
//     VecCopy(user->d, user->deltad);         // Copy user.d to user.deltad
//     VecAXPY(user->deltad, -1.0, user->d_n);  // Subtract user.d_n from user.deltad,  deltad = -1.0 d_n + deltad

//     user->deltalambda = user->lambda - user->lambda_n;
//     user->deltadisp = user->disp - user->disp_n;

//     VecDestroy(&sol);

//     VecAssemblyBegin(F);
//     VecAssemblyEnd(F);
    
//     return 0;
// }


// 预测步监视器：捕获预测解
// PetscErrorCode SolutionUpdateMonitor(SNES snes, PetscInt it, PetscReal norm, void *ctx) {
//     UserPointersLeeFrame *user = (UserPointersLeeFrame *)ctx;
    
//     //  Store iteration results
//     PetscInt d_size;
//     Vec d_temp;
//     VecGetSize(user->d, &d_size);
//     user->results.d_iter.push_back(user->d);
//     user->results.disp_iter.push_back(user->disp);
//     user->results.lambda_iter.push_back(user->lambda);

//     SNESGetSolutionFuntion

//     return PETSC_SUCCESS;
// }


// 线性求解器监视器：在线性方程组求解完后更新变量
// PetscErrorCode LinearSolveMonitor(KSP ksp, PetscInt n, PetscReal rnorm, void *ctx) {
//     UserPointersLeeFrame *user = (UserPointersLeeFrame *)ctx;
    
   


    
//     return 0;
// }

