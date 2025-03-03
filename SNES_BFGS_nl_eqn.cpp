#include <petscsnes.h>

typedef struct {
    PetscReal x;          // 材料参数x（15或25）
    PetscReal load_step;  // 载荷步长
    PetscInt  current_step; // 当前载荷步
} UserContext;

static PetscErrorCode FormFunction(SNES snes, Vec d, Vec F, void *ctx) {
    UserContext *user_ctx = (UserContext *)ctx;
    const PetscScalar *d_array;
    PetscScalar *f_array;

    VecGetArrayRead(d, &d_array);
    VecGetArray(F, &f_array);

    PetscReal d1 = d_array[0], d2 = d_array[1];
    PetscReal F1 = user_ctx->load_step * user_ctx->current_step;

    // 计算内部力
    PetscReal N1 = user_ctx->x * d1 / (10.0 - d1) - 0.5 * d2 * d2;
    PetscReal N2 = d2 - d1;

    // 残差计算
    f_array[0] = N1 - F1;
    f_array[1] = N2 - 0.0;

    VecRestoreArrayRead(d, &d_array);
    VecRestoreArray(F, &f_array);
    return PETSC_SUCCESS;
}

int main(int argc, char **argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    // 问题参数设置
    const PetscReal load_step = 0.25;
    const PetscReal load_max = 10.0;
    const PetscInt  load_num = (PetscInt)(load_max / load_step);
    const PetscReal x_values[2] = {15.0, 25.0}; // 分别计算x=15和x=25的情况

    // PETSc对象
    Vec         d, F;
    SNES        snes;
    UserContext user_ctx;
    PetscInt    indices[2] = {0, 1};

    // 初始化上下文
    user_ctx.load_step = load_step;

    // 创建向量
    VecCreate(PETSC_COMM_WORLD, &d);
    VecSetSizes(d, PETSC_DECIDE, 2);
    VecSetFromOptions(d);
    VecDuplicate(d, &F);

    // 主循环：分别处理x=15和x=25的情况
    for (int case_id = 0; case_id < 2; ++case_id) {
        user_ctx.x = x_values[case_id];
        
        // 初始化解向量
        PetscScalar d_initial[2] = {0.0, 0.0};
        VecSetValues(d, 2, indices, d_initial, INSERT_VALUES);
        VecAssemblyBegin(d);
        VecAssemblyEnd(d);

        // 创建SNES求解器
        SNESCreate(PETSC_COMM_WORLD, &snes);
        SNESSetType(snes, SNESQN);
        SNESQNSetType(snes, SNES_QN_BFGS); // 设置BFGS方法
        SNESSetFunction(snes, F, FormFunction, &user_ctx);
        
        // 设置收敛容差
        SNESSetTolerances(snes, 1e-4, PETSC_DEFAULT, PETSC_DEFAULT, 
                         PETSC_DEFAULT, PETSC_DEFAULT);

        // 分步加载
        PetscPrintf(PETSC_COMM_WORLD, "\n=== Solving for x = %g ===\n", user_ctx.x);
        for (user_ctx.current_step = 1; user_ctx.current_step <= load_num; ++user_ctx.current_step) {
            // 设置初始猜测为前一步的解
            SNESSolve(snes, NULL, d);

            // 获取并打印解
            const PetscScalar *sol;
            VecGetArrayRead(d, &sol);
            PetscPrintf(PETSC_COMM_WORLD, "Load step %3d: d1 = %8.5f, d2 = %8.5f\n",
                       user_ctx.current_step, (double)sol[0], (double)sol[1]);
            VecRestoreArrayRead(d, &sol);
        }
        SNESDestroy(&snes);
    }

    // 清理资源
    VecDestroy(&d);
    VecDestroy(&F);
    PetscFinalize();
    return 0;
}