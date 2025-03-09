#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"
#include "Gauss.hpp"
// #include "Res_Jac2nd.hpp"
#include "PolyBasis.hpp"
#include "HermiteShape.hpp"
#include "Eigen"
#include "UserPointersLeeFrame.hpp"
#include "LeeFrameTools.hpp"
#include "AssemblyK_SecondOrder.hpp"
#include "AssemblyN_SecondOrder.hpp"

// =========================================================================
// This is a numerical case for snap-through buckling of a hinged
// right-angle frame using Euler-Bernoulli beam element
//          |
//          V
//    ----------------------------o
//    |                          / \
//    |                         o---o
//    |
//    |
//    |
//    |
//    |
//    |
//    o
//   / \
//  o---o
//  Consider axial deformation, 
//  element nodal displacement column: d^e = [u1 w1 theta1 u2 w2 theta2]^T
//  element nodal force column: f^e = [fx1 fz1 m1 fx2 fz2 m2]^T
//  element stiffness matrix: k^e : 6 x 6

static char help[] = "Solves a nonlinear system with SNES.\n\n";


int main(int argc, char **args) {

    const double EPS = 1e-14;

    //  model data
    double L = 120;
    double E = 7.2e6;
    double I = 2;
    double A = 6;
    double nu = 0.3;

    double F_applied = -48000;
    double M = 0;
    double Q = 0;
    double g1 = 0;
    double g2 = 0;
    double g3 = 0;
    double g4 = 0;

    // Arc length
    double Deltalambda_init = 0.01;  // initial increment of lamdba
    int max_steps = 10000;           // maximum load steps
    int max_iter  = 1000;            // maximum iterate steps
    double max_ratio = 1.00;
    double tol       = 1e-05;        // tolerance of equilibrium

    // parameters of the FEM
    int n_sd = 2;              // space dimension
    int n_el = 20;             // number of elements
    int n_en = 2;              // number of element nodes
    int n_ed = 3;              // number of element degrees of freedom (per node)
    int deg  = n_en - 1;       // polynomial degree
    int n_np = n_el * deg + 1; // number of points
    int n_eq = n_np * n_ed - 4;       // number of equations
    int n_ee = n_ed * n_en;           // number of equations of an element

    //  quadrature rule
    int n_int = 2;              // number of quadrature points


    // Generate the mesh
    // generate the coordinates of the nodal points
    Eigen::VectorXd x_coor(n_np);
    Eigen::VectorXd y_coor(n_np);
    double hh  = L / (n_el / 2);
    int mid_pt = n_np / 2 + 1;

    for (int n = 0; n < n_np; ++n) {
        if (n <= mid_pt - 1) {
            x_coor(n) = 0;
            y_coor(n) = hh * n / deg;
        } else {
            x_coor(n) = x_coor(mid_pt - 1) + hh * (n - mid_pt + 1) / deg;
            y_coor(n) = hh * (mid_pt - 1) / deg;
        }
    }

    std::cout << "x_coor: " << x_coor.transpose() << std::endl;
    std::cout << "y_coor: " << y_coor.transpose() << std::endl;

    // ID and LM arrays
    Eigen::MatrixXd ID(n_ed, n_np);
    ID.setZero();
    int index = -1;
    for (int AA = 0; AA < n_np; ++AA) {
        for (int ii = 0; ii < n_ed; ++ii) {
            if ((AA == 0 && ii == 0) || (AA == 0 && ii == 1) || 
                (AA == n_np - 1 && ii == 0) || (AA == n_np - 1 && ii == 1)) {
                ID(ii, AA) = -1;
            } else {
                ID(ii, AA) = ++index;
            }
        }
    }
    std::cout << "ID:\n" << ID << std::endl;

    Eigen::MatrixXd IEN(n_en, n_el);
    for (int ee = 0; ee < n_el; ++ee) {
        for (int aa = 0; aa < n_en; ++aa) {
            IEN(aa, ee) = (ee) * (n_en - 1) + aa;
        }
    }
    std::cout << "IEN:\n" << IEN << std::endl;

    Eigen::MatrixXd LM(n_ee, n_el);
    for (int ee = 0; ee < n_el; ++ee) {
        for (int aa = 0; aa < n_en; ++aa) {
            for (int ii = 0; ii < n_ed; ++ii) {
                int pp = n_ed * aa + ii;
                LM(pp, ee) = ID.coeff(ii, IEN.coeff(aa, ee));
            }
        }
    }
    std::cout << "LM:\n" << LM << std::endl;

    std::vector<double> xi(n_int, 0.0);
    std::vector<double> weight(n_int, 0.0);
    Gauss(n_int, -1.0, 1.0, xi, weight, EPS);

    std::cout << "xi: \n";
    for(int ii = 0; ii < n_int; ++ii){
        std::cout << xi[ii] << "\t";
    }
    std::cout << std::endl;
    std::cout << "weight: \n";
    for(int ii = 0; ii < n_int; ++ii){
        std::cout << weight[ii] << "\t";
    }
    std::cout << std::endl;

    // Newton-Raphson Algorithm with using arc-length method
    // initial point value
    Eigen::VectorXd disp(n_np * n_ed);
    disp.setZero();
    std::cout << "disp: " << disp.transpose() << std::endl;
    Eigen::VectorXd disp_n = disp;

    double lambda = 0.0;
    double lambda_n = lambda;

    double deltalambda = 0.0;
    double Deltalambda = 0.0;

    // 初始化 Results 结构体
    Results results;
    results.d_step.reserve(1);      // 预留存储空间
    results.disp_step.reserve(1);
    results.lambda_step.reserve(1);
    results.d_iter.reserve(1);
    results.disp_iter.reserve(1);
    results.lambda_iter.reserve(1);
    results.iter_set.reserve(1);


    results.disp_step.push_back(disp);
    results.lambda_step.push_back(lambda);
    results.disp_iter.push_back(disp);
    results.lambda_iter.push_back(lambda);

    // 输出结果（可选）

    std::cout << "disp_step: " << results.disp_step[0].transpose() << std::endl;
    std::cout << "lambda_step: " << results.lambda_step[0] << std::endl;
    std::cout << "disp_iter: " << results.disp_iter[0].transpose() << std::endl;
    std::cout << "lambda_iter: " << results.lambda_iter[0] << std::endl;

    // F_ext
    Eigen::VectorXd Force_ext(n_np * n_ed);
    int Force_num = ((n_el / 2 + 1) + n_el/10)*3 - 1;
    Force_ext(Force_num - 1) = F_applied;
    std::cout << "Force_ext: " << Force_ext.transpose() << std::endl;



    // Struct Node
    std::vector<Node> Nodes(n_np);
    // 填充节点信息
    for (int ii = 0; ii < n_np; ++ii) {
        Nodes[ii].x = x_coor(ii);
        Nodes[ii].y = y_coor(ii);
        Nodes[ii].P = Force_ext(n_ed * ii);     // 索引转换为 0-based
        Nodes[ii].Q = Force_ext(n_ed * ii + 1);
        Nodes[ii].M = Force_ext(n_ed * ii + 2);

        Eigen::VectorXd dof(n_ed);
        for (int jj = 0; jj < n_ed; ++jj) {
            dof(jj) = n_ed * ii + jj;
        }
        Nodes[ii].dof = dof;
    }
    for (int ii = 0; ii < n_np; ++ii){
        std::cout << "Nodes[" << ii << "].x: " << Nodes[ii].x << "\t";
        std::cout << "Nodes[" << ii << "].y: " << Nodes[ii].y << "\t";
        std::cout << "Nodes[" << ii << "].P: " << Nodes[ii].P << "\t";
        std::cout << "Nodes[" << ii << "].Q: " << Nodes[ii].Q << "\t";
        std::cout << "Nodes[" << ii << "].M: " << Nodes[ii].M << "\t";
        std::cout << "Nodes[" << ii << "].dof: " << Nodes[ii].dof.transpose() << "\t";
        std::cout << std::endl;
    }

    // Updated Lagrangian 
    // Elem: save all the information of the element at the n-th loading step
    // 创建 Element 数组
    std::vector<Element> Elems;
    Elems.resize(n_el);

    // 初始化 Element
    for (int ee = 0; ee < n_el; ++ee) {
        int node1_idx = IEN(0, ee);
        int node2_idx = IEN(1, ee);

        // 获取节点信息
        Node Node1 = Nodes[node1_idx];
        Node Node2 = Nodes[node2_idx];

        // 计算节点向量
        Eigen::Vector3d vl(Node2.x - Node1.x, Node2.y - Node1.y, 0.0);
        double length = vl.norm();

        // 计算方向向量
        Eigen::Vector3d vx = vl / length;
        Eigen::Vector3d vz(0.0, 0.0, 1.0);
        Eigen::Vector3d vy = vz.cross(vx); // 叉积

        // 初始化 Element
        Elems[ee].Node1 = Node1;
        Elems[ee].Node2 = Node2;
        Elems[ee].eleLength = length;
        Elems[ee].eleAngle = atan2(Node2.y - Node1.y, Node2.x - Node1.x); // 弧度制

        // 初始化矩阵和力向量
        Elems[ee].eleFint = Eigen::VectorXd::Zero(n_ee);
        Elems[ee].eleTangentK = Eigen::MatrixXd::Zero(n_ee, n_ee);
        Elems[ee].eleElasticK = Eigen::MatrixXd::Zero(n_ee, n_ee);
    }

    for (int ee = 0; ee < n_el; ++ee){
        std::cout << "Elems[" << ee << "].Node1.x: " << Elems[ee].Node1.x << "\n";
        std::cout << "Elems[" << ee << "].Node2.x: " << Elems[ee].Node2.x << "\n";
        std::cout << "Elems[" << ee << "].eleLength: " << Elems[ee].eleLength << "\n";
        std::cout << "Elems[" << ee << "].eleAngle: " << Elems[ee].eleAngle << "\n";
        std::cout << "Elems[" << ee << "].Fint: " << Elems[ee].eleFint << "\n";
        std::cout << "Elems[" << ee << "].eleTangentK: " << Elems[ee].eleTangentK << "\n";
        std::cout << "Elems[" << ee << "].eleElasticK: " << Elems[ee].eleElasticK << "\n";
    }
    std::vector<Element> Elems_n = Elems;

    // Setup the stiffness matrix and load vector 
    UserPointersLeeFrame user; /* user-defined work context */

    //  model data
    user.L = L;
    user.E = E;
    user.I = I;
    user.A = A;
    user.nu = nu;

    user.F_applied = F_applied;
    user.M = M;
    user.Q = Q;
    user.g1 = g1;
    user.g2 = g2;
    user.g3 = g3;
    user.g4 = g4;

    user.Deltalambda_init = Deltalambda_init;  // initial increment of lamdba
    user.max_steps = max_steps;           // maximum load steps
    user.max_iter  = max_iter;            // maximum iterate steps
    user.max_ratio = max_ratio;
    user.tol       = tol;                 // tolerance of equilibrium

    // parameters of the FEM
    user.n_sd = n_sd;             // space dimension
    user.n_el = n_el;             // number of elements
    user.n_en = n_en;             // number of element nodes
    user.n_ed = n_ed;             // number of element degrees of freedom (per node)
    user.deg  = deg;              // polynomial degree
    user.n_np = n_np;             // number of points
    user.n_eq = n_eq;             // number of equations
    user.n_ee = n_ee;             // number of equations of an element

    //  quadrature rule
    user.n_int = n_int;           // number of quadrature points
    user.xi = xi;
    user.weight = weight;

    // mesh
    user.x_coor = x_coor;
    user.y_coor = y_coor;
    user.ID = ID;
    user.IEN = IEN;
    user.LM = LM;

    // solution
    user.disp = disp;
    user.disp_n = disp_n;
    user.lambda = lambda;
    user.lambda_n = lambda_n;
    user.deltalambda = deltalambda;
    user.Deltalambda = Deltalambda;
    user.Force_ext = Force_ext;
    user.Nodes = Nodes;
    user.Elems = Elems;
    user.Elems_n = Elems_n;

    SNES           snes;         /* nonlinear solver context */
    KSP            ksp;          /* linear solver context */
    // PC             pc;           /* preconditioner context */
    SNESLineSearch linesearch;   /* line search context */
    Vec            x,r;          /* solution, residual vectors */
    Mat            J, K0;            /* Jacobian matrix */
    PetscErrorCode ierr;
    // PetscMPIInt    size;
    PetscReal      rtol = 1.0e-8;  // 相对残差阈值
    PetscReal      atol = 1.0e-50;  // 绝对残差阈值
    PetscReal      stol = 1.0e-8;  // 范数残差
    PetscInt       maxit = 20;  // 最大迭代次数
    PetscInt       maxf = 16; // 最大函数评估次数

    PetscViewer    viewer;
    

    PetscInitialize(&argc,&args,(char*)0,help);
    // MPI_Comm_size(PETSC_COMM_WORLD,&size);
    // if (size > 1) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Example is only for sequential runs");

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer, PETSCVIEWERASCII); // ASCII 格式
    PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_COMMON);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Create nonlinear solver context
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    SNESCreate(PETSC_COMM_WORLD,&snes);
    SNESSetTolerances(snes, atol, rtol, stol, maxit, maxf);
    SNESSetType(snes, SNESNEWTONAL);
    // SNESSetType(snes, SNESNEWTONLS);
    SNESGetLineSearch(snes, &linesearch);
    SNESLineSearchSetType(linesearch, SNESLINESEARCHNONE);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Create matrix and vector data structures; set corresponding routines
        - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
        Create vectors for solution and nonlinear function
    */
    VecCreate(PETSC_COMM_WORLD,&x);
    VecSetSizes(x,PETSC_DECIDE,n_eq);
    VecSetFromOptions(x);
    VecDuplicate(x,&r);

    Vec d, d_n, Deltad, Deltad_n, deltad, q_bar, q_bar_n, R, N, Deltad_bar; 
    VecCreate(PETSC_COMM_WORLD,&d);
    VecSetType(d, VECSEQ);  
    VecSetSizes(d,PETSC_DECIDE,n_eq);
    VecSet(d,0.0);
    VecDuplicate(d,&d_n);
    VecDuplicate(d,&Deltad);      // increment of the iteration step
    VecDuplicate(d,&Deltad_n);
    VecDuplicate(d,&deltad);     // increment of the load step
    VecDuplicate(d,&q_bar); 
    VecDuplicate(d,&q_bar_n);
    VecDuplicate(d,&R);
    VecDuplicate(d,&N);
    VecDuplicate(d,&Deltad_bar);

    user.d = d;
    user.d_n = d_n;
    user.Deltad = Deltad;
    user.Deltad_n = Deltad_n;
    user.deltad = deltad;
    user.q_bar = q_bar;
    user.q_bar_n = q_bar_n;
    user.R = R;
    user.N = N;
    user.Deltad_bar = Deltad_bar;
    results.d_step.push_back(d);
    results.d_iter.push_back(d);
    user.results = results;
    std::cout << "user.d: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.d, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.d_n: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.d_n, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.Deltad: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.Deltad, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.Deltad_n: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.Deltad_n, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.deltad: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.deltad, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.q_bar: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.q_bar, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.q_bar_n: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.q_bar_n, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.R: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.R, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.N: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.N, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "user.Deltad_bar: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(user.Deltad_bar, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "d_step[0]: \n";
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(results.d_step[0], 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "d_iter[0]: \n";
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(results.d_iter[0], 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    // VecView(results.d_iter[0], PETSC_VIEWER_STDOUT_WORLD );

    Vec F_ext; 
    VecCreate(PETSC_COMM_WORLD,&F_ext);
    VecSetType(F_ext, VECSEQ); 
    VecSetSizes(F_ext,PETSC_DECIDE,n_eq);
    PetscInt F_num = ((n_el / 2 + 1) + n_el/10)*3 - 1 - 2 -1;
    VecSetValues(F_ext, 1, &F_num, &F_applied, INSERT_VALUES);
    PetscReal F_ext_norm;
    VecNorm(F_ext, NORM_2, &F_ext_norm);

    user.F_ext = F_ext;
    std::cout << "F_ext: \n"; 
    for (PetscInt i = 0; i < user.n_eq; ++i) {
        PetscScalar value;
        VecGetValues(F_ext, 1, &i, &value);
        PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
    }
    PetscPrintf(PETSC_COMM_WORLD, "\n");
    std::cout << "F_ext_norm: " << F_ext_norm << std::endl;

    /*
        Create Jacobian matrix data structure
    */
    MatCreate(PETSC_COMM_WORLD,&J);
    MatSetType(J, MATAIJ);
    MatSetSizes(J,PETSC_DECIDE,PETSC_DECIDE,n_eq,n_eq);
    MatSetFromOptions(J);
    MatSetUp(J);
    MatZeroEntries(J);
    for (PetscInt i = 0; i < n_eq; i++) {
        PetscInt col = i;
        PetscScalar value = 1.0; // 至少设置对角线元素
        MatSetValue(J, i, col, value, INSERT_VALUES);
    }
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

    MatCreate(PETSC_COMM_WORLD, &K0);
    MatSetType(K0, MATAIJ);
    MatSetSizes(K0, PETSC_DECIDE, PETSC_DECIDE, n_eq, n_eq);
    MatSetFromOptions(K0);
    MatSetUp(K0);
    MatZeroEntries(K0);
    for (PetscInt i = 0; i < n_eq; i++) {
        PetscInt col = i;
        PetscScalar value = 1.0; // 至少设置对角线元素
        MatSetValue(K0, i, col, value, INSERT_VALUES);
    }
    MatAssemblyBegin(K0, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K0, MAT_FINAL_ASSEMBLY);

    MatView(J, viewer);
    MatView(K0, viewer); 



    /*
    Set function evaluation routine and vector.
    */

    // SNESSetFunction(snes,r,FormFunction2nd,&user);

    /*
    Set Jacobian matrix data structure and Jacobian evaluation routine
    */
    // SNESSetJacobian(snes,J,J,FormJacobian2nd,&user);

    // 设置弧长约束函数
    // SNESNewtonALSetStepSize(snes, 480.0);
    // SNESNewtonALSetConsistencyParameter(snes, 0.0); 
    // SNESNewtonALSetFunction(snes, ComputeTangentLoad, &user);

    // 设置修正类型为法线修正
    // SNESNewtonALSetCorrectionType(snes, SNES_NEWTONAL_CORRECTION_EXACT);

    // 设置监视器，获取预测步的值
    // SNESMonitorSet(snes, SolutionUpdateMonitor, &user, NULL);

    // 设置 KSP 监视器
    // KSPMonitorSet(ksp, LinearSolveMonitor, &user, NULL);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Customize nonlinear solver; set runtime options
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
    Set linear solver defaults for this problem. By extracting the
    KSP and PC contexts from the SNES context, we can then
    directly call any KSP and PC routines to set various options.
    */
    PetscReal l_rtol = 1.0e-5;
    PetscReal l_atol = 1.0e-50;
    PetscReal l_dtol = 1.0e20;
    PetscInt l_maxits = 200;
    SNESGetKSP(snes,&ksp);
    KSPSetTolerances(ksp,l_rtol, l_atol, l_dtol, l_maxits);
    KSPSetFromOptions(ksp);
    // KSPGetPC(ksp, &pc);
    // PCSetType(pc, PCLU);
    // PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO); // 避免零主元
    KSPSetFromOptions(ksp);
    
    

    /*
    Set SNES/KSP/KSP/PC runtime options, e.g.,
        -snes_view -snes_monitor -ksp_type <ksp> -pc_type <pc>
    These options will override those specified above as long as
    SNESSetFromOptions() is called _after_ any other customization
    routines.
    */
   
    SNESSetFromOptions(snes);


    // start the process
    int n = 0;
    Eigen::VectorXd deltadisp(n_np * n_ed);
    user.deltadisp = deltadisp;
    while (n < max_steps){
        // Initialization of all loading step
        n += 1;

        VecCopy(user.d, user.d_n);
        user.lambda_n = user.lambda;
        disp_n = disp;

        for (int ee = 0; ee < n_el; ++ee){
            user.Elems_n[ee].eleLength = elemLength(user.Elems[ee], disp);
            user.Elems_n[ee].eleAngle  = user.Elems[ee].eleAngle;
            user.Elems_n[ee].eleFint  = user.Elems[ee].eleFint;
            user.Elems_n[ee].eleTangentK  = user.Elems[ee].eleTangentK;
            user.Elems_n[ee].eleElasticK  = user.Elems[ee].eleElasticK;
        }
        AssemblyK_SecondOrder(K0, &user);
        MatView(K0, viewer); 

        PetscScalar Mij;
        MatGetValue(K0, 0, 0, &Mij);
        PetscPrintf(PETSC_COMM_WORLD, "Mij: ");
        PetscPrintf(PETSC_COMM_WORLD, "%f ", Mij);
        PetscPrintf(PETSC_COMM_WORLD, "\n");

        // MatSolve(K0, user.F_ext, user.q_bar);
        KSPSetOperators(ksp, K0, K0);
        KSPSolve(ksp, user.F_ext, user.q_bar);

        VecView(user.q_bar, viewer);
        
        int s = 0;
        PetscScalar n2 = 0.0;
        if (n == 1){
            s = 1;
            user.Deltalambda = Deltalambda_init;
            VecDot(q_bar, q_bar, &n2);
        }
        else{
            PetscScalar m2;
            VecDot(q_bar_n, q_bar, &m2);
            PetscScalar GSP = n2 / m2;
            if (GSP < 0){
                s = -s;
            }
            // 柱面弧长控制公式
            PetscScalar n2_temp;
            PetscScalar d2_temp;
            VecDot(user.deltad, user.deltad, &d2_temp);
            VecDot(user.q_bar, user.q_bar, &n2_temp);
            user.Deltalambda = s * sqrt(d2_temp/ n2_temp);

            // 荷载比例限制
            if ((max_ratio > 0 && (user.lambda + user.Deltalambda) > max_ratio) || (max_ratio < 0 && (user.lambda + user.Deltalambda) < max_ratio)) {
                user.Deltalambda = max_ratio - user.lambda;
            }
            
        }


        //user.Deltad = user.Deltalambda * q_bar;
        VecCopy(user.q_bar, user.Deltad);         // Copy q_bar to Deltad   
        VecScale(user.Deltad, user.Deltalambda);  // Scale Deltad by Deltalambda

        std::cout << "user.d: \n"; 
        for (PetscInt i = 0; i < user.n_eq; ++i) {
            PetscScalar value;
            VecGetValues(user.d, 1, &i, &value);
            PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");

        std::cout << "user.d_n: \n"; 
        for (PetscInt i = 0; i < user.n_eq; ++i) {
            PetscScalar value;
            VecGetValues(user.d_n, 1, &i, &value);
            PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");
        

        VecAXPY(user.d, 1.0, user.Deltad);    // 位移更新, d = d + Deltad
        user.lambda += user.Deltalambda;       // 荷载因子更新

        std::cout << "user.d: \n"; 
        for (PetscInt i = 0; i < user.n_eq; ++i) {
            PetscScalar value;
            VecGetValues(user.d, 1, &i, &value);
            PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");

        std::cout << "user.d_n: \n"; 
        for (PetscInt i = 0; i < user.n_eq; ++i) {
            PetscScalar value;
            VecGetValues(user.d_n, 1, &i, &value);
            PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");


        // 处理边界条件（例如disp的拼接）
        const PetscScalar *d_array;
        PetscInt d_size;
        VecGetSize(user.d, &d_size);
        VecGetArrayRead(user.d, &d_array); 
        Eigen::VectorXd d_eigen = Eigen::Map<const Eigen::VectorXd>(d_array, d_size);
        VecRestoreArrayRead(user.d, &d_array);
        // 分块赋值
        user.disp << 0.0,                     // g_BC[0] 对应 MATLAB 的 g_BC(1)
                     0.0,                     // g_BC[1] 对应 MATLAB 的 g_BC(2)
                     d_eigen.head(d_size - 1),  // d(1:end-1)
                     0.0,                     // g_BC[2] 对应 MATLAB 的 g_BC(3)
                     0.0,                     // g_BC[3] 对应 MATLAB 的 g_BC(4)
                     d_eigen.tail(1);          // d(end)
        std::cout << "user.disp: " << user.disp.transpose() << std::endl;
        

        // 计算增量
        // user.deltad = user.d - user.d_n;
        VecCopy(user.d, user.deltad);         // Copy user.d to user.deltad
        VecAXPY(user.deltad, -1.0, user.d_n);  // Subtract user.d_n from user.deltad,  deltad = -1.0 d_n + deltad

        user.deltalambda = user.lambda - user.lambda_n;
        user.deltadisp = user.disp - user.disp_n;

        std::cout << "user.deltadisp: " << user.deltadisp.transpose() << std::endl;

        // Stored the current load step initial increment for the next load step
        VecCopy(user.Deltad, user.Deltad_n);
        double Deltalambda_n = Deltalambda;
        VecCopy(q_bar, q_bar_n);       // will be used for the next GSP calculation

        std::cout << "user.d: \n"; 
        for (PetscInt i = 0; i < user.n_eq; ++i) {
            PetscScalar value;
            VecGetValues(user.d, 1, &i, &value);
            PetscPrintf(PETSC_COMM_WORLD, "%f ", PetscRealPart(value));
        }
        PetscPrintf(PETSC_COMM_WORLD, "\n");

        // Solve
        // SNESSolve(snes,NULL,user.d);

        int iter_step = 0;
        PetscInt conv = 0;
        double tol = 1e-6;

        while (!conv && iter_step < max_iter) {
            AssemblyN_SecondOrder(user.N, &user);
            VecWAXPY(user.R, -user.lambda, user.F_ext, user.N); // R = -lambda * F_ext + N
            VecScale(user.R, -1.0);              // R = -R =  lambda * F - N
            PetscReal R_norm;
            VecNorm(user.R, NORM_2, &R_norm);

            // Store iteration results
            user.results.d_iter.push_back(user.d);
            user.results.disp_iter.push_back(user.disp);
            user.results.lambda_iter.push_back(user.lambda);

            // Check for convergence
            if (R_norm / F_ext_norm < tol) conv = 1;
            if (conv == 1) break;

            iter_step += 1;
            AssemblyK_SecondOrder(K0, &user);

            KSPSetOperators(ksp, K0, K0);
            KSPSolve(ksp, user.F_ext, user.q_bar);

            KSPSetOperators(ksp, K0, K0);
            KSPSolve(ksp, user.R, user.Deltad_bar);
            

            // Constant arc-length cylindrical
            PetscScalar a, b, c, tempScalar;
            Vec tempVec;
            VecDuplicate(d,&tempVec);
            // a
            VecDot(user.q_bar, user.q_bar, &a);
            VecWAXPY(tempVec, 1.0, user.Deltad_bar, user.deltad); // tempVec = Deltad_bar + deltad
            // b
            VecDot(user.q_bar, tempVec, &b);
            VecWAXPY(tempVec, 2.0, user.deltad, user.Deltad_bar); // tempVec = Deltad_bar + 2 * deltad
            // c
            VecDot(user.Deltad_bar, tempVec, &c);
            // s_iter
            VecDot(user.deltad, user.q_bar, &tempScalar);
            PetscScalar s_iter;
            if (tempScalar >= 0) s_iter = 1.0;
            else s_iter = -1.0;

            // 计算 Deltalambda
            PetscScalar sqrt_term;
            sqrt_term = (b / a) * (b / a) - c / a;
            if (sqrt_term < 0) {
                std::cerr << "Error: sqrt_term is negative! Setting Deltalambda to 0.\n";
                conv = -1;
                break;
            } else {
                user.Deltalambda = -b / a + s_iter * sqrt(sqrt_term);
            }

            // Deltad^(i) computation
            VecWAXPY(user.Deltad, user.Deltalambda, user.q_bar, user.Deltad_bar);

            // Update total values
            VecAXPY(user.d, 1.0, user.Deltad);    // 位移更新, d = d + Deltad
            user.lambda += user.Deltalambda;       // 荷载因子更新

            // 处理边界条件（例如disp的拼接）
            const PetscScalar *d_array;
            PetscInt d_size;
            VecGetSize(user.d, &d_size);
            VecGetArrayRead(user.d, &d_array); 
            Eigen::VectorXd d_eigen = Eigen::Map<const Eigen::VectorXd>(d_array, d_size);
            VecRestoreArrayRead(user.d, &d_array);
            // 分块赋值
            user.disp << 0.0,                     // g_BC[0] 对应 MATLAB 的 g_BC(1)
                        0.0,                     // g_BC[1] 对应 MATLAB 的 g_BC(2)
                        d_eigen.head(d_size - 1),  // d(1:end-1)
                        0.0,                     // g_BC[2] 对应 MATLAB 的 g_BC(3)
                        0.0,                     // g_BC[3] 对应 MATLAB 的 g_BC(4)
                        d_eigen.tail(1);          // d(end)
            

            // 计算增量
            // user.deltad = user.d - user.d_n;
            VecCopy(user.d, user.deltad);         // Copy user.d to user.deltad
            VecAXPY(user.deltad, -1.0, user.d_n);  // Subtract user.d_n from user.deltad,  deltad = -1.0 d_n + deltad

            user.deltalambda = user.lambda - user.lambda_n;
            user.deltadisp = user.disp - user.disp_n;

            VecDestroy(&tempVec);
        }

        results.iter_set.push_back(iter_step);

        if (conv == 0) {
            std::cout << "Convergence not reached!\nStep = " << n << "\nIter = " << iter_step << std::endl;
            return 1;
        }
        else if (conv == -1){
            std::cout << "Unable to compute load increment!\nStep = " << n << "\nIter = " << iter_step << std::endl;
            return 1;
        }

        // Store the solutions: equilibrium configuration
        user.results.d_step.push_back(user.d);
        user.results.disp_step.push_back(user.disp);
        user.results.lambda_step.push_back(user.lambda);
        

        std::cout << "Step: " << n << " | Ratio: " << user.lambda << " | Iter: " << iter_step << std::endl;
        if (user.lambda >= 0.999 * max_ratio) break;

    }

    std::cout << "user.deltadisp: " << user.deltadisp.transpose() << std::endl;

    // PetscViewerDestroy(&viewer);
    VecDestroy(&x); VecDestroy(&r); VecDestroy(&F_ext);
    VecDestroy(&d); VecDestroy(&d_n); VecDestroy(&Deltad);
    VecDestroy(&deltad); VecDestroy(&q_bar); VecDestroy(&q_bar_n);
    VecDestroy(&R); VecDestroy(&N); VecDestroy(&Deltad_bar);
    MatDestroy(&J); MatDestroy(&K0); 
    SNESDestroy(&snes);      // 只需要销毁snes
    ierr = PetscFinalize();

    std::cout << "user.deltadisp: " << user.deltadisp.transpose() << std::endl;
    return ierr;

}

