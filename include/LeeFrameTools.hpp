#ifndef LEEFRAMETOOLS_HPP
#define LEEFRAMETOOLS_HPP

#include <iostream>
#include <cmath>
#include <vector>
#include "petsc.h"
#include "Eigen"

// 定义结构体 Results
struct Results {
    std::vector<Vec> d_step;        // 等价于 Results.d_step
    std::vector<Eigen::VectorXd> disp_step;     // 等价于 Results.disp_step
    std::vector<double> lambda_step;     // 等价于 Results.lambda_step

    std::vector<Vec> d_iter;        // 等价于 Results.d_iter
    std::vector<Eigen::VectorXd> disp_iter;     // 等价于 Results.disp_iter
    std::vector<double> lambda_iter;     // 等价于 Results.lambda_iter
    std::vector<int> iter_set;
};


// 定义结构体 Nodes
struct Node {
    double x;
    double y;
    double P;
    double Q;
    double M;
    Eigen::VectorXd dof;
};

// 定义结构体 Element
struct Element {
    Node Node1;
    Node Node2;
    double eleLength;
    double eleAngle;
    Eigen::VectorXd eleFint;
    Eigen::MatrixXd eleTangentK;
    Eigen::MatrixXd eleElasticK;
};

inline double arclengthFun(const Eigen::VectorXd& deltaU, double deltaLam, double F_ext_norm, double delta_a) {
    double term1 = deltaU.dot(deltaU); // deltaU' * deltaU
    double term2 = deltaLam * deltaLam * F_ext_norm * F_ext_norm; // (deltaLam * F_ext_norm)^2
    double term3 = delta_a * delta_a; // delta_a squared

    return term1 + term2 - term3;
}


inline double elemLength( const Element &elem, const Eigen::VectorXd &disp){
    Eigen::VectorXd disp_eva = disp.eval();
    // 获取节点1的位移
    double dx1 = disp_eva[elem.Node1.dof(0)];  // x方向位移 (0-based)
    double dy1 = disp_eva[elem.Node1.dof(1)];  // y方向位移

    // 获取节点2的位移
    double dx2 = disp_eva[elem.Node2.dof(0)];
    double dy2 = disp_eva[elem.Node2.dof(1)];

    // 计算新坐标
    double x1_new = elem.Node1.x + dx1;
    double y1_new = elem.Node1.y + dy1;
    double x2_new = elem.Node2.x + dx2;
    double y2_new = elem.Node2.y + dy2;

    // 计算新长度
    double delta_x = x2_new - x1_new;
    double delta_y = y2_new - y1_new;
    return std::sqrt(delta_x * delta_x + delta_y * delta_y);

    // 或者使用Eigen向量计算
    // Eigen::Vector2d vec(delta_x, delta_y);
    // return vec.norm();
}

inline double elemRigidRot(const Element &elem, const Eigen::VectorXd &disp, const Eigen::VectorXd &disp_n){
    Eigen::VectorXd disp_eva = disp.eval();
    Eigen::VectorXd disp_n_eva = disp_n.eval();
    // Element nodes
    Node Node1 = elem.Node1;
    Node Node2 = elem.Node2;

    // Initial load step Nodal displacements
    double dx1_0 = disp_n_eva[Node1.dof(0)];
    double dy1_0 = disp_n_eva[Node1.dof(1)];
    double dx2_0 = disp_n_eva[Node2.dof(0)];
    double dy2_0 = disp_n_eva[Node2.dof(1)];

    // Initial load step Nodal coordinates
    double x1_0 = Node1.x + dx1_0;
    double y1_0 = Node1.y + dy1_0;
    double x2_0 = Node2.x + dx2_0;
    double y2_0 = Node2.y + dy2_0;

    // Initial load step element axial direction
    Eigen::Vector3d vx_0 = Eigen::Vector3d(x2_0 - x1_0, y2_0 - y1_0, 0.0);

    // Current iterate step Nodal displacements
    double dx1_1 = disp_eva[Node1.dof(0)];
    double dy1_1 = disp_eva[Node1.dof(1)];
    double dx2_1 = disp_eva[Node2.dof(0)];
    double dy2_1 = disp_eva[Node2.dof(1)];

    // Current iterate step Nodal coordinates
    double x1_1 = Node1.x + dx1_1;
    double y1_1 = Node1.y + dy1_1;
    double x2_1 = Node2.x + dx2_1;
    double y2_1 = Node2.y + dy2_1;

    // Current iterate step element axial direction
    Eigen::Vector3d vx_1 = Eigen::Vector3d(x2_1 - x1_1, y2_1 - y1_1, 0.0);

    // Increment of the angle
    Eigen::Vector3d dir = vx_0.cross(vx_1);
    double ang_incr = atan2(dir.norm(), vx_0.dot(vx_1));

    // Direction of the rotation
    double s;
    if (dir(2) >= 0.0)  s = 1.0;
    else s = -1.0;
    ang_incr = s * ang_incr;

    return ang_incr;
}

inline Eigen::MatrixXd RotationMatrix(const double &angle){
    // Calculate direction cosines
    double cos_theta = cos(angle);
    double sin_theta = sin(angle);

    // Construct the rotation matrix
    Eigen::MatrixXd T(6, 6);
    T << cos_theta, -sin_theta, 0, 0, 0, 0,
         sin_theta, cos_theta,  0, 0, 0, 0,
         0,         0,          1, 0, 0, 0,
         0,         0,          0, cos_theta, -sin_theta, 0,
         0,         0,          0, sin_theta, cos_theta, 0,
         0,         0,          0, 0, 0, 1;

    return T;
}

inline double elemAngle(const Element &elem, const Eigen::VectorXd &disp){
    // Element nodes
    Node Node1 = elem.Node1;
    Node Node2 = elem.Node2;

    // Nodal displacements
    double dx1 = disp.eval()[Node1.dof(0)];
    double dy1 = disp.eval()[Node1.dof(1)];
    double dx2 = disp.eval()[Node2.dof(0)];
    double dy2 = disp.eval()[Node2.dof(1)];

    // New nodal coordinates
    double x1 = Node1.x + dx1;
    double y1 = Node1.y + dy1;
    double x2 = Node2.x + dx2;
    double y2 = Node2.y + dy2;

    // Compute the new element angle
    double angle = atan2(y2 - y1, x2 - x1);

    return angle;
}





#endif