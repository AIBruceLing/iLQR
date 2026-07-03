#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "matplotlibcpp.h"
 
namespace plt = matplotlibcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;
 
// 车辆模型类
class Vehicle {
public:
    double x, y, theta, v;
 
    // 车辆构造函数，初始化车辆的状态
    Vehicle(double x_init = 0.0, double y_init = 0.0, double theta_init = 0.0, double v_init = 0.0)
        : x(x_init), y(y_init), theta(theta_init), v(v_init) {}
 
    // 更新车辆状态，a为加速度，omega为角速度，dt为时间步长
    void update(double a, double omega, double dt) {
        x += v * cos(theta) * dt;    // 更新x坐标
        y += v * sin(theta) * dt;    // 更新y坐标
        theta += omega * dt;         // 更新航向角theta
        v += a * dt;                 // 更新速度v
    }
};
 
// 轨迹类，用于存储参考轨迹信息
class Trajectory {
public:
    std::vector<double> cx, cy, theta, v;
 
    // 轨迹构造函数，初始化参考轨迹
    Trajectory() {
        // 生成参考轨迹的x和y坐标
        for (double i = 0; i < 50; i += 0.1) {
            cx.push_back(i);
            cy.push_back(sin(i / 5.0) * i / 2.0);
        }
        // 计算每个轨迹点的航向角theta
        for (size_t i = 0; i < cx.size() - 1; ++i) {
            theta.push_back(atan2(cy[i+1] - cy[i], cx[i+1] - cx[i]));
        }
        theta.push_back(0.0);  // 最后一个点的theta设置为0
        v.assign(cx.size(), 3.0);  // 设置所有点的目标速度为3 m/s
        //assign：先清空向量 v 原有的所有内容，然后将 v 的大小重新调整为 cx.size()，最后把这 cx.size() 个格子全部填充为双精度浮点数 3.0
    }
 
    // 获取指定索引处的参考状态向量[x, y, theta, v]
    VectorXd get_reference(size_t index) const {
        VectorXd ref(4);
        ref << cx[index], cy[index], theta[index], v[index];
        return ref;
    }
};
 
// iLQR（迭代线性二次调节器）控制器类
class iLQRController {
public:
    int N;             // 控制步数
    int max_iter;      // 最大迭代次数
    double dt;         // 时间步长
    MatrixXd Q, R, Qf; // 代价函数的权重矩阵，Q和R是stage cost，Qf是terminal cost
 
    // iLQR控制器构造函数
    iLQRController(int N_input = 50, int max_iter_input = 10, double dt_input = 0.1)
        : N(N_input), max_iter(max_iter_input), dt(dt_input) {
        // 初始化状态代价矩阵Q，控制代价矩阵R，终端状态代价矩阵Qf
        Q = MatrixXd::Zero(4, 4);//这里状态量有四个，分别是x,y,theta,v，所以这里的Q矩阵大小是4*4
        Q(0, 0) = 1.0; Q(1, 1) = 1.0; Q(2, 2) = 0.5; Q(3, 3) = 0.1;//
        R = MatrixXd::Identity(2, 2) * 0.1;//这里控制量有两个，分别是a,omega，所以这里的R矩阵大小是2*2
        Qf = Q * 10.0;//终端状态代价矩阵Qf是Q的10倍
    }
 
    // iLQR算法的实现，返回当前时刻的最优控制输入
    VectorXd ilqr(Vehicle &vehicle, const Trajectory &trajectory, size_t index) {
        int x_dim = 4;  // 状态维度 [x, y, theta, v]
        int u_dim = 2;  // 控制维度 [a, omega]
        std::vector<VectorXd> xs(N + 1, VectorXd::Zero(x_dim));  // 状态序列
        std::vector<VectorXd> us(N, VectorXd::Zero(u_dim));      // 控制序列
 
        // 初始化状态
        xs[0] << vehicle.x, vehicle.y, vehicle.theta, vehicle.v;
 
        // 迭代更新控制输入和状态
        for (int iter = 0; iter < max_iter; ++iter) {
            // 前向传播，基于当前控制序列预测状态序列
            for (int i = 0; i < N; ++i) {
                xs[i + 1] = dynamics(xs[i], us[i], dt);
            }
 
            // 反向传播，更新控制增量和反馈矩阵
            VectorXd Vx = Qf * (xs[N] - trajectory.get_reference(index + N));// 终端状态的价值函数梯度，是Qf*(xN-x_refN)
            MatrixXd Vxx = Qf;                                               // 终端状态的价值函数Hessian矩阵
 
            std::vector<VectorXd> k_control_list(N, VectorXd::Zero(u_dim));          // 控制增量
            std::vector<MatrixXd> K_feedback_list(N, MatrixXd::Zero(u_dim, x_dim));  // 反馈矩阵
 
            for (int i = N - 1; i >= 0; --i) {
                VectorXd x_ref = trajectory.get_reference(index + i);
                std::pair<MatrixXd, MatrixXd> linearized = linearize_dynamics(xs[i], us[i], dt);
                MatrixXd fx = linearized.first;   // 线性化动力学模型的状态矩阵，雅可比矩阵
                MatrixXd fu = linearized.second;  // 线性化动力学模型的控制矩阵，雅可比矩阵
 
                // 计算代价梯度和Hessian矩阵
                VectorXd lx = Q * (xs[i] - x_ref);
                VectorXd lu = R * us[i];
                MatrixXd lxx = Q;
                MatrixXd luu = R;
                MatrixXd lux = MatrixXd::Zero(u_dim, x_dim);
 
                // 计算Q函数的二次近似
                VectorXd Qx = lx + fx.transpose() * Vx;
                VectorXd Qu = lu + fu.transpose() * Vx;
                MatrixXd Qxx = lxx + fx.transpose() * Vxx * fx;
                MatrixXd Quu = luu + fu.transpose() * Vxx * fu;
                MatrixXd Qux = lux + fu.transpose() * Vxx * fx;
 
                // 计算控制增量和反馈矩阵
                MatrixXd Quu_inv = Quu.inverse();
                VectorXd k_control = -Quu_inv * Qu;
                MatrixXd K_feedback = -Quu_inv * Qux;
 
                // 更新价值函数
                Vx = Qx + K_feedback.transpose() * Quu * k_control + K_feedback.transpose() * Qu + Qux.transpose() * k_control;
                Vxx = Qxx + K_feedback.transpose() * Quu * K_feedback + K_feedback.transpose() * Qux + Qux.transpose() * K_feedback;
 
                k_control_list[i] = k_control;
                K_feedback_list[i] = K_feedback;
            }
 
            // 更新控制序列并进行前向模拟
            std::vector<VectorXd> xs_new(N + 1, xs[0]);
            std::vector<VectorXd> us_new;
 
            for (int i = 0; i < N; ++i) {
                VectorXd du = k_control_list[i] + K_feedback_list[i] * (xs_new[i] - xs[i]);
                us_new.push_back(us[i] + du);
                xs_new[i + 1] = dynamics(xs_new[i], us_new.back(), dt);
            }
            xs = xs_new;
            us = us_new;
 
            // 判断是否收敛
            double cost = compute_total_cost(xs, us, trajectory, index);
            std::cout << "Iteration " << iter << ", Cost: " << cost << std::endl;
            if (cost < 1e-6) {
                break;
            }
        }
        return us[0];  // 返回当前时刻的最优控制输入
    }
 
    // 车辆动力学模型
    VectorXd dynamics(const VectorXd &x, const VectorXd &u, double dt) {
        VectorXd x_next = x;
        x_next[0] += x[3] * cos(x[2]) * dt;  // 更新x坐标
        x_next[1] += x[3] * sin(x[2]) * dt;  // 更新y坐标
        x_next[2] += u[1] * dt;              // 更新theta
        x_next[3] += u[0] * dt;              // 更新速度v
        return x_next;
    }
 
    // 线性化车辆动力学模型
    std::pair<MatrixXd, MatrixXd> linearize_dynamics(const VectorXd &x, const VectorXd &u, double dt) {
        MatrixXd fx = MatrixXd::Identity(4, 4);
        fx(0, 2) = -x[3] * sin(x[2]) * dt;
        fx(0, 3) = cos(x[2]) * dt;
        fx(1, 2) = x[3] * cos(x[2]) * dt;
        fx(1, 3) = sin(x[2]) * dt;
 
        MatrixXd fu = MatrixXd::Zero(4, 2);
        fu(2, 1) = dt;
        fu(3, 0) = dt;
 
        return {fx, fu};// 返回线性化+离散化（前向欧拉离散）的状态矩阵和控制矩阵，雅可比矩阵
    }
 
    // 计算总成本
    double compute_total_cost(const std::vector<VectorXd> &xs, const std::vector<VectorXd> &us, const Trajectory &trajectory, size_t index) {
        double cost = 0.0;
        for (int i = 0; i < N; ++i) {
            VectorXd dx = xs[i] - trajectory.get_reference(index + i);
            cost += (dx.transpose() * Q * dx)(0, 0) + (us[i].transpose() * R * us[i])(0, 0);
        }
        VectorXd dx_terminal = xs[N] - trajectory.get_reference(index + N);
        cost += (dx_terminal.transpose() * Qf * dx_terminal)(0, 0);
        return cost;
    }
};
 
// 主函数
int main() {
    Vehicle vehicle;  // 初始化车辆
    Trajectory trajectory;  // 初始化轨迹
    iLQRController controller(50, 10, 0.1);  // 初始化iLQR控制器
    double dt = 0.1;
    std::vector<double> x_history, y_history;  // 记录车辆轨迹
 
    for (size_t t = 0; t < trajectory.cx.size() - controller.N - 1; ++t) {
        VectorXd u_opt = controller.ilqr(vehicle, trajectory, t);  // 获取最优控制输入
        vehicle.update(u_opt[0], u_opt[1], dt);  // 更新车辆状态
        x_history.push_back(vehicle.x);  // 记录x坐标
        y_history.push_back(vehicle.y);  // 记录y坐标
 
        // 绘制车辆轨迹和参考轨迹
        plt::clf();
        plt::named_plot("Reference Trajectory", trajectory.cx, trajectory.cy, "-r");
        plt::named_plot("Vehicle Trajectory", x_history, y_history, "-b");
        plt::legend();
        plt::xlim(0, 50);
        plt::ylim(-20, 25);
        plt::title("iLQR Trajectory Tracking");
        plt::xlabel("x [m]");
        plt::ylabel("y [m]");
        plt::grid(true);
        plt::pause(0.001);  // 短暂暂停以动态展示
    }
 
    plt::show();  // 展示最终图形
    return 0;
}