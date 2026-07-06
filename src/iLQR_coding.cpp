#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class VehicleModel{
    public:
    double x_, y_, theta_, v_;//转态量

    // 车辆构造函数，初始化车辆的状态
    VehicleModel(double x_init = 0.0, double y_init = 0.0, double theta_init = 0.0, double v_init = 0.0)
        : x_(x_init), y_(y_init), theta_(theta_init), v_(v_init) {}
    
    void UpdateVehicleState(double a, double w, double dt){
        x_ += v_*cos(theta_)*dt;
        y_ += v_*sin(theta_)*dt;
        theta_ += w*dt;
        v_ += a*dt;
    }  
};


class Trajectory{
    public:
    std::vector<double> cx_, cy_, theta_, v_;

    Trajectory(){
        for(double i = 0;i < 50; i += 0.1){
            cx_.push_back(i);
            cy_.push_back(sin(i/5.0)*i/2.0);
        }

        for(size_t i = 0; i < cx_.size(); i++){
            theta_.push_back(atan2(cy_[i+1] - cy_[i],cx_[i+1]-cx_[i]));
        }
        theta_.push_back(theta_.back());
        v_.assign(cx_.size(), 3.0);
    }

    VectorXd get_reference(size_t index) const{
        VectorXd ref(4);
        ref << cx_[index], cy_[index], theta_[index], v_[index];
        return ref;
    }
};


class iLQRController{
    public:
    int N;
    int max_iter;
    double dt;
    MatrixXd Q, R, Qf;

    iLQRController(size_t N_input = 50, size_t max_iter_input = 10, double dt_input = 0.1)
        : N(N_input), max_iter(max_iter_input), dt(dt_input){
            Q = MatrixXd::Zero(4, 4);
            Q(0,0) = 1.0;
            Q(1,1) = 1.0;
            Q(2,2) = 0.5;
            Q(3,3) = 0.1;
            R = MatrixXd::Identity(2,2)*0.1;
            Qf = Q * 10;
        }


        VectorXd iLQR(VehicleModel& vehicle, const Trajectory& trajectory, size_t index){ 
            const size_t state_dim = 4;
            const size_t control_dim = 2;
            std::vector<VectorXd> state_sequence(N + 1, VectorXd::Zero(state_dim));
            std::vector<VectorXd> control_sequence(N, VectorXd::Zero(control_dim));

            //车辆初始状态
            state_sequence[0] << vehicle.x_, vehicle.y_, vehicle.theta_, vehicle.v_;
            double pre_cost = 0.0;

            //迭代更新控制输入和状态
            for(size_t iter = 0; iter < max_iter; iter++){
                // 前向传播，rollout，得到标称轨迹state_sequence（Nominal Trajectory）
                for(int i = 0; i < N; i++){
                    state_sequence[i + 1] = dynamics(state_sequence[i],control_sequence[i], dt);
                }

                // 反向传播，backward
                /*
                (状态)值函数是state Value，只和状态值有关
                这里是反向传播（Backward Pass）起点时的初始化计算

                因为最后一步没有控制量了，所以最后一步(终端)的值函数等于代价函数
                终端代价函数（Terminal Cost）的一阶导数（梯度）和二阶导数（海森矩阵/曲率）
                */
                VectorXd Vx = Qf * (state_sequence[N] - trajectory.get_reference(index + N)); // 终端值函数的梯度 = Qf * dx
                MatrixXd Vxx = Qf;                                                            // 终端值函数的曲率 = 终端状态权重矩阵

                // 初始化前馈向量数组和反馈矩阵数组
                std::vector<VectorXd> k_control_list(N, VectorXd::Zero(control_dim));               
                std::vector<MatrixXd> k_feedback_list(N, MatrixXd::Zero(control_dim, state_dim));    
                 
                //得到每一个时间步的前馈向量d_k和反馈矩阵K_k，并更新价值函数梯度Vx和Vxx
                for(int i = N - 1; i >= 0; i--){
                    VectorXd state_ref = trajectory.get_reference(index + i);
                    std::pair<MatrixXd,MatrixXd> linearized = linearize_dynamics(state_sequence[i],
                                                                                 control_sequence[i], dt);
                    MatrixXd A = linearized.first;
                    MatrixXd B = linearized.second;

                    //计算二次型代价梯度和hessian矩阵
                    VectorXd lx = Q * (state_sequence[i] - state_ref);
                    VectorXd lu = R * control_sequence[i];                // 这里的参考轨迹中的参考控制量为0，所以du[i] = control_sequence[i]
                    MatrixXd lxx = Q; 
                    MatrixXd luu = R;
                    MatrixXd lux = MatrixXd::Zero(control_dim,state_dim); // 先对 x 求偏导再对 u 求偏导，结果为0

                    //计算Q函数(action value)的二阶泰勒展开系数——5个偏导：
                    VectorXd Qx = lx + A.transpose()*Vx;
                    VectorXd Qu = lu + B.transpose()*Vx;
                    MatrixXd Qxx = lxx + A.transpose()*Vxx*A;
                    MatrixXd Quu = luu + B.transpose()*Vxx*B;
                    MatrixXd Qux = lux + B.transpose()*Vxx*A;

                    //计算控制增量和反馈矩阵
                    // 最优控制增量：du = -(d_k + K_k * dx)
                    MatrixXd Quu_inv = Quu.inverse();
                    VectorXd k_control = -Quu_inv*Qu;     // 前馈向量  d_k = -Quu?? Qu
                    MatrixXd K_feedback = -Quu_inv*Qux;   // 反馈矩阵  K_k = -Quu?? Qux

                                           
                    k_control_list[i] = k_control;
                    k_feedback_list[i] = K_feedback;

                    //更新价值函数
                    Vx = Qx + K_feedback.transpose()*Quu*k_control + 
                         K_feedback.transpose()*Qu+Qux.transpose()*k_control;
                    Vxx = Qxx + K_feedback.transpose() * Quu * K_feedback + 
                         K_feedback.transpose() * Qux + Qux.transpose() * K_feedback;
                }


                //基于上一步得到的反馈和前馈更新每个时间步(轨迹步)的控制序列，并进行前向rollout模拟
                std::vector<VectorXd> state_sequence_new(N + 1, state_sequence[0]);
                std::vector<VectorXd> control_sequence_new;
                for(int i = 0; i < N; i++){
                    // 最优控制增量：du = -(d_k + K_k * dx)
                    VectorXd du = k_control_list[i] + k_feedback_list[i]*(state_sequence_new[i] - state_sequence[i]);
                    control_sequence_new.push_back(control_sequence[i] + du);
                    state_sequence_new[i+1] = dynamics(state_sequence_new[i], control_sequence_new.back(),dt);
                }
                state_sequence = state_sequence_new;
                control_sequence = control_sequence_new;

                //判断是否收敛
                double cur_cost = compute_total_cost(state_sequence,control_sequence,trajectory,index);
                double cost_error = std::fabs(cur_cost - pre_cost);
                std::cout<<"Iteration "<< iter << ",cost_error:" << cost_error << std::endl;
                if(cost_error < 1e-6 ){
                    break;
                }
                pre_cost = cur_cost;
            }
            std::cout<<"The best control: a is "<<control_sequence[0][0]<<", angle is "<<control_sequence[0][1]<<std::endl;

            return control_sequence[0];//返回控制序列的第一个控制量
        }
        

        //车辆运动学模型，这里是单轮模型
        VectorXd dynamics(const VectorXd& current_state, const VectorXd& current_control, double dt){
            VectorXd next_state = current_state;
            next_state[0] += current_state[3]*cos(current_state[2])*dt;
            next_state[1] += current_state[3]*sin(current_state[2])*dt;
            next_state[2] += current_control[1]*dt;
            next_state[3] += current_control[0]*dt;
            return next_state;
        }

        
        //对非线性运动学模型进行线性化，得到线性误差状态空间方程
        /*
        状态雅可比矩阵： A = [1 0 -v*sin(theta) cos(theta)
                           0 1  v*cos(theta) sin(theta)
                           0 0        1          0
                           0 0        0          1     ]
        
        
        状态雅可比矩阵： B = [0 0
                           0 0
                           0 dt
                           dt 0]
        */
       std::pair<MatrixXd, MatrixXd> linearize_dynamics(const VectorXd& state, const VectorXd& control, double dt){
            MatrixXd A = MatrixXd::Identity(4,4);
            A(0,2) = -state[3]*sin(state[2])*dt;
            A(0,3) = cos(state[2])*dt;
            A(1,2) = state[3]*cos(state[2])*dt;
            A(1,3) = sin(state[2])*dt;

            MatrixXd B = MatrixXd::Zero(4,2);
            B(2,1) = dt;
            B(3,0) = dt;
            return {A, B};
       }


       // 计算total cost
       // 注意这里的trajectory是全局参考轨迹，不是局部参考轨迹
       // 注意：默认“最完美、最省能”的参考控制量是 0，所以这里计算cost的时候直接用的绝对量
       double compute_total_cost(const std::vector<VectorXd>& state_sequence, 
                                 const std::vector<VectorXd>& control_sequence,
                                 const Trajectory& trajectory, size_t index){
            double total_cost = 0.0;
            for(int i = 0; i < N; ++i){
                VectorXd state_error = state_sequence[i] - trajectory.get_reference(index + i);

                //(dx.transpose() * Q * dx)得到是是1X1的矩阵，不是标量，所以要取（0，0）
                total_cost += (state_error.transpose()*Q*state_error).value() + 
                        (control_sequence[i].transpose()*R*control_sequence[i]).value();
            }
            VectorXd terminal_state_error = state_sequence[N] - trajectory.get_reference(index+N);
            total_cost += (terminal_state_error.transpose()*Qf*terminal_state_error).value();

            return total_cost;
        }

};

int main(){
    VehicleModel vehicle;
    Trajectory trajectory;
    iLQRController controller(50, 10, 0.1);
    double dt = 0.1;
    std::vector<double> x_history, y_history;

    for(size_t t = 0; t < trajectory.cx_.size() - controller.N - 1; t++){
        VectorXd u_opt = controller.iLQR(vehicle, trajectory, t);
        vehicle.UpdateVehicleState(u_opt[0], u_opt[1], dt);
        x_history.push_back(vehicle.x_);
        y_history.push_back(vehicle.y_);

       // 绘制车辆轨迹和参考轨迹
        plt::clf();
        plt::named_plot("Reference Trajectory", trajectory.cx_, trajectory.cy_, "-r");
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
    plt::show();
    return 0;
}
