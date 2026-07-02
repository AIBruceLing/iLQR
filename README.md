# iLQR 轨迹跟踪控制器

基于迭代线性二次调节器 (iterative Linear Quadratic Regulator, iLQR) 的车辆轨迹跟踪控制器，使用 C++ 和 Eigen 实现，并通过 matplotlib-cpp 进行实时可视化。

## 算法概述

iLQR 是一种非线性最优控制算法，通过在当前轨迹附近反复线性化动力学模型并求解 LQR 子问题，迭代地优化控制序列，使系统状态尽可能跟踪参考轨迹。

### 核心思想

标准 LQR 只能处理线性系统，而 iLQR 通过以下策略将其扩展到非线性系统：

1. 在当前名义轨迹上**线性化**非线性动力学模型
2. 利用线性化后的模型执行 **Backward Pass**，求解最优反馈增益
3. 利用反馈增益执行 **Forward Pass**，更新控制序列与状态轨迹
4. **重复迭代**直到代价函数收敛

## 系统建模

### 状态空间

状态向量为 4 维：

$$\mathbf{x} = [x, \; y, \; \theta, \; v]^T$$

| 符号 | 含义 | 单位 |
|------|------|------|
| $x$ | 车辆在全局坐标系下的 x 位置 | m |
| $y$ | 车辆在全局坐标系下的 y 位置 | m |
| $\theta$ | 车辆航向角 | rad |
| $v$ | 车辆速度 | m/s |

### 控制输入

控制向量为 2 维：

$$\mathbf{u} = [a, \; \omega]^T$$

| 符号 | 含义 | 单位 |
|------|------|------|
| $a$ | 纵向加速度 | m/s^2 |
| $\omega$ | 航向角速度（转向速率） | rad/s |

### 运动学模型（离散时间）

```
x_{k+1}     = x_k     + v_k * cos(theta_k) * dt
y_{k+1}     = y_k     + v_k * sin(theta_k) * dt
theta_{k+1} = theta_k + omega_k * dt
v_{k+1}     = v_k     + a_k * dt
```

对应代码中的 `dynamics()` 函数和 `Vehicle::update()` 方法。

### 线性化

在名义状态 $(\mathbf{x}_k, \mathbf{u}_k)$ 处对动力学模型求雅可比矩阵：

**状态转移矩阵** $A = \frac{\partial f}{\partial x}$：

```
A = I + dt * | 0  0  -v*sin(theta)  cos(theta) |
             | 0  0   v*cos(theta)  sin(theta) |
             | 0  0   0             0           |
             | 0  0   0             0           |
```

**控制输入矩阵** $B = \frac{\partial f}{\partial u}$：

```
B = dt * | 0  0 |
         | 0  0 |
         | 0  1 |
         | 1  0 |
```

对应代码中的 `linearize_dynamics()` 函数。

## iLQR 算法流程

### 代价函数

算法最小化如下二次代价函数：

$$J = \sum_{k=0}^{N-1} \left[ (\mathbf{x}_k - \mathbf{x}_k^{ref})^T Q (\mathbf{x}_k - \mathbf{x}_k^{ref}) + \mathbf{u}_k^T R \mathbf{u}_k \right] + (\mathbf{x}_N - \mathbf{x}_N^{ref})^T Q_f (\mathbf{x}_N - \mathbf{x}_N^{ref})$$

其中权重矩阵：

| 矩阵 | 含义 | 取值 |
|------|------|------|
| $Q$ | 状态跟踪代价 | diag(1.0, 1.0, 0.5, 0.1) |
| $R$ | 控制输入代价 | diag(0.1, 0.1) |
| $Q_f$ | 终端状态代价 | 10 * Q |

- $Q$ 中 x、y 位置权重最高（1.0），航向角次之（0.5），速度最低（0.1）
- $R$ 惩罚过大的控制输入，避免加速度和转向速率过大
- $Q_f$ 放大终端代价，保证规划末端也接近参考轨迹

### 算法步骤

```
输入: 当前车辆状态, 参考轨迹, 起始索引
输出: 当前时刻的最优控制量 u*

1. 初始化控制序列 u[0..N-1] = 0，状态序列 x[0] = 当前状态

2. FOR iter = 0 to max_iter-1:

   ===== Forward Pass (前向传播) =====
   3. 用当前控制序列前向模拟动力学模型:
      FOR i = 0 to N-1:
          x[i+1] = dynamics(x[i], u[i], dt)

   ===== Backward Pass (反向传播) =====
   4. 初始化终端值函数梯度和 Hessian:
      Vx  = Qf * (x[N] - x_ref[N])
      Vxx = Qf

   5. 从 i = N-1 反向递推到 i = 0:
      a. 线性化动力学: A, B = linearize_dynamics(x[i], u[i])
      b. 计算阶段代价的梯度和 Hessian:
         lx  = Q * (x[i] - x_ref[i])
         lu  = R * u[i]
         lxx = Q, luu = R, lux = 0
      c. 计算 Q-function 的二阶展开:
         Qx  = lx  + A^T * Vx
         Qu  = lu  + B^T * Vx
         Qxx = lxx + A^T * Vxx * A
         Quu = luu + B^T * Vxx * B
         Qux = lux + B^T * Vxx * A
      d. 计算最优反馈增益:
         k[i] = -Quu^{-1} * Qu       (前馈项)
         K[i] = -Quu^{-1} * Qux      (反馈项)
      e. 更新值函数:
         Vx  = Qx + K^T * Quu * k + K^T * Qu + Qux^T * k
         Vxx = Qxx + K^T * Quu * K + K^T * Qux + Qux^T * K

   ===== 控制更新 =====
   6. 使用反馈增益更新控制序列:
      FOR i = 0 to N-1:
          du = k[i] + K[i] * (x_new[i] - x[i])
          u_new[i] = u[i] + du
          x_new[i+1] = dynamics(x_new[i], u_new[i], dt)

   7. 检查收敛: 若总代价 < 1e-6 则终止迭代

3. 返回 u[0] 作为当前时刻的最优控制量
```

## 参考轨迹

参考轨迹为一条正弦调制曲线：

```
x(t) = t,              t ∈ [0, 50), 步长 0.1
y(t) = sin(t / 5) * t / 2
```

目标速度统一设为 3.0 m/s，航向角由相邻点差分计算 `atan2(dy, dx)` 得出。

## 代码结构

```
src/
├── iLQR.cpp          # 主程序（包含所有类和 main 函数）
├── matplotlibcpp.h   # matplotlib C++ 绑定（第三方头文件）
└── CMakeLists.txt    # 构建配置
```

### 类说明

| 类名 | 职责 |
|------|------|
| `Vehicle` | 车辆运动学模型，维护状态 [x, y, theta, v] 并提供 `update()` 方法 |
| `Trajectory` | 生成并存储参考轨迹（位置、航向角、目标速度），提供 `get_reference()` 查询 |
| `iLQRController` | iLQR 算法核心，包含动力学模型、线性化、代价计算和迭代求解 |

### 关键方法

| 方法 | 所属类 | 功能 |
|------|--------|------|
| `dynamics()` | iLQRController | 非线性运动学前向递推 |
| `linearize_dynamics()` | iLQRController | 在给定状态处计算雅可比矩阵 A、B |
| `ilqr()` | iLQRController | iLQR 主循环：前向模拟 -> 反向求增益 -> 更新控制 |
| `compute_total_cost()` | iLQRController | 计算完整轨迹的二次代价函数值 |
| `update()` | Vehicle | 使用欧拉积分更新车辆状态 |
| `get_reference()` | Trajectory | 返回指定索引处的参考状态向量 |

### 关键参数

| 参数 | 默认值 | 含义 |
|------|--------|------|
| `N` | 50 | 预测步数（控制时域长度） |
| `max_iter` | 10 | iLQR 最大迭代次数 |
| `dt` | 0.1 s | 离散化时间步长 |

## 依赖

- **Eigen3** - 线性代数库（矩阵运算）
- **Python3** + **NumPy** - matplotlib-cpp 可视化所需
- **CMake** >= 3.8

## 编译与运行

```bash
cd src
mkdir -p build && cd build
cmake ..
make
./iLQR
```

运行后会实时显示参考轨迹（红色）与车辆实际轨迹（蓝色）的跟踪效果图。

## 算法特点与局限

### 特点

- 利用二阶信息（Hessian）加速收敛，通常 5-10 次迭代即可获得良好的控制序列
- 每个时间步仅应用 `u[0]`，实现滚动时域控制（类 MPC 策略）
- 同时输出前馈项 k 和反馈增益 K，兼顾轨迹优化与扰动鲁棒性

### 局限

- 未实现正则化（regularization）：当 $Q_{uu}$ 接近奇异时，直接求逆可能数值不稳定
- 未实现线搜索（line search）：Forward Pass 中缺少步长缩放，大偏差下可能不收敛
- 未施加控制量约束：加速度和转向速率无上下限，实际应用中需增加约束处理
- 每步重新从零初始化控制序列，未利用上一步的解作为暖启动（warm start）
