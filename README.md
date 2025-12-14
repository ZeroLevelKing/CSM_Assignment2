# CSM Assignment 2 - MATLAB 有限元分析项目

这是一个基于 MATLAB 的有限元分析 (FEA) 项目，包含两个独立的计算固体力学 (CSM) 问题求解器。该项目旨在通过编程实现有限元方法，解决二维线弹性问题和 Mindlin 板弯曲问题。

## 项目概览

本项目包含两个主要部分：
1.  **Problem 1**: 二维线性弹性力学分析 (平面应力/应变)。
2.  **Problem 2**: Mindlin 板弯曲分析 (考虑剪切变形的中厚板)。

每个问题都拥有独立的源代码目录 (`src/`) 和结果输出目录 (`result/`)。

---

## Problem 1: 二维线性弹性分析 (2D Linear Elasticity)

该求解器用于模拟二维平面结构在载荷作用下的位移和应力分布。

### 主要功能
*   **分析类型**: 二维平面应力 (Plane Stress)。
*   **单元类型**:
    *   **Tri3**: 3节点常应变三角形单元 (CST)。
    *   **Quad4**: 4节点双线性四边形单元。
*   **可视化**:
    *   位移云图。
    *   Von Mises 应力云图。
    *   变形前后网格对比。

### 如何运行
1.  打开 MATLAB 并进入 `problem1/src/` 目录。
2.  打开 `main.m` 文件。
3.  (可选) 修改参数：
    *   设置 `flag = 1` 使用三角形单元，或 `flag = 2` 使用四边形单元。
    *   在 `generate_mesh.m` 中修改 `nx`, `ny` 以调整网格密度。
4.  运行 `main.m`。
5.  结果将保存在 `problem1/result/` 目录下。

---

## Problem 2: Mindlin 板弯曲分析 (Mindlin Plate Bending)

该求解器基于 Mindlin-Reissner 板理论 (一阶剪切变形理论)，适用于求解中厚板的弯曲问题。

### 主要功能
*   **物理模型**: Mindlin 板理论 (包含剪切效应)。
*   **自由度**: 每个节点 3 个自由度 ($w, \theta_x, \theta_y$)。
*   **单元技术**:
    *   支持 Tri3 和 Quad4 单元。
    *   **剪切自锁控制**: 对剪切刚度矩阵采用**减缩积分 (Reduced Integration)** 策略 (Quad4 单元使用 1x1 高斯积分)，以防止剪切自锁。
*   **载荷**: 高斯分布的横向载荷 $q(x,y) = C \cdot e^{-R^2/2\sigma^2}$ (通过数值积分计算等效节点力)。
*   **边界条件**: 四边简支 (Pinned edges, $w=0$)。

### 如何运行
1.  打开 MATLAB 并进入 `problem2/src/` 目录。
2.  运行 `main.m`。
3.  程序将自动执行**网格收敛性研究**：
    *   依次运行 Tri3 和 Quad4 单元。
    *   自动测试不同的网格密度 (如 5x5, 10x10, 50x50)。
4.  结果将保存在 `problem2/result/` 目录下，按案例分类存储 (例如 `result/Quad4_50x50/`)。

---

## 项目结构

```
CSM_Assignment2/
├── README.md                # 项目说明文档
├── Assignment_discription/  # 作业详细说明和图片
├── problem1/                # 问题 1 文件夹
│   ├── src/                 # 源代码 (main.m, K_matrix.m 等)
│   └── result/              # 结果输出
└── problem2/                # 问题 2 文件夹
    ├── src/                 # 源代码
    └── result/              # 结果输出
```

### 通用代码架构
两个问题的代码结构高度一致，遵循标准的 FEM 流程：
*   `main.m`: 主程序入口，控制流程。
*   `generate_mesh.m`: 生成节点坐标和单元连接关系。
*   `K_matrix.m`: 组装全局刚度矩阵。
*   `F_vector.m`: 组装全局载荷向量。
*   `Boundary_conditions.m`: 定义边界条件。
*   `Enforce_BC.m`: 施加边界条件并求解方程组。
*   `plot_results.m`: 后处理与绘图。

## 环境要求
*   MATLAB (推荐 R2018b 或更高版本)。
