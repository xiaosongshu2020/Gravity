# 月球引力场计算项目 (Lunar Gravity Field Calculation)

![项目封面](cover.png)

这是一个用于计算月球引力场的Fortran项目，基于球谐函数模型实现高精度引力场计算。

## 项目概述

本项目提供了一个完整的月球引力场计算框架，包括：

- **引力场模型定义**：支持球谐系数加载和引力加速度计算
- **模块化设计**：采用面向对象风格的Fortran模块设计
- **高精度计算**：支持高达80阶80次的球谐展开
- **测试程序**：包含完整的测试用例验证计算准确性

## 文件结构

```
lunarGravity/
├── module_Gravity.f90      # 引力场计算核心模块
├── test_lunar_gravity.f90  # 测试程序
├── grail.txt              # 月球引力场系数数据文件
├── README.md              # 项目说明文档
└── .gitignore            # Git忽略文件配置
```

## 核心功能

### 引力场对象 (GravityField)

引力场对象包含以下属性：
- 引力场名称和中心天体
- GM值（引力常数 × 质量）
- 参考半径
- 最大阶数和次数
- 球谐系数矩阵 C 和 S
- 系数加载状态

### 主要方法

1. **初始化** (`initialize_gravity_field`)：创建引力场对象
2. **读取系数** (`read_gravity_coefficients`)：从文件加载球谐系数
3. **计算加速度** (`compute_gravity_acceleration`)：计算给定位置的引力加速度
4. **显示信息** (`get_gravity_field_info`)：输出引力场基本信息

## 编译和运行

### 编译方法

使用Fortran编译器编译项目：

```bash
# 编译模块
gfortran -c module_Gravity.f90 -o module_Gravity.o

# 编译测试程序
gfortran -c test_lunar_gravity.f90 -o test_lunar_gravity.o

# 链接生成可执行文件
gfortran module_Gravity.o test_lunar_gravity.o -o test_lunar_gravity.exe
```

### 运行测试

```bash
./test_lunar_gravity.exe
```

## 数据文件格式

`grail.txt` 文件格式：
- 第一行：参考半径、GM值、参考经度、参考纬度
- 后续行：阶数l、次数m、C系数、S系数、C不确定性、S不确定性

## 技术特点

1. **球谐函数计算**：使用完全归一化的缔合勒让德多项式
2. **数值稳定性**：采用递归算法确保计算稳定性
3. **模块化设计**：易于扩展支持其他天体引力场
4. **高精度**：支持高阶球谐展开，满足科研需求

## 应用场景

- 月球轨道动力学分析
- 月球探测器轨道设计
- 月球重力场科学研究
- 天体力学教学演示

## 依赖要求

- Fortran编译器（gfortran、ifort等）
- 基本的线性代数运算库

## 许可证

本项目采用开源许可证，具体信息请参考LICENSE文件。

## 贡献

欢迎提交Issue和Pull Request来改进项目。

## 联系方式

如有问题请联系项目维护者。
