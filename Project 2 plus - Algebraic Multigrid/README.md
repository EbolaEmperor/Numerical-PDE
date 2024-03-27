# 代数多重网格

实现了一个简单的代数多重网格，参考了William L. Briggs等人写的多重网格经典入门书籍《A Multigrid Tutorial, 2nd Edition》

给定稀疏矩阵A和方程右端向量b，使用如下代码调用求解器以求解Ax=b

```cpp
    amgSolver solver;
    solver.generateGrid(A);
    ColVector sol = solver.solve(b, "FMG", 20, 1e-12);
```

上面的`FMG`是Cycle的类型，可以选择`FMG`或`V`；`20`是最大迭代次数，`1e-12`是终止误差（以相对误差的无穷范数计），二者达成任何一个条件即终止迭代。

当网格较大时，请关闭系统的栈空间限制，否则会发生Segmentation Fault. 在Linux系统中，您可以运行命令：
 
```bash
    ulimit -s unlimited 
```

编译方式：假设你的文件名是`youCode.cpp`，使用如下命令编译：

```bash
g++ youCode.cpp src/*.cpp -Iinclude -O2 -O3 -Ofast
```
