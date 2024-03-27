# 数值积分求解器

请确保您已经安装`gnuplot`，若没有安装，请通过命令：

```bash
    sudo apt install gnuplot-qt
```

来完成安装。然后您可以运行命令：

```bash
    make story
```

将会开始自动生成图片与报告（约需三分钟，含求解器的使用方法），但不会复现报告中的数值测试。若您想验证报告中数值测试的正确性，请运行：

```bash
    make run
```

它将会复现所有的测试。您可以通过命令

```bash
    make clean
```

来清除所有生成的文件（但求解器与报告将被保留）

# 热方程与平流方程

在项目根目录运行

```bash
make heat
```

来编译热方程求解器。它会在MOL文件夹中生成一个heat，使用

```bash
./heat xxx.json
```

来运行。其中`xxx.json`是输入文件名，可参考样例`MOL/samples/heat.json`。

在项目根目录运行

```bash
make advection
```

来编译平流方程求解器。它会在MOL文件夹中生成一个advection，使用

```bash
./advection xxx.json
```

来运行。其中`xxx.json`是输入文件名，可参考样例`MOL/samples/advection.json`。