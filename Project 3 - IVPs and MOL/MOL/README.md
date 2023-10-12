# 关于热方程的特别说明

我们为热方程专门提供了 Fourier Solver 和 DFT Solver。前者利用分离变量法后的Fourier级数展开，直接解析求解，但是只支持齐次边界条件。后者利用DFT将热方程转化为频域上若干个独立的线性ODE，但是只支持周期边界条件，另外由于DFT的实现基于FFT，所以空间划分数必须是2的幂次。

具体使用方法请见 `samples/heat-fourier.json` 与 `samples/heat-dft.json`