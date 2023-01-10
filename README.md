# 快速傅里叶变换的逻辑实现
傅里叶变换可以将时域信号转化为频域信号，从而能够得到信号的频谱结构和变化规律。傅里叶变换最基础的部分是离散傅里叶变换（Discrete Fourier Transform，DFT），是不可或缺的一部分。离散傅里叶变换是一种将傅里叶变换离散化的一种方法，即将时域中的连续信号采样后得到的信号转化为频域中的信号。由于离散傅里叶变换的运算量太大，因此直接使用 DFT 算法进行实时的信号处理就会花费不必要的消耗，快速傅里叶变换应运而生。

快速傅里叶变换（Fast Fourier Transform，FFT）是离散傅里叶变换（Discrete Fourier Transform，DFT）的快速算法，快速傅里叶变换降低了运算量，加快了 DFT 的运算速度的方法是将长序列的 DFT 分解成为短序列的 DFT，从而为数字信号处理技术实现多种信号实时处理提供了充足的条件，使得用户更加容易实现和应用数字信号处理，故 FFT 是数字信号处理领域广泛使用的最有效的算法之一。因其计算量小的优势，FFT 被广泛的应用于语音识别、图像处理、无线通信、频谱分析、雷达处理及数字通信等领域中。

FFT 的基本思想是将一个信号分解为若干个正弦波和余弦波的加权和，从而得到信号的频谱。通过这种方式，可以在信号中提取出特定频率分量，并且可以在频域中进行操作。下面将介绍傅里叶级数以及傅里叶变换。

## 连续周期信号的傅里叶级数
### 正弦余弦型傅里叶级数（三角型傅里叶级数）
以 $T$ 为周期的连续周期信号 $f_T(t)$ ，如果满足狄利克雷条件：
  1. 在一个周期内只有有限个间断点；
  2. 在一个周期内有有限个极值点；
  3. 在一个周期内函数绝对可积，即 $\int_{t_0}^{t_0+T}|f_T(t)|\\ dt < \infty$

则可展开成完备的正交三角函数集 $\{\cdots,cos(n\omega_0 t),\cdots,sin(n\omega_0 t),\cdots \}(n=0,1,2,\cdots)$ 的线性组合
$$f_T(t)=\frac{a_0}{2}+\sum_{n=1}^{\infty}[a_n cos(n\omega_0 t)+b_n sin(n\omega_0 t)],\\ n=0,1,2,\cdots \tag{1.1}$$

上式称为正弦余弦形式的三角型傅里叶级数。其中， $\omega_0=\frac{2\pi}{T}$ 是基波角频率， $a_0,a_n,b_n$ 是傅里叶级数系数。根据正弦余弦函数的正交条件，可得
$$a_n=\frac{2}{T}\int_{t_0}^{t_0+T}f_T(t)cos(n\omega_0 t)\\ dt \tag{1.2}$$
$$b_n=\frac{2}{T}\int_{t_0}^{t_0+T}f_T(t)sin(n\omega_0 t)\\ dt \tag{1.3}$$

显然，傅里叶系数 $a_n$ 和 $b_n$ 都是 $n\omega_0$ 的函数， $a_n$ 是 $n$ 的偶函数， $b_n$ 是 $n$ 的奇函数。

### 指数型傅里叶级数
三角型傅里叶级数物理含义明确，但是运算不便，因此常用指数型的傅里叶级数。

以 $T$ 为周期的周期信号 $f_T(t)$ 的指数型傅里叶级数定义为
$$f_T(t) \overset{def}{=} \sum_{n=-\infty}^{\infty}F_n e^{jn\omega_0 t},\\ \omega_0=\frac{2\pi}{T} \tag{1.4}$$

其中
$$F_n \overset{def}{=} \frac{1}{T} \int_{-\frac{T}{2}}^{\frac{T}{2}}f_T(t) e^{-jn\omega_0 t}\\ dt,\\ n=0,\pm 1,\pm 2, \cdots \tag{1.5}$$
称为指数型傅里叶级数的系数，又称为傅里叶复系数。

傅里叶复系数除了用式 $(1.5)$ 求得，也可以根据三角型傅里叶系数导出，根据欧拉公式，有
$$cos(n\omega_0 t)=\frac{e^{jn\omega_0 t}+e^{-jn\omega_0 t}}{2},\\ sin(n\omega_0 t)=\frac{e^{jn\omega_0 t}-e^{-jn\omega_0 t}}{2j} \tag{1.6}$$

代入式 $(1.1)$ ，得

$$\begin{align*}
  f_T(t) &= \frac{a_0}{2}+\sum_{n=1}^{\infty}a_n \frac{e^{jn\omega_0 t}+e^{-jn\omega_0 t}}{2}+\sum_{n=1}^{\infty}b_n\frac{e^{jn\omega_0 t}-e^{-jn\omega_0 t}}{2j} \\
  &= \frac{a_0}{2}+\sum_{n=1}^{\infty}\frac{a_n-jb_n}{2}e^{jn\omega_0 t}+\sum_{n=1}^{\infty}\frac{a_n+jb_n}{2}e^{-jn\omega_0 t} \\
  &= \frac{a_0}{2}+\sum_{n=1}^{\infty}F_n e^{jn\omega_0 t}+\sum_{n=1}^{\infty}F_{-n} e^{-jn\omega_0 t} \\
  &= F_0+\sum_{n=1}^{\infty}F_n e^{jn\omega_0 t}+\sum_{n=-\infty}^{-1}F_n e^{jn\omega_0 t} \\
  &= \sum_{n=-\infty}^{\infty}F_n e^{jn\omega_0 t}
\end{align*}$$

即
$$f_T(t)=\sum_{n=-\infty}^{\infty}F_n e^{jn\omega_0 t}$$
其中

$$\begin{cases}
  F_0=\frac{a_0}{2}=\frac{A_0}{2} \\
  F_n=\frac{a_n-jb_n}{2}=\frac{A_n}{2}e^{j\varphi_n} \\
  F_{-n}=F_{n}^{*}=\frac{a_n+jb_n}{2}
\end{cases}$$

## 连续周期信号的频谱
周期信号是一系列相互正交的正弦分量 $A_n cos(n\omega_0 t+\varphi_n)$ 或复指数分量 $F_n e^{jn\omega_0 t}=|F_n|e^{j\varphi_n}e^{jn\omega_0 t}$ 的线性组合。对实信号而言，复指数信号分量 $F_n e^{jn\omega_0 t}$ 与 $F_{-n}e^{-jn\omega_0 t}$ 是成对出现的。傅里叶系数的幅度 $|F_n|$ 或 $A_n$ 随角频率 $n\omega_0$ 变化的规律称为信号的**幅度频谱**，简称幅度谱；傅里叶系数的相位 $\varphi_n$ 随角频率 $n\omega_n$ 变化的规律称为信号的**相位频谱**，简称相位谱；幅度谱和相位谱总称为信号的**幅相谱**。

### 三角频谱（单边谱）
三角频谱是指余弦形式的傅里叶级数的振幅 $A_n$ 随角频率 $n\omega_0$ 变化的规律，称为振幅频谱；余弦形式的傅里叶级数的初相 $\varphi_n$ 随角频率 $n\omega_0$ 变化的规律，称为**相位频谱**。由于三角型傅里叶级数总有 $n\ge 0$ ，谱线只出现在 $A_n (-n\omega_0)$ 或 $\varphi_n (-n \omega_0)$ 平面的右半平面，故称为**单边频谱**。

### 指数频谱（双边谱）
指数频谱是指傅里叶复系数 $F_n$ 的幅度（模） $|F_n|$ 随角频率 $n\omega_0$ 变化的规律，以及它的相位 $\varphi_n$ 随角频率 $n\omega_0$ 变化的规律，前者称振幅频谱，后者称相位频谱。由于指数型傅里叶级数总有 $-\infty <0< \infty$ ，所以 $n\omega_0$ 的取值是正负整倍基波角频率，故称作双边频谱。

## 连续非周期信号的频谱——傅里叶变换
当周期 $T$ 趋于无穷大时，周期信号 $f_T(t)$ 变为非周期信号 $f(t)$ 。从频谱分析的观点来看，当 $T$ 增加时，基波频率 $\omega_0$ 变小，离散谱线变密，频谱幅度变小，但频谱包络线的形状保持不变。在极限情况下，当周期 $T$ 趋于无穷大时，其谱线间隔与幅度将会趋于无穷小，原来由许多谱线组成的周期信号的离散频谱就会连成一片，称为面频谱，并且从 $|F_n|(-n\omega_0)$ 平面消失。这时，非周期信号的 $|F_n|$ 成为无限小量而存在却再也看不见，傅里叶级数在此没有任何意义。为了描述非周期信号的频谱特性，引入频谱密度函数。

定义
$$F(j\omega) \overset{def}{=} \lim_{\Delta \to 0}\frac{F_n}{\Delta F}=\lim_{T \to \infty}TF_n \tag{2.1}$$
为频谱密度函数。 $F(j\omega)$ 表示了单位频带 $\Delta F$ 上的复振幅 $F_n$ 分布状况。 $F(j\omega)$ 可理解为一种频谱密度，它表达了非周期信号在任意频率 $\omega$ 处的频谱密度函数。

根据周期信号傅里叶级数定义式，有
$$TF_n=\int_{-\frac{T}{2}}^{\frac{T}{2}}f_T(t)e^{-jn\omega_0 t}\\ dt \tag{2.2}$$
$$f_T(t)=\sum_{n=-\infty}^{\infty}F_n e^{jn\omega_0 t} \tag{2.3}$$
考虑到当周期 $T$ 趋近于无穷大时，上式中各变量将作如下变化： $\omega_0=\frac{2\pi}{T}$ 趋于无穷小， $\omega_0 \to \Delta \omega \to d\omega$ ，取其为 $d\omega$ ；离散变量 $n\omega_0$ 则变为连续变量 $\omega$ ； $F_n$ 变为 $\frac{d\omega}{2\pi}F(j\omega)$ ；离散和变为连续和。于是式 $(2.2)$ 变为
$$F(j\omega)=\lim_{T \to \infty}TF_n=\lim_{T \to \infty}\frac{T}{T}\int_{-\frac{T}{2}}^{\frac{T}{2}}f_T(t)e^{-jn\omega_0 t}\\ dt=\int_{-\infty}^{\infty}f(t)e^{-j\omega t}\\ dt \tag{2.4}$$
式 $(2.3)$ 变为

$$\begin{align*}
f(t) &= \lim_{T \to \infty}f_T(t)=\lim_{T \to \infty}\sum_{n=-\infty}^{\infty}F_n e^{jn\omega_0 t}\\
&= \lim_{T \to \infty}\sum_{n=-\infty}^{\infty}\frac{1}{2\pi}F(j\omega)e^{jn\omega_0 t}\Delta \omega=\frac{1}{2\pi}\int_{-\infty}^{\infty}F(j\omega)e^{j\omega t}\\ d\omega
\end{align*}$$

即
$$F(j\omega) \overset{def}{=} \int_{-\infty}^{\infty}f(t)e^{-j\omega t}\\ dt \tag{2.5}$$
$$f(t) \overset{def}{=}\frac{1}{2\pi}\int_{-\infty}^{\infty}F(j\omega)e^{j\omega t}\\ d\omega \tag{2.6}$$

式 $(2.5)$ 和 $(2.6)$ 称为傅里叶变换对，其中式 $(2.5)$ 称为傅里叶正变换（FT），简称傅里叶变换；而式 $(2.6)$ 称为傅里叶反变换（IFT）。 $F(j\omega)$ 是 $f(t)$ 的**频谱密度函数**即频谱函数，而 $f(t)$ 是 $F(j\omega)$ 的**原函数**。

## 离散傅里叶变换
设 $f(k)$ 是一个有限长序列，其长度为 $N$ ，即在区间 $0 \le k \le N-1$ 以外， $f(k)$ 为 0 。将 $f(k)$ 以周期 $N$ 延拓而成的周期序列记为 $f_N(k)$ ，则有
$$f_N(k)=\sum_{m=-\infty}^{\infty}f_1(k-mN),\\ m 为整数 \tag{3.1}$$
其中
$$f_1(k)=f(k)g_N(k) \tag{3.2}$$

周期序列 $f_N(k)$ 的离散时间傅里叶级数表达式为
$$F_n=\frac{1}{N}\sum_{k=(N)}f_N(k)e^{-jn\Omega_0k} \tag{3.3}$$
$$f_N(k)=\sum_{n=(N)}F_ne^{jn\Omega_0k} \tag{3.4}$$

如果将 $NF_n$ 表示成 $F(n)$ ，且在 $[0,N-1]$ 内 $f(k)=f_N(k)$ 并令 $W=e^{-j\Omega_0}=e^{-j\frac{2\pi}{N}}$ ，则上式可改写为
$$F(n) \overset{def}{=}\sum_{k=0}^{N-1}f(k)W^{nk} \tag{3.5}$$
$$f(k) \overset{def}{=}\frac{1}{N}\sum_{n=0}^{N-1}F(n)W^{-nk} \tag{3.6}$$

式 $(3.5)$ 和 $(3.6)$ 所定义的变换关系称为离散傅里叶变换（DFT）。

DFT 表明，时域 $N$ 点有限长序列 $f(k)$ 可以变换为频域的 $N$ 点有限长序列 $F(n)$ 。

式 $(3.5)$ 和 $(3.6)$ 写为矩阵形式，有

$$\begin{bmatrix}
F(0)\\
F(1)\\
\vdots \\
F(N-1)
\end{bmatrix}
=\begin{bmatrix}
W^0 & W^0 & W^0 & \cdots & W^0 \\
W^0 & W^{1\times 1} & W^{2\times 1} & \cdots & W^{(N-1)\times 1} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
W^0 & W^{1\times (N-1)} & W^{2\times (N-1)} & \cdots & W^{(N-1)\times (N-1)}
\end{bmatrix}
\begin{bmatrix}
f(0) \\
f(1) \\
\vdots \\
f(N-1)
\end{bmatrix} \tag{3.7}$$

$$\begin{bmatrix}
f(0)\\
f(1)\\
\vdots \\
f(N-1)
\end{bmatrix}
=\frac{1}{N}
\begin{bmatrix}
W^0 & W^0 & W^0 & \cdots & W^0 \\
W^0 & W^{-1\times 1} & W^{-2\times 1} & \cdots & W^{-(N-1)\times 1} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
W^0 & W^{-1\times (N-1)} & W^{-2\times (N-1)} & \cdots & W^{-(N-1)\times (N-1)}
\end{bmatrix}
\begin{bmatrix}
F(0) \\
F(1) \\
\vdots \\
F(N-1)
\end{bmatrix} \tag{3.8}$$

简记为
$$\mathbf{F}(n)=\mathbf{W}^{kn}\mathbf{f}(k) \tag{3.9}$$
$$\mathbf{f}(k)=\frac{1}{N}\mathbf{W}^{-kn}\mathbf{F}(n) \tag{3.10}$$

其中， $\mathbf{W}^{kn}$ 和 $\mathbf{W}^{-kn}$ 都是 $N\times N$ 对称矩阵。

## 快速傅里叶变换
