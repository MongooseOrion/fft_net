# 使用 Verilog 编写一维快速傅里叶变换

## 基于矩阵运算的 FFT 算法

基于矩阵运算的FFT算法主要有以下几个步骤：

  1. 对输入序列进行递归分解，直到分解到长度为 1 的子序列；
  2. 计算旋转因子 $w_n = e^{-2\pi \frac{i}{n}}$；
  3. 计算长度为 $n/2$ 的偶数部分和奇数部分，分别进行 FFT 计算；
  4. 将偶数部分和奇数部分的 FFT 结果组合起来，得到 FFT 的结果。

基于矩阵运算的 FFT 算法相对于蝶形运算算法更易于理解，但是它的运算量和空间复杂度较高，因此不如蝶形运算算法快速和高效。

### 旋转因子

旋转因子是用于 DFT 计算中的一个复数系数，它的值为 $w_n = e^{-2πi/n}$ ，其中 $i$ 是虚数单位，$n$ 是 DFT 的长度。旋转因子在 FFT 算法中被频繁地使用。旋转因子的工作原理是通过将输入序列分解为偶数项和奇数项，并在计算它们的 DFT 时使用旋转因子 $w_n$ 进行合并。

假设有两个长度为 $n$ 的向量 $a$ 和 $b$ ，其中 $b$ 是 $a$ 的偶数项， $c$ 是 $a$ 的奇数项，则可以将 $a$ 表示为：

$$a = (c_0, c_1, ..., c_{n/2-1}, c_{n/2}, c_{n/2+1}, ..., c_{n-1})$$

将 $a$ 的 DFT 表示为 $X$ ，即：

$$X_k = \sum_{j=0}^{n-1} a_j e^{-2πij\frac{k}{n}}$$

其中， $k$ 表示 DFT 的系数， $j$ 表示 $a$ 的索引。将 $a$ 分解为偶数项和奇数项，则 $X_k$ 表示为：

$$X_k = \sum_{j=0}^{\frac{n}{2}-1} c_j \cdot e^{-2πij\frac{k}{n}} + e^{-2πi\frac{k}{n}} \sum_{j=0}^{\frac{n}{2}-1} c_{j+\frac{n}{2}} \cdot e^{-2πij\frac{k}{n}}$$

通过使用旋转因子 $w_n = e^{-2πi/n}$ ，可以将第二项中的旋转因子表示为 $w_n^k$ ，即：

$$X_k = \sum_{j=0}^{\frac{n}{2}-1} c_j \cdot e^{-2πijk/n} + w_n^k \sum_{j=0}^{\frac{n}{2}-1} c_{j+\frac{n}{2}} \cdot e^{-2πij\frac{k}{n}}$$

### 计算方法

下述为一段使用矩阵运算 FFT 的伪代码：
```
function FFT(x)
    n = length(x)
    if n == 1 then
        return x
    end if

    w_n = exp(-2 * pi * i / n)                  // 定义旋转因子
    w = 1
    x_even = FFT([x[0], x[2], ..., x[n-2]])     // 偶数项
    x_odd = FFT([x[1], x[3], ..., x[n-1]])      // 奇数项
    X = new array of size n
    for k from 0 to n/2 - 1 do
        X[k] = x_even[k] + w * x_odd[k]
        X[k + n/2] = x_even[k] - w * x_odd[k]
        w = w * w_n
    end for
    return X
end function
```

对于矩阵运算的 FFT，我们可以将输入的序列看成一个列向量，将傅里叶变换看成是一个矩阵乘法操作。那么，对于一个长度为 $N$ 的序列，FFT 可以表示为如下形式：

$$y = F_Nx$$

其中 $x$ 是 $N\times 1$ 的列向量， $y$ 也是 $N\times 1$ 的列向量，而 $F_N$ 则是一个 $N\times N$ 的矩阵，它的第 $k$ 行第 $n$ 列的元素为：

$$F_N[k,n]=e^{-j\frac{2\pi kn}{N}}$$

其中 $j$ 为虚数单位， $e$ 为自然对数的底数。上述公式是 DFT 的矩阵形式，而 FFT 可以看作是 DFT 的快速实现。

因此，对于一个长度为 $N$ 的序列 $x$，我们可以将其 FFT 的计算过程表示为如下形式：

$$y_n=\sum_{k=0}^{N-1}F_N[n,k]x_k$$

根据矩阵乘法的定义，我们可以进一步将 $y=F_Nx$ 表示为如下形式：

$$\begin{align*}
\vec{y}&= \begin{bmatrix} 
\vec{y_0} \ \vec{y_1} \ \cdots \ \vec{y}_{N-1} 
\end{bmatrix} \\
&= \begin{bmatrix} 
F_N[0,0] & F_N[0,1] & \cdots & F_N[0,N-1] \\
F_N[1,0] & F_N[1,1] & \cdots & F_N[1,N-1] \\ 
\vdots & \vdots & \ddots & \vdots \\ 
F_N[N-1,0] & F_N[N-1,1] & \cdots & F_N[N-1,N-1] 
\end{bmatrix} 
\begin{bmatrix} 
\vec{x_0} \ \vec{x_1} \ \cdots \ \vec{x}_{N-1} 
\end{bmatrix} \\
&= FX
\end{align*}$$



## 基于蝶形运算方法的 FFT 算法
