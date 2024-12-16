**<center><font face="黑体" size=6>真空球对称爱因斯坦方程的简单求解</font></center>**
**<center>周子杰</center>**

**<center>化学化工学院 化学强基2301</center>**

**摘要:** 爱因斯坦方程是广义相对论的基本方程。在广义相对论框架下，求解爱因斯坦方程后可以得到时空的度规张量，而度规张量决定了引力场下的物体该如何运动，因此爱因斯坦方程在广义相对论中有着重要的作用。真空球对称时空是最容易求解的一种情况，该情况下爱因斯坦方程的解被称为施瓦西解。施瓦西解可以解决许多牛顿力学无法完全解决的问题，如水星近日点的进动、光线弯折等物理问题。本文会给出施瓦西解的一个简单推导，并会简单讨论施瓦西解的一些性质与应用。

**关键词:** 真空球对称爱因斯坦方程；施瓦西解；水星近日点进动；光线弯折.

## 1.引言
&emsp;&emsp;广义相对论是描述宇宙四大基本相互作用之一——引力的一门学科。与狭义相对论不同，广义相对论框架下的背景时空并不平直，时空曲率不可忽略。广义相对论认为，引力的本质是时空的弯曲，在弯曲时空下，自由运动的物体不再做匀速直线运动，而是会沿着时空的测地线（令时空空中两点之间的线长取极值的轨迹）运动。在广义相对论框架下，时空的弯曲可以用度规张量来表示。不考虑宇宙学常数项，度规张量的运动方程可以表示为（为简化公式，本文采用抽象指标记号、几何单位制和爱因斯坦求和约定）

$$\begin{equation} R_{ab}-\frac{1}{2}g_{ab}R=8\pi T_{ab} \end{equation}$$

这就是爱因斯坦方程。其中$R_{ab}$是里奇张量，是以度规张量为自变量的函数；$g_{ab}$是度规张量，代表时空的弯曲；$R=g^{ab}R_{ab}$定义为标量曲率；$T_{ab}$为能动张量，代表时空中物质的分布情况。这是一个复杂的二阶非线性偏微分方程，通过求解该方程，我们可以得到给定物质分布的时空的度规信息，从而确定物体在该时空下的运动状况。

&emsp;&emsp;在一般情况下，方程（1）是十分难以求得解析解的，但是对于一些简单情况，我们可以求得方程（1）的解析解，本文要处理的真空球对称引力场就是为数不多的一类可以直接求得解析解的情况，这种情况下的爱英斯坦方程的解被称为施瓦西解。施瓦西解的应用非常广泛，比如可以用作微扰理论和格林函数的零阶近似，可以用来近似描述恒星周围的引力场等。本文的第二节将会给出真空球对称引力场的爱因斯坦方程的一个简单求解过程，并在第三节会对该解析解进行一些简单的讨论。

## 2.爱因斯坦方程的求解

&emsp;&emsp;由于爱因斯坦方程是一个复杂的二阶非线性偏微分方程，即使是真空球对称引力场这样的简单情况，爱因斯坦方程的求解都是非常困难且繁琐的。由于篇幅限制，本文在推导的过程中会省略一些步骤，读者对这些步骤感兴趣的话可以去阅读相关读物[^1]。

#### 2.1 前置知识

**1. 测地线方程**

&emsp;&emsp;给定度规张量$g_{ab}$，在该坐标系下的测地线方程可以写为

$$\begin{equation}\frac{d^2x^\sigma}{dt^2}+\Gamma^{\sigma}_{\mu\nu}\frac{dx^\mu}{dt}\frac{dx^\nu}{dt}=0\end{equation}$$

其中$x^\mu$是坐标的第$\mu$分量，$\Gamma_{\mu\nu}^{\sigma}$被称为克里斯托弗联络（简称克氏符）的分量形式，可以用度规张量表示为

$$\begin{equation}\Gamma^{\sigma}_{\mu\nu}=\frac{1}{2}g^{\sigma\rho}(g_{\rho\mu,\nu}+g_{\nu\rho,\mu}-g_{\mu\nu,\rho})\end{equation}$$

其中$g_{\mu\nu,\rho}=\frac{\partial g_{\mu\nu}}{\partial x^{\rho}}$代表度规张量的分量对坐标的导数。通过求解该测地线方程，可以解得给定度规张量的时空的测地线族，进而可以确定该度规下自由质点的运动轨迹。

**2. 真空爱因斯坦方程**

&emsp;&emsp;在真空条件下，物质的能动张量为零，方程（1）可简化为
$$\begin{equation}R_{ab}-\frac{1}{2}g_{ab}R=0\end{equation}$$

对方程（4）两边同时乘上$g^{ab}$，并对ab指标进行缩并，根据度规张量缩并的性质$g^{ab}g_{ab}=1$可得

$$\begin{equation}g^{ab}R_{ab}-\frac{1}{2}g^{ab}g_{ab}R=\frac{1}{2}R=0\end{equation}$$

得到真空条件下的标量曲率等于零，最后可以写出真空爱因斯坦方程为

$$\begin{equation}R_{ab}=0\end{equation}$$

所以在真空条件下，里奇张量恒为零。里奇张量是度规的函数，被定义为黎曼曲率张量的对第二、第四指标的缩并，其坐标分量可以用克氏符表示为

$$\begin{equation}R_{\mu\sigma}=R_{\mu\nu\sigma}^{\nu}=\Gamma^\nu_{\mu\sigma,\nu}-\Gamma^{nu}_{\nu\sigma,\mu}+\Gamma^{\lambda}_{\mu\sigma}\Gamma^{nu}_{\lambda\nu}-\Gamma^{\lambda}_{\nu\sigma}\Gamma^{nu}_{\lambda\mu}\end{equation}$$

而根据方程（3）克氏符也是与度规张量有关，所以方程（6）是对于度规张量的偏微分方程。

**3. 静态球对称度规**

&emsp;&emsp;设该时空下的线元为$ds^2$，可以得到线元在给定坐标系下的表达式

$$\begin{equation}ds^2=g_{\mu\nu}dx^{\mu}dx^{\nu}\end{equation}$$

在这里我们采用球坐标系，使用$t,r,\theta,\phi$作为坐标，即令$x^0=t,x^1=r,x^2=\theta,x^3=\phi$。

&emsp;&emsp;由于我们讨论的背景时空是静态的，所以度规不随时间的变化而变化，故度规张量各分量都不含时。又因为线元具有时间反演不变性，故$g_{0\mu}dtdx^{\mu}=-g_{0\mu}dtdx^\mu$，所以$g_{0\mu}=0$，故度规的时间坐标和空间坐标的交叉项为0。

&emsp;&emsp;由于我们讨论的背景时空是球对称的，所以在等 $t$ 等 $r$ 的情况下，度规就是以球面度规，即

$$\begin{equation}d\hat{s}^2=r^2(d\theta^2+sin^2\theta d\phi^2)\end{equation}$$

综上所述，我们可以得到静态球对称时空下线元的一般表达式

$$\begin{equation}\begin{aligned}ds^2&=g_{00}dt^2+g_{11}dr^2+d\hat{s}^2
\\\ &=g_{00}dt^2+g_{11}dr^2+r^2(d\theta^2+sin^2\theta d\phi^2)\end{aligned}\end{equation}$$

其中$g_{00}$和$g_{11}$是未知函数，可以通过爱因斯坦方程（1）确定。考虑到时空的球对称性，可以认为$g_{00}$和$g_{11}$不含角度变量。为了便于下一步的求解，可以令$g_{00}=-e^{2A(r)},g_{11}=e^{2B(r)}$，即

$$\begin{equation}ds^2=-e^{2A(r)}dt^2+e^{2B(r)}dr^2+r^2(d\theta^2+sin^2\theta d\phi^2)\end{equation}$$

**4. 线元（11）所对应的克氏符与里奇张量**

&emsp;&emsp;在之后的推导中，我们需要用到线元（11）所对应的克氏符与里奇张量，其计算可以直接将度规代入克氏符和里奇张量的计算式（3）（7）来进行计算。此处使用了Mathematica辅助公式推导，现将计算所得非零克氏符列举如下

$$\begin{equation}
\begin{aligned}
&\Gamma^{0}_{01}=\Gamma^{0}_{10}=A'\\
&\Gamma^{1}_{00}=A'e^{2(A-B)}\\
&\Gamma^{1}_{11}=B'\\
&\Gamma^{1}_{22}=-re^{-2B}\\
&\Gamma^{1}_{33}=-r\sin^2 \theta e^{-2B}\\
&\Gamma^{2}_{12}=\Gamma^{2}_{21}=\Gamma^{3}_{13}=\Gamma^{3}_{31}=\frac{1}{r}\\
&\Gamma^{2}_{33}=-\sin\theta \cos \theta \\
&\Gamma^{3}_{23}=\Gamma^{3}_{32}=\cot\theta
\end{aligned}
\end{equation}$$

将（12）代入（7）中，计算得到的非零里奇张量有
$$\begin{equation}
\begin{aligned}
&R_{00}=-e^{2(A-B)}(-A''+A'B'-A'^2-2r^{-1}A')\\
&R_{11}=A''+  A'^2-A'B'-\frac{2 B'}{r}\\
&R_{22}=e^{-2 B } \left(r\left(A'-B'\right)+1\right)-1\\
&R_{33}=   \left(-e^{-2 B}(r \left(B'-A'\right)-1)+1\right)\sin ^2\theta\\
\end{aligned}
\end{equation}$$
#### 2.2 施瓦西解的导出

&emsp;&emsp;在有了上面那些前置知识后，我们可以正式开始求解静态球对称时空下的真空爱因斯坦方程(6)了。写出方程（6）的坐标分量形式，即

$$\begin{equation}R_{\mu\nu}=0\end{equation}$$

将（13）代入方程（14），可得

$$\begin{equation}-A''+A'B'-A'^2-\frac{2A'}{r}=0\end{equation}$$

$$\begin{equation}A''+  A'^2-A'B'-\frac{2 B'}{r}=0\end{equation}$$

$$\begin{equation}e^{-2 B } \left(r\left(A'-B'\right)+1\right)-1=0\end{equation}$$

(14)加（15）可得

$$\begin{equation}A'+B'=0\end{equation}$$

故

$$\begin{equation}A=-B+C\end{equation}$$

其中$C$为任意常数。将（19）代入方程（17），化简可得

$$\begin{equation}1-2rB'=e^{2B}\end{equation}$$

将（20）移项，两边同时积分，可得

$$\begin{equation}e^{2B}=(1+\frac{D}{r})^{-1}\end{equation}$$

将（21）代入（19）可得

$$\begin{equation}e^{2A}=(1+\frac{D}{r})e^{2C}\end{equation}$$

由于广义相对性原理，进行坐标变换并不会改变物理规律，故总可选取新的时间坐标，使得$t=e^{2C}t$，将（21）（22）代入（11），得

$$\begin{equation}ds^2=-(1+\frac{D}{r})dt^2+(1+\frac{D}{r})^{-1}dr^2+r^2(d\theta^2+\sin^2\theta d\phi^2)\end{equation}$$

由于牛顿引力论是广义相对论的弱场近似，所以当引力场很弱（$r$很大）的时候，广义相对论表述就会退化成牛顿引力论就表述。根据广义相对论的线性近似[^2]，可以确定常数$D=-2M$，其中 $M$ 为引力源质量。综上所述，静态球对称时空下的真空爱因斯坦方程的解为

$$\begin{equation}ds^2=-(1-\frac{2M}{r})dt^2+(1-\frac{2M}{r})^{-1}dr^2+r^2(d\theta^2+\sin^2\theta d\phi^2)\end{equation}$$

这就是施瓦西解。

## 3.对施瓦西解的简单讨论

## 4.总结

#### 参考文献:

[^1]: 梁灿彬，周彬. 微分几何入门与广义相对论上[M]. 第二版. 北京：科学出版社，**2006**.
[^2]: 赵峥，刘文彪. 广义相对论基础[M]. 第一版. 北京：清华大学出版社，**2010**.
[^3]: 刘辽，费保俊，张允中. 狭义相对论[M]. 第二版. 北京：科学出版社，**2008**.