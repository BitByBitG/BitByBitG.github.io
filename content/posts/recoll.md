---
title: "P15649 [省选联考 2026] 找寻者：从概率 DP 到连续点值分治与实测优化"
date: 2026-10-19T09:41:00+08:00
draft: false
tags: ["随笔"]
---

本文给出两种完整做法：

| 做法 | 目标 | 时间复杂度 | 空间复杂度 |
| --- | --- | --- | --- |
| 连续点值、线性泛函与长链分治 | 保证渐进复杂度 | $O(n\log^2 n)$ | $O(n\log n)$ |
| 压缩分布、混合卷积与同配置区间批处理 | 优化原题 $n\le 5000$ 下的实际速度 | 保守上界 $O(n^2\log n)$ | $O(n\log n)$ |

第二种做法对应此前讨论中的“最优解版”。这里的含义是针对原题规模进行实测优化的版本，并不表示已经证明它在所有实现、所有机器、所有测试数据上都最快。

两种做法使用完全相同的概率 DP。差别在于，第一种方法用统一的代数结构批量处理长链，第二种方法根据分布宽度和侧儿子的具体情况选择计算方式。

文末附两份完整的 GNU C++17 程序。

## 一、题意与符号

[题目链接：P15649 找寻者](https://www.luogu.com.cn/problem/P15649)。原题有多组测试，每组 $n\le5000$，测试组数至多为 $5$。输入首行的测试点编号只用于题目说明，算法不依赖它。

给定一棵以 $1$ 为根的树。每个非叶结点都选择一个儿子作为重儿子，其余儿子对应的边为轻边。

选择过程自底向上进行。对于结点 $u$，假设各儿子向下延伸的重链长度分别为 $L_v$，则选择儿子 $v$ 的条件概率为

$$
\frac{L_v}{\sum_{w\in\operatorname{ch}(u)}L_w}.
$$

要求所有结点到根的路径上轻边数量的期望之和，对 $p=998244353$ 取模。

使用以下符号：

| 符号 | 含义 |
| --- | --- |
| $\operatorname{ch}(u)$ | $u$ 的儿子集合 |
| $\operatorname{sz}_u$ | $u$ 的子树大小 |
| $h_u$ | $u$ 子树的高度，按路径上的结点数计算 |
| $f_u[k]$ | 完成 $u$ 子树内的选择后，从 $u$ 向下延伸的重链长度为 $k$ 的概率 |
| $F_u(z)$ | 分布多项式 $\sum_k f_u[k]z^k$ |
| $M(m)$ | 长度为 $O(m)$ 的多项式乘法复杂度，使用 NTT 时为 $O(m\log(m+2))$ |

下文的概率等式先按有理数理解，再映射到模域。实际出现的链长之和小于模数，因此相应分母可逆。

## 二、将答案按边拆开

如果连接 $u$ 与其儿子 $v$ 的边是轻边，那么 $v$ 子树内的每一个结点都会受到一次贡献。

由期望的线性性，答案为

$$
\operatorname{Ans}
=\sum_u\sum_{v\in\operatorname{ch}(u)}
\operatorname{sz}_v\Pr(v\text{ is not selected by }u).
$$

因此，只要计算出每个儿子被选择的概率，就能得到答案。

需要注意：不同儿子子树内的随机选择互相独立，但一般有

$$
\mathbb E\left[\frac{L_v}{\sum_w L_w}\right]
\ne
\frac{\mathbb E[L_v]}{\sum_w\mathbb E[L_w]}.
$$

所以只维护链长的期望是不够的，我们维护完整的链长分布。

叶子的分布显然为

$$
F_u(z)=z.
$$

## 三、辅助长链剖分

### 3.1 固定一个长儿子

对每个非叶结点 $u$，固定高度最大的儿子 $a_u$ 作为辅助长儿子，其余儿子组成侧儿子集合 $S_u$。

这个选择只用于组织算法，和题目中的随机重儿子没有关系。

沿辅助长儿子形成若干条长链。在处理一条长链之前，先处理它的所有侧儿子所在的长链。于是，沿当前长链自底向上转移时，侧儿子的分布已经全部确定。

定义

$$
d_u=\sum_{v\in S_u}h_v.
$$

由于

$$
h_u=1+h_{a_u},
$$

对非叶结点有

$$
d_u=\sum_{v\in\operatorname{ch}(u)}h_v-h_u+1.
$$

叶子也满足上式，因为它没有儿子，且 $h_u=1$。

将所有结点相加，得到

$$
\begin{aligned}
\sum_u d_u
&=\sum_{v\ne 1}h_v-\sum_u h_u+n\\
&=n-h_1.
\end{aligned}
$$

因此

$$
\boxed{\sum_u d_u<n.}
$$

这是后续复杂度分析最重要的等式。

### 3.2 一条辅助长链的基本转移

设结点 $u$ 的长儿子为 $a$，长儿子的分布记为 $f[k]$。定义侧儿子分布的乘积

$$
A(z)=\prod_{v\in S_u}F_v(z)
=\sum_{s=0}^{d_u}A_s z^s.
$$

因为各侧儿子子树内的随机选择互相独立，$A_s$ 就是所有侧儿子的链长之和恰好为 $s$ 的概率。

如果存在侧儿子，每条侧链都非空，所以 $A_0=0$。

固定长儿子的链长为 $k$、侧儿子的链长总和为 $s$，选择长儿子的概率为 $k/(k+s)$。于是选择长儿子的输出为

$$
\boxed{
f^{\mathrm{keep}}[k+1]
=k f[k]\sum_{s=1}^{d_u}\frac{A_s}{k+s}.
}
$$

我们把继续选择辅助长儿子的转移称为“保留”，把选择任意侧儿子的转移称为“改选”。

### 3.3 侧儿子的输出和选择概率

定义查询

$$
\boxed{
r_s=\sum_k\frac{f[k]}{k+s},\qquad 1\le s\le d_u.
}
$$

对一个侧儿子 $v$，令

$$
B_v(z)=\prod_{\substack{w\in S_u\\w\ne v}}F_w(z).
$$

当 $v$ 的链长为 $\ell$、其他侧儿子的链长总和为 $j$ 时，选择 $v$ 的概率因子为

$$
\frac{\ell}{k+\ell+j}.
$$

将长儿子的链长 $k$ 求和，得到

$$
b_v[\ell]
=\sum_j[z^j]B_v(z)\,r_{\ell+j}.
$$

因此，选择 $v$ 且其链长为 $\ell$ 的概率为

$$
\ell f_v[\ell]b_v[\ell].
$$

所有改选产生的分布为

$$
\boxed{
J[\ell+1]
=\sum_{v\in S_u}\ell f_v[\ell]b_v[\ell].
}
$$

侧儿子 $v$ 被选择的概率为

$$
\pi_v=\sum_\ell\ell f_v[\ell]b_v[\ell].
$$

将保留与改选相加，得到 $u$ 的完整分布：

$$
f_u=f^{\mathrm{keep}}+J.
$$

### 3.4 当前结点的答案贡献

长儿子被选择的概率为 $1-\sum_{v\in S_u}\pi_v$。因此

$$
\begin{aligned}
C_u
&=\operatorname{sz}_a\sum_{v\in S_u}\pi_v
+\sum_{v\in S_u}\operatorname{sz}_v(1-\pi_v)\\
&=\sum_{v\in S_u}\operatorname{sz}_v
+\sum_{v\in S_u}
(\operatorname{sz}_a-\operatorname{sz}_v)\pi_v.
\end{aligned}
$$

即

$$
\boxed{
C_u=\sum_{v\in S_u}\operatorname{sz}_v
+\sum_{v\in S_u}(\operatorname{sz}_a-\operatorname{sz}_v)\pi_v.
}
$$

如果没有侧儿子，长儿子必选，答案贡献为零，分布只需整体右移一位。

### 3.5 用五点树检查公式

考虑边集 $(1,2),(1,3),(2,4),(2,5)$。

结点 $2$ 的两个儿子都是叶子，所以无论选谁，从 $2$ 向下的重链长度都是 $2$，即 $F_2(z)=z^2$。

在根结点处，选择长儿子 $2$ 的概率为 $2/3$，选择侧儿子 $3$ 的概率为 $1/3$。于是

$$
F_1(z)=\frac13z^2+\frac23z^3.
$$

根结点的答案贡献为

$$
C_1=1+(3-1)\cdot\frac13=\frac53.
$$

结点 $2$ 的两条儿子边中恰有一条轻边，贡献为 $1$，总答案为

$$
\operatorname{Ans}=\frac83\equiv 665496238\pmod{998244353}.
$$

## 四、侧儿子乘积树及其反传

如果对每个侧儿子都重新计算 $B_v$，会重复进行大量乘法。我们建立侧儿子分布的乘积树，根多项式为 $A$。

假设一个内部结点表示 $B(z)C(z)$，收到一组权重 $\lambda_j$，含义是要计算系数的加权和

$$
\sum_j\lambda_j[z^j](B(z)C(z)).
$$

固定 $B$ 的系数下标 $i$，它所对应的权重为

$$
\lambda_B[i]=\sum_j C_j\lambda_{i+j}.
$$

类似地，传给 $C$ 的权重为

$$
\lambda_C[j]=\sum_i B_i\lambda_{i+j}.
$$

两者都是相关卷积，将一个序列反转后即可使用普通卷积计算。

根处令 $\lambda_s=r_s$。反传到侧儿子 $v$ 时，得到的向量正好是 $b_v[\ell]$，从而同时计算 $J$、$\pi_v$ 和答案贡献。

每个结点只建立、反传一次侧儿子乘积树，费用为

$$
O(d_u\log^2(d_u+2)).
$$

全树总费用为

$$
O\left(\sum_u d_u\log^2(d_u+2)\right)
=O(n\log^2 n).
$$

所以，对于第一种做法，剩余问题是怎样在同样的复杂度内完成所有长链上的查询与分布传播。

## 五、做法一：$O(n\log^2 n)$ 的连续点值分治

### 5.1 为什么不能逐层维护完整分布

一条长链上，每经过一个结点，链长最多增加一。即使每层只有一个侧叶子，完整分布的宽度也可能逐渐增长到线性规模。

逐层扫描会产生平方级开销。

我们的处理方式是：先批量计算一段内“始终保留”的贡献，再通过分治处理区间内新发生的改选。

### 5.2 统一长链上的坐标

一条辅助长链最底部的叶子先初始化为分布 $f[1]=1$。将其上方的非叶结点自底向上编号为

$$
t=0,1,\ldots,q-1.
$$

在第 $t$ 层进入转移之前，若此前一直保留，就把链长写为 $X+t$。其中 $X$ 在这一段保留过程中保持不变。

令第 $t$ 层的侧儿子乘积为

$$
A_t(z)=\sum_{s=0}^{d_t}A_{t,s}z^s.
$$

该层保留的概率函数为

$$
\rho_t(X)
=(X+t)\sum_{s=1}^{d_t}\frac{A_{t,s}}{X+t+s}.
$$

定义

$$
q_t(X)=\prod_{s=1}^{d_t}(X+t+s),
$$

$$
p_t(X)
=(X+t)\sum_{s=1}^{d_t}
A_{t,s}\frac{q_t(X)}{X+t+s}.
$$

因为每个分母都是 $q_t$ 的因子，所以 $p_t$ 是多项式，并且

$$
\deg p_t\le d_t,\qquad \deg q_t=d_t.
$$

于是

$$
\rho_t(X)=\frac{p_t(X)}{q_t(X)}.
$$

没有侧儿子时，单独定义 $p_t=q_t=1$，保留概率为 $1$，但链长仍然增加一。

对于区间 $I=[l,r)$，记

$$
P_I=\prod_{t=l}^{r-1}p_t,\qquad
Q_I=\prod_{t=l}^{r-1}q_t,\qquad
D_I=\sum_{t=l}^{r-1}d_t.
$$

若区间入口的实际链长为 $k$，则对应的固定坐标为 $X=k-l$。整个区间都保留的概率为

$$
\boxed{\frac{P_I(k-l)}{Q_I(k-l)}.}
$$

### 5.3 连续点值表示

通常用系数表示多项式，但这里改为保存连续整数点上的值。

对于次数小于 $m$ 的多项式 $B$，定义

$$
\operatorname{val}_m(B)=(B(0),B(1),\ldots,B(m-1)).
$$

由于 $m<p$，这些点在模域中互不相同，$m$ 个值可以唯一确定 $B$。

定义映射 $E_{m,H,a}$：输入 $\operatorname{val}_m(B)$，输出

$$
(B(a),B(a+1),\ldots,B(a+H-1)).
$$

这里的 $a$ 可以是负整数。我们只使用有界整数区间，所有点数和阶乘参数都远小于模数。

下面证明，这个映射及其转置都可以通过卷积完成。

### 5.4 连续网格转移的卷积公式

设

$$
V_m(x)=\prod_{j=0}^{m-1}(x-j),
$$

$$
c_j=\frac{(-1)^{m-1-j}}{j!(m-1-j)!}.
$$

Lagrange 插值公式告诉我们，当 $x$ 不在原网格 $\{0,1,\ldots,m-1\}$ 中时，

$$
B(x)=V_m(x)\sum_{j=0}^{m-1}\frac{c_jB(j)}{x-j}.
$$

令 $x=a+i$，中间的和变为

$$
\sum_{j=0}^{m-1}\frac{c_jB(j)}{a+i-j}.
$$

反转序列 $c_jB(j)$ 后，这就是一个分母下标为两个序列下标之和的相关卷积。

另外，$V_m(x)$ 可以用阶乘、逆阶乘在常数时间求值：

$$
V_m(x)=
\begin{cases}
\dfrac{x!}{(x-m)!},&x\ge m,\\[6pt]
0,&0\le x<m,\\[6pt]
(-1)^m\dfrac{(m-1-x)!}{(-1-x)!},&x<0.
\end{cases}
$$

当目标点恰好落在原网格中时，直接复制原值，不使用上面的除法公式。

因此

$$
\boxed{E_{m,H,a}\text{ can be computed in }O(M(m+H)).}
$$

特别注意，复杂度只取决于输入、输出长度，不需要建立从 $0$ 到 $|a|$ 的额外长数组。

### 5.5 转置映射也只需要一次相关卷积

将 $E_{m,H,a}$ 看成一个 $H\times m$ 的矩阵。设输出位置的权重为 $w_i$，希望求出输入位置对应的权重。

这就是矩阵转置乘向量 $E_{m,H,a}^{\mathsf T}w$，不需要显式构造矩阵。

先将目标位置分成两类。对于不与原网格重合的点，定义

$$
t_i=w_iV_m(a+i).
$$

对于与原网格重合的点，令 $t_i=0$，随后单独补回直接复制产生的贡献。

第 $j$ 个输入位置得到的权重为

$$
\boxed{
(E_{m,H,a}^{\mathsf T}w)_j
=c_j\sum_{i=0}^{H-1}\frac{t_i}{a+i-j}
+\begin{cases}
w_{j-a},&0\le j-a<H,\\
0,&\text{otherwise}.
\end{cases}
}
$$

卷积内的零分母位置临时填零即可，因为对应的 $t_i$ 已经是零；真正的重合点贡献由后一项补回。这不表示数学上允许对零求逆。

第一项仍然是相关卷积，所以

$$
\boxed{E_{m,H,a}^{\mathsf T}\text{ can be computed in }O(M(m+H)).}
$$

代码中的 `Grid::eval`、`Grid::adjoint` 分别实现正向和转置转移；`Grid::cauchy` 负责分母形如 $b+i+j$ 的相关卷积。

### 5.6 在点值形式中构造区间乘积

区间 $I$ 保存

$$
\operatorname{val}_{D_I+1}(P_I),\qquad
\operatorname{val}_{D_I+1}(Q_I).
$$

必须用 $D_I+1$ 个点，因为 $P_I,Q_I$ 的次数可能等于 $D_I$。

如果 $I=L\cup R$，将左右儿子的点值都扩展到网格 $0,1,\ldots,D_I$，再对应位置相乘，即可得到父亲的点值。

一次合并的费用是 $O(M(D_I+1))$。

叶子也可以直接构造点值。对于 $0\le j\le d_t$，有

$$
q_t(j)=\frac{(j+t+d_t)!}{(j+t)!},
$$

$$
p_t(j)=(j+t)q_t(j)
\sum_{s=1}^{d_t}\frac{A_{t,s}}{j+t+s}.
$$

后一个求和又是相关卷积，因此叶子构造只需要 $O(M(d_t+1))$，不需要通用的多点求值。

### 5.7 把区间查询写成同一种加权求和

现在考虑入口分布 $f[1..H]$ 对区间 $I=[l,r)$ 的贡献，只统计进入区间后始终保留的历史。

固定入口链长 $k$，令

$$
x_k=k-l,\qquad w_k=\frac{f[k]}{Q_I(x_k)}.
$$

第 $t$ 层的查询需要把入口质量先传播到 $t$，再除以该层的链长与侧链长之和：

$$
r_{t,s}^{\mathrm{input}}
=\sum_{k=1}^H
\frac{f[k]}{x_k+t+s}
\prod_{a=l}^{t-1}\frac{p_a(x_k)}{q_a(x_k)}.
$$

统一乘上分母 $Q_I$，定义

$$
B_{t,s}(X)=
\frac{
\left(\prod_{a=l}^{t-1}p_a(X)\right)
\left(\prod_{a=t}^{r-1}q_a(X)\right)
}{X+t+s}.
$$

因为 $X+t+s$ 是 $q_t$ 的因子，所以 $B_{t,s}$ 是多项式，并且

$$
\deg B_{t,s}<D_I.
$$

于是所有查询都有相同形式：

$$
\boxed{
r_{t,s}^{\mathrm{input}}=\sum_k w_kB_{t,s}(x_k).
}
$$

不同查询之间，$x_k,w_k$ 都不变，只有多项式 $B_{t,s}$ 不同。

### 5.8 “线性泛函”在这里究竟是什么

定义一个操作

$$
\mathcal L_I(B)=\sum_k w_kB(x_k).
$$

它输入一个多项式，输出一个加权求和结果。因为

$$
\mathcal L_I(\alpha B+\gamma C)
=\alpha\mathcal L_I(B)+\gamma\mathcal L_I(C),
$$

这个操作被称为线性泛函。这里不需要额外的泛函分析知识，可以直接把它理解成一个固定的线性加权求和器。

例如，当 $B(X)=a+bX$ 时，

$$
\mathcal L_I(B)
=a\underbrace{\sum_k w_k}_{\mu_0}
+b\underbrace{\sum_k w_kx_k}_{\mu_1}.
$$

提前保存 $\mu_0,\mu_1$，之后换任何 $a,b$，都不用再次遍历所有 $k$。

一般地，次数小于 $D_I$ 的多项式构成一个 $D_I$ 维的线性空间，所以只需要 $D_I$ 个数就能描述这个求和器。

我们为了方便连续点值转移，保存权重向量 $\beta_I$，使得对任意次数小于 $D_I$ 的多项式 $B$，都有

$$
\boxed{
\mathcal L_I(B)=\sum_{j=0}^{D_I-1}\beta_I[j]B(j).
}
$$

也就是把原来的 $H$ 个加权求值点，换成 $D_I$ 个固定网格点上的权重，同时保证所有需要的多项式查询结果不变。

由转置的定义，直接得到

$$
\boxed{
\beta_I\mathrel{+}=E_{D_I,H,1-l}^{\mathsf T}w.
}
$$

这里用加法，是因为同一个区间可能接收多个来源的概率质量。由于操作是线性的，它们的权重可以直接累加。

### 5.9 `feed` 操作

定义 `feed(I,f)`，完成两件事：

1. 把入口分布按“区间内始终保留”传播到出口。
2. 将入口分布对区间内所有查询的贡献，暂存在区间根的 $\beta_I$ 中。

实现步骤如下：

1. 在连续点 $x_k=k-l$ 上求出 $P_I(x_k),Q_I(x_k)$。
2. 用批量求逆计算 $w_k=f[k]/Q_I(x_k)$。
3. 将 $E_{D_I,H,1-l}^{\mathsf T}w$ 加到 $\beta_I$。
4. 返回无改选出口

$$
\boxed{
f^{\mathrm{out}}[k+r-l]=w_kP_I(x_k).
}
$$

此时不递归地把 $\beta_I$ 传到所有叶子。

若 $D_I=0$，整个区间不存在侧儿子，直接平移分布，不建立泛函。

一次 `feed` 的费用为

$$
O\bigl(M(H+D_I+|I|+1)\bigr).
$$

### 5.10 泛函怎样向左右儿子传递

设 $I=L\cup R$。对于属于左区间的查询，其多项式在父区间中额外乘上 $Q_R$；对于属于右区间的查询，其多项式额外乘上 $P_L$。

因此，左右儿子对应的求和器满足

$$
\mathcal L_L(B)=\mathcal L_I(BQ_R),
$$

$$
\mathcal L_R(B)=\mathcal L_I(P_LB).
$$

注意次数仍然合法。例如 $\deg B<D_L$ 时，

$$
\deg(BQ_R)<D_L+D_R=D_I.
$$

用点值来实现这两次乘法，就得到

$$
\boxed{
\beta_L\mathrel{+}=
E_{D_L,D_I,0}^{\mathsf T}
\left(\beta_I\odot\operatorname{val}_{D_I}(Q_R)\right),
}
$$

$$
\boxed{
\beta_R\mathrel{+}=
E_{D_R,D_I,0}^{\mathsf T}
\left(\beta_I\odot\operatorname{val}_{D_I}(P_L)\right).
}
$$

其中 $\odot$ 表示逐项乘法。式中的点值长度为 $D_I$，这里只取对应网格上的值，不要求乘数本身的次数小于 $D_I$；乘数值由已保存的完整点值表示求出。

某个儿子的 $D$ 为零时，直接跳过。一次下传的费用为 $O(M(D_I+1))$。

### 5.11 到叶子后恢复查询

叶子 $t$ 的查询多项式为

$$
\frac{q_t(X)}{X+t+s}.
$$

所以

$$
\begin{aligned}
r_{t,s}
&=\mathcal L_t\left(\frac{q_t(X)}{X+t+s}\right)\\
&=\sum_{j=0}^{d_t-1}
\frac{\beta_t[j]q_t(j)}{j+t+s},
\qquad 1\le s\le d_t.
\end{aligned}
$$

这又是一次相关卷积，费用为 $O(M(d_t+1))$。

得到了完整的 $r_{t,s}$，就可以使用第四节的侧儿子乘积树，求出当前改选分布 $J_t$ 和答案贡献。

### 5.12 通过分治处理区间内新发生的改选

记 $U_I$ 为区间 $I$ 内始终保留的分布转移。

对于 $I=[l,r)$，定义

$$
G_I=\sum_{t=l}^{r-1}U_{[t+1,r)}J_t.
$$

它表示：最后一次改选发生在 $I$ 内，之后一直保留，传播到区间出口的概率质量。

如果 $I=L\cup R$，则

$$
\boxed{G_I=U_RG_L+G_R.}
$$

这给出如下分治流程：

```text
base = feed(root, delta_1)

solve(I):
    如果 D_I = 0：
        返回空分布

    如果 I 为叶子 t：
        由 beta_I 求出完整的 r_{t,s}
        计算当前结点答案贡献以及 J_t
        返回 J_t

    将 beta_I 向左右儿子各传一层

    G_L = solve(L)
    kept = feed(R, G_L)
    G_R = solve(R)

    返回 kept + G_R

整条长链的输出分布 = base + solve(root)
```

其中 `delta_1` 表示链底叶子的分布，即长度 $1$ 上的概率为 $1$。

这里有一个重要的时序：进入右区间前，父区间传下来的权重、左区间新改选质量产生的权重，都已经累加到右区间根。因此每个区间只需要向下传一次。

### 5.13 正确性证明

先考虑某个目标结点之前的随机历史。如果此前从未改选辅助长儿子，那么它属于链底初始分布的保留贡献，由最开始的 `feed(root,delta_1)` 负责。

否则，设此前最后一次改选发生在 $t$。这条历史包含在 $J_t$ 中，之后到目标结点之前全部是保留转移。

因此，以最后一次改选的位置分类，可以把完整分布唯一拆开，不会重计或漏计。

分治中，左区间产生的 $G_L$ 对右区间所有查询的贡献由一次 `feed(R,G_L)` 计算。右区间内部新发生的改选则由递归处理。每个“此前最后一次改选位置、目标查询位置”的配对，恰好在分治树上将二者分开的那次操作中处理。

所有部分来源对查询的贡献都是线性的。因此，可以先在区间根累加泛函，再一起向下传递。到达叶子时，获得的查询就是完整分布对应的查询，随后计算出的改选分布与答案贡献也都是正确的。

每个结点的答案常数项 $\sum_{v\in S_u}\operatorname{sz}_v$ 只在叶子最终结算时加入一次，不能对每个部分来源重复加入。

### 5.14 复杂度证明

对区间 $I$，定义

$$
W_I=|I|+D_I.
$$

若最后一次改选发生在 $t\in I$，所选侧儿子的链长为 $\ell$，那么传播到区间末端时，链长为

$$
\ell+r-t\le D_I+|I|=W_I.
$$

因此，`solve(I)` 返回的分布只需要 $O(W_I)$ 个位置，不会携带区间外输入的长分布。

在一个分治结点内，区间点值构造、泛函下传、左向右 `feed`、分布合并，总费用均可控制为

$$
O(M(W_I+1))=O(W_I\log(W_I+2)).
$$

按链上结点数二分，同一层区间互不相交，满足

$$
\sum_{I\text{ on one level}}W_I\le W,
$$

其中 $W=q+\sum_t d_t$。分治深度为 $O(\log(q+1))$，所以一条链的总费用是

$$
O(W\log^2(W+2)).
$$

全树有

$$
\sum_{\text{chains}}q\le n,
\qquad
\sum_{\text{chains}}\sum_t d_t=n-h_1.
$$

因此

$$
\sum_{\text{chains}}W=O(n).
$$

再加入侧儿子乘积树的费用，得到

$$
\boxed{O(n\log^2 n)\text{ time},\qquad O(n\log n)\text{ space}.}
$$

空间上，每层保存的点值总规模为 $O(n)$，乘以分治深度就是 $O(n\log n)$。中间分布同样按区间支撑上界保存，并及时释放临时数组。

### 5.15 实现时需要守住的边界

- $P_I,Q_I$ 需要 $D_I+1$ 个点值，泛函只需要 $D_I$ 个权重。二者不能混淆。
- 辅助长儿子按真实树高选择，不按模域中最高非零系数选择。模意义下最高系数可能偶然为零。
- 在实际入口点上，$Q_I$ 的每个因子都是 $k+(t-l)+s>0$，因此可以求逆；不需要对 $P_I$ 求逆。
- 插值网格重合时直接复制，转置时直接补回权重，不能真的对零求逆。
- 没有侧儿子时，只移动链长；改选分布为空，答案贡献为零。
- 分治返回值只包含区间内部改选产生的质量。若把原入口分布也塞进返回值，支撑和复杂度证明都会失效。
- `SideTree::apply` 只在该结点的查询完整后调用一次，避免重复累计答案常数项。

## 六、做法二：针对原题规模的实测优化版

第一种做法对所有长链使用统一的多项式框架，渐进复杂度很好，但原题的 $n$ 只有 $5000$。对很窄的分布使用多次 NTT、插值、转置，实际常数并不划算。

实测优化版保留共同 DP，但改为以下策略：

1. 分布只保存有效区间，确定的链长只占一个系数。
2. 只有一个儿子时，只修改下标偏移。
3. 所有侧儿子的链长都确定时，直接使用更简单的公式。
4. 一般转移中，小规模直接计算，大规模改用 NTT。
5. 连续多层的确定侧链配置相同时，使用生成函数一次处理整段。
6. 在直接计算内核中减少取模，并可选使用 AVX2。

这些分支是对同一精确 DP 的不同计算方式，不包含随机近似。

### 6.1 压缩分布的存储区间

把分布写成

$$
F(z)=z^L\overline F(z),
$$

其中 $L$ 是保留区间的左端点，数组只保存 $\overline F$ 的系数，即

$$
\overline F[i]=f[L+i].
$$

例如，确定长度为 $1000$ 的分布 $F(z)=z^{1000}$，只保存 `lo=1000` 和数组 `{1}`，不保存前面一千个零。

实现中的区间允许保留少量零系数，不要求每次都裁到模域意义下最紧的支撑。只要真实非零系数全部包含在区间中，转移就正确。

两个分布相乘时

$$
F(z)G(z)=z^{L_F+L_G}\overline F(z)\overline G(z).
$$

因此偏移直接相加，系数数组照常卷积。

对于只有一个儿子的结点，重儿子必选，只需把 `lo` 加一，数组本身无需移动。

### 6.2 所有侧儿子的链长都确定

设侧儿子 $v$ 的确定链长为 $\ell_v$，定义

$$
S=\sum_{v\in S_u}\ell_v.
$$

此时侧儿子乘积为单项式

$$
A(z)=z^S.
$$

只需要一个查询

$$
r=\sum_k\frac{f[k]}{k+S}.
$$

侧儿子 $v$ 被选择的概率为

$$
\pi_v=\ell_v r,
$$

改选输出为

$$
J[\ell+1]
=\left(\sum_{v:\ell_v=\ell}\ell_v\right)r.
$$

保留输出则可以改写为

$$
\begin{aligned}
f^{\mathrm{keep}}[k+1]
&=\frac{k}{k+S}f[k]\\
&=f[k]-S\frac{f[k]}{k+S}.
\end{aligned}
$$

因此，遍历长儿子的分布时，先计算

$$
z_k=\frac{f[k]}{k+S},
$$

一方面把 $z_k$ 加到 $r$，另一方面得到 $f[k]-Sz_k$。若 $S=1$ 或 $S=2$，后一项甚至可以只用加减法计算。

当前结点的答案为

$$
\boxed{
C_u=\sum_{v\in S_u}\operatorname{sz}_v
+\left(
\operatorname{sz}_a S
-\sum_{v\in S_u}\operatorname{sz}_v\ell_v
\right)r.
}
$$

整个分支不需要建立侧儿子乘积树。这里要求的是链长确定，不要求侧子树本身是一条链；例如所有叶子深度相同的某些分叉树，其根链长也可能确定。

### 6.3 一般侧儿子：直接计算与 NTT 混合

设压缩后的长儿子分布为

$$
f[L+i]=a_i,\qquad 0\le i<H,
$$

压缩后的侧儿子乘积为

$$
A(z)=z^B\sum_{j=0}^{D-1}c_jz^j.
$$

所需查询为

$$
r_j=\sum_{i=0}^{H-1}\frac{a_i}{L+B+i+j},
\qquad 0\le j<D.
$$

保留输出为

$$
f^{\mathrm{keep}}[L+i+1]
=(L+i)a_i
\sum_{j=0}^{D-1}\frac{c_j}{L+B+i+j}.
$$

如果 $HD$ 较小，直接双重循环。可以共用中间量

$$
z_{i,j}=\frac{a_i}{L+B+i+j}.
$$

它同时贡献给查询 $r_j$ 和当前链长的保留输出。

如果规模较大，定义核序列

$$
K_t=\frac1{L+B+t},\qquad 0\le t<H+D-1.
$$

记反转后的序列为 $a^{\mathrm{rev}},c^{\mathrm{rev}}$，则

$$
\boxed{r_j=(a^{\mathrm{rev}}*K)_{H-1+j},}
$$

$$
\boxed{
\sum_{j=0}^{D-1}\frac{c_j}{L+B+i+j}
=(c^{\mathrm{rev}}*K)_{D-1+i}.
}
$$

两次 NTT 卷积即可得到全部结果。

随后，依然使用侧儿子乘积树反传，只是所有数组都带有偏移。根的 $r_j$ 对应实际侧链总长 $B+j$；向下反传时，偏移由兄弟区间的最低次数自动相加，系数数组的相关卷积公式不变。

程序使用的直接计算条件是 $HD\le 80000$ 或 $\min(H,D)\le4$。这些阈值只影响速度，不影响正确性，也不表示在所有机器上都最优。

### 6.4 连续相同的确定侧链配置

考虑长链上连续的 $q$ 层，满足下面两个条件：

1. 每一层的所有侧儿子，其链长都确定。
2. 每一层的侧链长度多重集合完全相同。

设每一层恰有 $c_\ell$ 个长度为 $\ell$ 的侧链，定义

$$
S=\sum_{\ell\ge1}c_\ell\ell,
\qquad
m_{\max}=\max\{\ell:c_\ell>0\}.
$$

注意，不要求侧子树同构，也不要求侧子树大小相同。只要链长配置相同，分布转移就相同；不同的子树大小只影响最后的答案系数。

把这段内部的层编号为 $t=0,1,\ldots,q-1$。设入口分布为 $f[k]$，第 $t$ 层入口的完整分布为 $f_t[k]$。

定义

$$
r_t=\sum_k\frac{f_t[k]}{k+S}.
$$

这里 $r_t$ 本身不是总改选概率，总改选概率为 $Sr_t$。长度为 $\ell$ 的侧链共向输出位置 $\ell+1$ 加入

$$
c_\ell\ell r_t.
$$

由于每层的分布转移相同，按此前最后一次改选的位置分类后，层间的影响只取决于时间差。这就能够建立固定卷积递推。

### 6.5 连续保留的概率

从长度 $k$ 开始，连续保留 $t$ 次的概率为

$$
\begin{aligned}
\prod_{a=0}^{t-1}\frac{k+a}{k+a+S}
&=\frac{(k+t-1)!}{(k-1)!}
\frac{(k+S-1)!}{(k+t+S-1)!}.
\end{aligned}
$$

当 $t=0$ 时，空乘积为 $1$，右边也等于 $1$。

如果此后要计算第 $t$ 层的查询，还需要再乘上

$$
\frac1{k+t+S}.
$$

为了统一这些阶乘比，定义

$$
T_j=\frac{(j-1)!}{(j+S)!},\qquad j\ge1,
$$

以及入口分布的缩放系数

$$
b_k=f[k]\frac{(k+S-1)!}{(k-1)!}.
$$

### 6.6 原入口分布对各层查询的贡献

对于入口长度 $k$，在此前始终保留的情况下，它对第 $t$ 层查询的贡献为

$$
\begin{aligned}
&f[k]
\frac{(k+t-1)!}{(k-1)!}
\frac{(k+S-1)!}{(k+t+S-1)!}
\frac1{k+t+S}\\
&\qquad=b_kT_{k+t}.
\end{aligned}
$$

于是定义

$$
\boxed{a_t=\sum_k b_kT_{k+t}.}
$$

全部 $a_t$ 是一个相关卷积，可以批量求出。

若入口压缩区间是 $[L_0,H_0]$，数组保存 $b_{L_0+i}$，那么将这个数组反转后与 $T$ 卷积，取下标 $H_0+t$，得到的正好是 $a_t$。这个下标也与代码中的 `c[H+t]` 对应。

### 6.7 一次改选对后续查询的影响

假设在第 $i$ 层改选了长度为 $\ell$ 的侧链，输出长度为 $\ell+1$，产生的质量为

$$
c_\ell\ell r_i.
$$

考虑第 $t=i+j$ 层，且 $j\ge1$。此前要保留 $j-1$ 次。到达目标查询时，链长为 $\ell+j$。

因此，这部分质量对目标查询的贡献为

$$
\begin{aligned}
&c_\ell\ell r_i
\left(\prod_{a=1}^{j-1}\frac{\ell+a}{\ell+a+S}\right)
\frac1{\ell+j+S}\\
&\qquad=
r_i\left(c_\ell\ell\frac{(\ell+S)!}{\ell!}\right)
\frac{(\ell+j-1)!}{(\ell+j+S)!}.
\end{aligned}
$$

定义

$$
w_\ell=c_\ell\ell\frac{(\ell+S)!}{\ell!}.
$$

得到时间差为 $j$ 时的卷积核

$$
\boxed{
K_j=\sum_{\ell\ge1}w_\ell T_{\ell+j},\qquad j\ge1.
}
$$

全部 $K_j$ 同样只需要一次相关卷积。

必须单独规定

$$
K_0=0.
$$

因为当前层发生的改选只能影响后面的层，不能再影响当前层已经用于决定这次改选的查询。不能把上面的求和式不加区分地用于 $j=0$。

### 6.8 一次多项式求逆求出整段查询

由最后一次改选位置的分类，有

$$
\boxed{
r_t=a_t+\sum_{i=0}^{t-1}r_iK_{t-i}.
}
$$

定义截断生成函数

$$
\mathcal A(z)=\sum_{t=0}^{q-1}a_tz^t,
\qquad
\mathcal K(z)=\sum_{j=1}^{q-1}K_jz^j,
\qquad
\mathcal R(z)=\sum_{t=0}^{q-1}r_tz^t.
$$

递推等价于

$$
\mathcal R(z)
\equiv\mathcal A(z)+\mathcal K(z)\mathcal R(z)
\pmod{z^q}.
$$

因此

$$
\boxed{
\mathcal R(z)
\equiv\frac{\mathcal A(z)}{1-\mathcal K(z)}
\pmod{z^q}.
}
$$

分母的常数项为 $1$，可以使用 Newton 迭代求逆。

若 $B$ 已经是 $1-\mathcal K$ 在低阶上的逆，那么将精度翻倍的迭代为

$$
B_{\mathrm{new}}
\equiv B\bigl(2-(1-\mathcal K)B\bigr)
\pmod{z^{2m}}.
$$

全部迭代的费用为 $O(M(q))$。

于是，这段内所有查询可以通过两次相关卷积、一次多项式求逆和一次乘法得到。

### 6.9 同时恢复完整出口分布

只求出 $r_t$ 还不够，因为上面的长链仍然需要当前区间的完整输出分布。

定义

$$
v_L=\frac{(L-1)!}{(L+S-1)!}.
$$

首先考虑原入口中从未改选的质量。若入口长度为 $k$，出口长度就是 $L=k+q$，质量为

$$
\begin{aligned}
f[k]\prod_{a=0}^{q-1}\frac{k+a}{k+a+S}
&=b_k\frac{(k+q-1)!}{(k+q+S-1)!}\\
&=b_kv_L.
\end{aligned}
$$

接着考虑最后一次改选发生在第 $i$ 层，选择长度为 $\ell$ 的侧链。出口长度为

$$
L=\ell+q-i.
$$

这部分输出质量为

$$
\begin{aligned}
&c_\ell\ell r_i
\prod_{a=1}^{q-i-1}\frac{\ell+a}{\ell+a+S}\\
&\qquad=r_iw_\ell v_L.
\end{aligned}
$$

定义

$$
W(z)=\sum_{\ell\ge1}w_\ell z^\ell,
$$

$$
\mathcal R^{\mathrm{rev}}(z)
=\sum_{i=0}^{q-1}r_i z^{q-1-i}.
$$

在乘积 $W(z)\mathcal R^{\mathrm{rev}}(z)$ 中，下标 $L-1$ 对应的正是所有满足 $L=\ell+q-i$ 的配对。

因此整个出口分布可以统一写为

$$
\boxed{
f^{\mathrm{out}}[L]
=v_L\left(
b_{L-q}
+[z^{L-1}]W(z)\mathcal R^{\mathrm{rev}}(z)
\right).
}
$$

其中 $b_k$ 在入口分布范围外定义为零。

这样只需再做一次卷积，就能恢复所有出口系数，不会因为要输出完整分布重新退化成逐层扫描。

辅助长儿子按高度选择，所以本段入口的高度上界 $H_0$ 不小于最大侧链长度 $m_{\max}$。出口长度不超过 $H_0+q$，批处理数组也按这个上界分配。

### 6.10 整段的答案贡献

侧链长度配置相同，但子树大小可以因层而异。设第 $t$ 层的辅助长儿子为 $a_t$，当前侧儿子集合为 $S_t$。

这一层的答案贡献为

$$
\boxed{
C_t=
\sum_{v\in S_t}\operatorname{sz}_v
+\left(
\operatorname{sz}_{a_t}S
-\sum_{v\in S_t}\operatorname{sz}_v\ell_v
\right)r_t.
}
$$

已经求出了所有 $r_t$，所以沿区间线性扫描一次，就能把全部答案贡献加入。

### 6.11 一个特例：每层挂一个叶子的梳子树

当每层只有一个长度为 $1$ 的侧叶子，且入口分布是长度 $1$ 时，有

$$
S=1,\qquad
a_t=\frac1{(t+1)(t+2)},
$$

$$
K_j=\frac2{(j+1)(j+2)},\qquad j\ge1.
$$

因此

$$
r_t=
\frac1{(t+1)(t+2)}
+\sum_{i=0}^{t-1}
\frac{2r_i}{(t-i+1)(t-i+2)}.
$$

这是前面一般批处理公式的直接特例，也解释了为什么长梳子树可以从逐层平方扫描变成一次近线性的多项式计算。

### 6.12 批处理的启用条件

不能认为满足相同配置就一定值得批处理。很短的区间，或者原本分布就很窄的区间，直接循环通常更快。

代码先将每层的确定侧链长度排序，作为该层的 `key`。如果存在非确定侧链，或者没有侧儿子，就不将该层视为可批处理的非空配置。

对连续相同的非空 `key`，设区间长度为 $q$、入口分布的最高保留下标为 $H_0$。当前实现仅在

$$
q\ge128
$$

且

$$
qH_0+\frac{q^2}{2}\ge200000
$$

时启用批处理。代码中使用整数运算表达该估算式。

它估计的是逐层扫描的工作量，而不是正确性条件。阈值取不同数值，算法仍然精确，只是运行时间不同。

若批处理区间的入口最高下标为 $H_0$、最大侧链长度为 $m_{\max}$，整段费用为

$$
O(M(H_0+q+m_{\max})+q+\text{number of side children}).
$$

这里使用长度上界 $H_0$ 而不是入口宽度，是因为实现中用于阶乘核和出口的数组按实际长度坐标建立。由于只在足够大的区间启用，这个开销在原题规模下通常划算。

### 6.13 减少取模次数

模乘是直接计算中的主要开销之一。当前模数满足

$$
p<10^9.
$$

两个标准模域元素的乘积小于 $10^{18}$，而

$$
16(p-1)^2<2^{64}.
$$

所以，朴素卷积和相关卷积中，可以用无符号 $64$ 位整数累计最多 $16$ 个乘积，然后统一取模。

直接 DP 中，每个 $z_{i,j}$ 已经是模域元素。对固定的 $j$，最多累加 $5000$ 个这样的元素，使用 $64$ 位整数同样安全，可以在这一列结束后统一取模。

这两种累加不能混淆：未取模的乘积只能按安全的小块累计，已经取模的元素则可以累计更多次。

### 6.14 AVX2 模乘内核

对于确定侧链分支，程序可以一次处理四个 $z_k=f[k]/(k+S)$。

内核用双精度浮点数估计整数商

$$
\widehat q\approx\left\lfloor\frac{ab}{p}\right\rfloor,
$$

但余数由精确的 $64$ 位整数运算计算：

$$
r=ab-\widehat q p.
$$

由于 $a,b<p<2^{30}$，双精度计算的绝对商误差远小于 $1$。估计商与正确的向下取整商至多相差 $1$，所以

$$
-p<r<2p.
$$

随后，如果 $r<0$ 就加 $p$，如果 $r\ge p$ 就减 $p$，得到精确模乘结果。

因此，这个分支并不是用浮点数近似概率或近似答案；浮点数只用于估计整数商，最终模域结果仍然精确。

程序运行时检测 AVX2 支持，不支持时使用普通循环，编译时不必额外添加 `-mavx2`。该内核不使用 `popcount`，因此单独开启 `avx512vpopcntdq` 不会改善这部分算法。

### 6.15 实测优化版的正确性

所有分支都等价于第三节的基本 DP：

- 单儿子分支利用该儿子必选的事实，只移动长度。
- 确定侧链分支将侧儿子乘积简化为单项式，直接计算同一个查询和同一个答案公式。
- 一般分支直接或通过 NTT 计算完全相同的相关卷积。
- 批处理分支按此前最后一次改选的位置分类，将同一转移写成固定卷积递推，并同时恢复完整出口分布。
- 延迟取模和 SIMD 只改变算术实现，不改变模域运算结果。

按辅助长链的依赖关系自底向上处理，所有侧儿子分布都已正确得到；沿每条长链，逐步或按区间进行等价转移。因此最终根分布以及所有结点的答案贡献都正确。

### 6.16 一个保守的复杂度上界

这一版本的目标是实测速度，不追求最佳的统一渐进上界。可以给出如下保守分析。

设某个分叉结点处，长儿子的压缩宽度为 $H_u$，侧儿子乘积的压缩宽度为 $D_u$。有

$$
H_u\le n,\qquad D_u\le d_u+1.
$$

于是

$$
\sum_{u\text{ branching}}H_uD_u
\le n\sum_{u\text{ branching}}(d_u+1)
=O(n^2).
$$

直接计算不会超过这部分上界。使用 NTT 时，单次费用为

$$
O((H_u+D_u)\log(n+2))
=O(H_uD_u\log(n+2)),
$$

因为 $H_u,D_u\ge1$。

所有侧儿子乘积树及反传仍然共需 $O(n\log^2 n)$。不同批处理区间互不重叠，区间数量至多 $n$，每个区间处理的数组长度为 $O(n)$，因此其总费用也有 $O(n^2\log n)$ 的保守上界。

所以当前实现可以保证

$$
\boxed{O(n^2\log n)\text{ time},\qquad O(n\log n)\text{ space}.}
$$

这个时间上界很松，不代表原题中真的会执行这么多运算。

空间方面，辅助长链上只保留当前分布，已消耗的长儿子数组及时释放；尚未合并的链头分布对应互不重叠的辅助长链，其总存储长度为 $O(n)$。侧儿子乘积树和反传临时数组占 $O(n\log n)$，批处理和 NTT 临时数组占 $O(n)$。

## 七、两种做法怎样对应到代码

### 7.1 渐进复杂度版

| 代码部分 | 数学含义 |
| --- | --- |
| `Grid::cauchy` | 分母为 $b+i+j$ 的相关卷积 |
| `Grid::eval` | 连续点值映射 $E_{m,H,a}$ |
| `Grid::adjoint` | 转置映射 $E_{m,H,a}^{\mathsf T}$ |
| `SideTree::build` | 侧儿子分布乘积树 |
| `SideTree::push` | 系数泛函反传，计算侧儿子选择概率及改选分布 |
| `Chain::build` | 用连续点值构造区间 $P_I,Q_I$ |
| `Chain::feed` | 无改选出口及区间根泛函累加 |
| `Chain::push` | 区间泛函向左右儿子传一层 |
| `Chain::solve` | 按最后一次改选来源进行分治 |
| `Chain::run` | 合并链底初始来源和区间改选来源 |

### 7.2 实测优化版

| 代码部分 | 数学含义 |
| --- | --- |
| `Dist` | 偏移加系数数组的压缩分布 |
| `det_update` | 确定侧链时，同时求查询与保留输出 |
| `Solver::build`、`push` | 压缩后的侧儿子乘积树及反传 |
| `Solver::step` | 单结点的直接计算或 NTT 转移 |
| `inverse_series` | Newton 迭代求幂级数逆 |
| `Solver::batch` | 同配置区间的固定卷积递推和完整输出恢复 |
| `keys` | 每层确定侧链长度的有序多重集合 |

两份程序都先以 $1$ 为根计算子树大小和高度，再按辅助长链的依赖关系自底向上处理，不需要深度为 $n$ 的递归。

## 八、实现检查与实际表现

### 8.1 边界情况

- $n=1$：没有边，答案为 $0$，分布为 $z$。
- 整棵树是一条链：每个非叶结点只有一个儿子，所有边必为重边，答案为 $0$。
- 根连着 $n-1$ 个叶子的星形树：根恰好选一个叶子，其他 $n-2$ 条边是轻边，答案为 $n-2$。
- 非确定侧链：不能使用单项式或固定配置分支，必须回到一般 DP。
- 相同长度配置但不同子树大小：允许批处理分布，但答案系数仍要逐层按真实子树大小计算。
- 模域中某个系数为零：不影响多项式运算，但辅助树高不能由最高非零模系数代替。
- 多组测试：图、分布和答案按测试组重建；优化版 NTT 单位根表可以跨组复用。

### 8.2 对拍与运行时间

在此前的实现验证中，使用了独立朴素树 DP，比较的不只是最终答案，还包括根结点的完整链长分布。

验证内容包括：2142 棵小树、173 棵结构或混合大树、批处理强制启用测试、SIMD 与普通内核的随机及边界值对比，以及 AddressSanitizer、UBSan、STL 边界检查。

此前同一运行环境下，使用 GNU C++17、`-O2`，对 23 类结构、共 115 棵 $n=5000$ 的构造树重复计时，记录总耗时中位数：

| 实现 | 总耗时中位数 |
| --- | --- |
| 本文做法一的 $O(n\log^2 n)$ 实现 | 约 $7.335$ 秒 |
| 讨论过程中的旧分块实现 | 约 $1.362$ 秒 |
| 本文做法二的实测优化实现 | 约 $0.298$ 秒 |

这些数据来自本地构造测试，并不是洛谷隐藏数据的运行时间，也不是对排行榜名次的保证。具体时间会随处理器、编译器和测试结构变化。

### 8.3 编译方式

两份代码都可使用以下命令编译：

```bash
g++ -std=c++17 -O2 -Wall -Wextra P15649.cpp -o P15649
g++ -std=c++17 -O2 -Wall -Wextra P15649_fast.cpp -o P15649_fast
```

优化版中的 AVX2 函数使用独立的目标属性，并进行运行时检测；正常情况下无需把整个程序都编译成只支持 AVX2 的程序。

## 附录 A：$O(n\log^2 n)$ 完整代码


```cpp
#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using poly=vector<int>;
constexpr int mod=998244353;
int add(int x,int y){return x+y>=mod?x+y-mod:x+y;}
int sub(int x,int y){return x<y?x-y+mod:x-y;}
int mul(int x,int y){return (ll)x*y%mod;}
int power(int x,int y)
{
	int z=1;
	for(;y;y>>=1,x=mul(x,x))if(y&1)z=mul(z,x);
	return z;
}
void ntt(poly& a,bool inv)
{
	int n=a.size();
	for(int i=1,j=0;i<n;i++)
	{
		int k=n>>1;
		for(;j&k;k>>=1)j^=k;
		j^=k;
		if(i<j)swap(a[i],a[j]);
	}
	for(int len=2;len<=n;len<<=1)
	{
		int w=power(3,(mod-1)/len);
		if(inv)w=power(w,mod-2);
		for(int l=0;l<n;l+=len)
			for(int j=0,z=1;j<len/2;j++,z=mul(z,w))
			{
				int x=a[l+j],y=mul(a[l+j+len/2],z);
				a[l+j]=add(x,y);a[l+j+len/2]=sub(x,y);
			}
	}
	if(inv)
	{
		int z=power(n,mod-2);
		for(int& x:a)x=mul(x,z);
	}
}
poly multiply(const poly& a,const poly& b)
{
	if(a.empty()||b.empty())return {};
	int n=a.size(),m=b.size();
	if(min(n,m)<=8||(ll)n*m<=1024)
	{
		poly c(n+m-1);
		for(int i=0;i<n;i++)for(int j=0;j<m;j++)c[i+j]=add(c[i+j],mul(a[i],b[j]));
		return c;
	}
	int k=1;
	while(k<n+m-1)k<<=1;
	poly x=a,y=b;x.resize(k);y.resize(k);
	ntt(x,0);ntt(y,0);
	for(int i=0;i<k;i++)x[i]=mul(x[i],y[i]);
	ntt(x,1);x.resize(n+m-1);
	return x;
}
void trim(poly& a){while(!a.empty()&&!a.back())a.pop_back();}
void append(poly& a,const poly& b)
{
	if(a.size()<b.size())a.resize(b.size());
	for(int i=0;i<(int)b.size();i++)a[i]=add(a[i],b[i]);
}
poly transpose(const poly& a,const poly& b,int n)
{
	poly c=multiply(a,poly(b.rbegin(),b.rend())),r(n);
	for(int i=0;i<n;i++)r[i]=c[i+b.size()-1];
	return r;
}
struct Grid
{
	poly fac,ifac,inv;
	Grid(int n):fac(n+1,1),ifac(n+1,1),inv(n+1)
	{
		inv[1]=1;
		for(int i=2;i<=n;i++)inv[i]=mod-mul(mod/i,inv[mod%i]);
		for(int i=1;i<=n;i++)fac[i]=mul(fac[i-1],i);
		ifac[n]=power(fac[n],mod-2);
		for(int i=n;i;i--)ifac[i-1]=mul(ifac[i],i);
	}
	int inverse(int x)const{return x<0?mod-inv[-x]:inv[x];}
	int falling(int m,int x)const
	{
		if(x>=0&&x<m)return 0;
		if(x>=m)return mul(fac[x],ifac[x-m]);
		int z=mul(fac[m-1-x],ifac[-1-x]);
		return m&1?mod-z:z;
	}
	int coef(int m,int j)const
	{
		int z=mul(ifac[j],ifac[m-1-j]);
		return (m-1-j)&1?mod-z:z;
	}
	poly cauchy(const poly& a,int n,int b)const
	{
		int m=a.size();poly r(n);
		if(!m||!n)return r;
		if((ll)m*n<=1024)
		{
			for(int i=0;i<n;i++)for(int j=0;j<m;j++)
				r[i]=add(r[i],mul(a[j],inverse(b+i+j)));
			return r;
		}
		poly k(n+m-1);
		for(int i=0;i<n+m-1;i++)k[i]=inverse(b+i);
		poly c=multiply(poly(a.rbegin(),a.rend()),k);
		for(int i=0;i<n;i++)r[i]=c[m-1+i];
		return r;
	}
	poly eval(const poly& v,int n,int a=0)const
	{
		int m=v.size();
		if(!n)return {};
		if(m==1)return poly(n,v[0]);
		if(a>=0&&a+n<=m)return poly(v.begin()+a,v.begin()+a+n);
		poly c(m);
		for(int j=0;j<m;j++)c[m-1-j]=mul(v[j],coef(m,j));
		poly r=cauchy(c,n,a-m+1);
		for(int i=0;i<n;i++)
		{
			int x=a+i;
			r[i]=x>=0&&x<m?v[x]:mul(r[i],falling(m,x));
		}
		return r;
	}
	poly adjoint(const poly& w,int m,int a=0)const
	{
		int n=w.size();poly r(m);
		if(!m||!n)return r;
		if(m==1)
		{
			for(int x:w)r[0]=add(r[0],x);
			return r;
		}
		if(a>=0&&a+n<=m)
		{
			copy(w.begin(),w.end(),r.begin()+a);return r;
		}
		poly v(n);
		for(int i=0;i<n;i++)v[i]=mul(w[i],falling(m,a+i));
		poly c=cauchy(v,m,a-m+1);
		for(int j=0;j<m;j++)
		{
			r[j]=mul(c[m-1-j],coef(m,j));
			if(j>=a&&j<a+n)r[j]=add(r[j],w[j-a]);
		}
		return r;
	}
};
struct SideTree
{
	struct Node{poly a;int l=-1,r=-1,sz=0;};
	vector<Node> tr;
	int root=-1,sum=0,h=0;
	int build(const vector<int>& v,int l,int r,vector<poly>& f,const vector<int>& sz)
	{
		int u=tr.size();tr.emplace_back();
		if(r-l==1)
		{
			int x=v[l];tr[u].a=move(f[x]);tr[u].sz=sz[x];
			sum+=sz[x];h=max(h,(int)tr[u].a.size()-1);
		}
		else
		{
			int m=(l+r)/2,x=build(v,l,m,f,sz),y=build(v,m,r,f,sz);
			tr[u].l=x;tr[u].r=y;tr[u].a=multiply(tr[x].a,tr[y].a);
		}
		return u;
	}
	void init(const vector<int>& v,vector<poly>& f,const vector<int>& sz)
	{
		if(v.empty())return;
		tr.reserve(2*v.size());root=build(v,0,v.size(),f,sz);
	}
	void push(int u,const poly& b,poly& out,int hs,int& ans)const
	{
		const auto& z=tr[u];
		if(z.l<0)
		{
			int p=0;
			for(int k=1;k<(int)z.a.size();k++)
			{
				int x=mul(mul(k,z.a[k]),b[k]);
				out[k+1]=add(out[k+1],x);p=add(p,x);
			}
			ans=add(ans,mul((hs-z.sz+mod)%mod,p));return;
		}
		poly x=transpose(b,tr[z.r].a,tr[z.l].a.size());
		poly y=transpose(b,tr[z.l].a,tr[z.r].a.size());
		push(z.l,x,out,hs,ans);push(z.r,y,out,hs,ans);
	}
	poly apply(const poly& b,int hs,int& ans)const
	{
		if(root<0)return {};
		poly out(h+2);ans=add(ans,sum);
		push(root,b,out,hs,ans);trim(out);return out;
	}
};
struct Step{SideTree side;int hs=0;};
struct Chain
{
	struct Node{int l=0,r=0,d=0;poly p,q,b;};
	vector<Step> st;
	vector<Node> tr;
	const Grid& g;
	int& ans;
	Chain(vector<Step> s,const Grid& g,int& ans):st(move(s)),tr(4*st.size()),g(g),ans(ans)
	{build(1,0,st.size());}
	void build(int u,int l,int r)
	{
		auto& z=tr[u];z.l=l;z.r=r;
		if(r-l==1)
		{
			const auto& s=st[l].side;
			if(s.root<0){z.p=z.q={1};return;}
			const poly& A=s.tr[s.root].a;
			int d=z.d=A.size()-1;z.b.resize(d);z.p.resize(d+1);z.q.resize(d+1);
			poly a(A.begin()+1,A.end()),v=g.cauchy(a,d+1,l+1);
			for(int j=0;j<=d;j++)
			{
				z.q[j]=mul(g.fac[j+l+d],g.ifac[j+l]);
				z.p[j]=mul(mul(j+l,z.q[j]),v[j]);
			}
			return;
		}
		int m=(l+r)/2;build(u*2,l,m);build(u*2+1,m,r);
		z.d=tr[u*2].d+tr[u*2+1].d;z.b.resize(z.d);
		int n=z.d+1;
		z.p=g.eval(tr[u*2].p,n);z.q=g.eval(tr[u*2].q,n);
		poly p=g.eval(tr[u*2+1].p,n),q=g.eval(tr[u*2+1].q,n);
		for(int i=0;i<n;i++)z.p[i]=mul(z.p[i],p[i]),z.q[i]=mul(z.q[i],q[i]);
	}
	poly feed(int u,poly f)
	{
		trim(f);if(f.empty())return {};
		auto& z=tr[u];int H=f.size()-1,len=z.r-z.l;
		poly out(f.size()+len);
		if(!z.d){copy(f.begin(),f.end(),out.begin()+len);return out;}
		poly p=g.eval(z.p,H,1-z.l),q=g.eval(z.q,H,1-z.l),w(H),pre(H+1,1);
		for(int i=0;i<H;i++)pre[i+1]=mul(pre[i],q[i]);
		int x=power(pre[H],mod-2);
		for(int i=H-1;i>=0;i--)
		{
			w[i]=mul(f[i+1],mul(pre[i],x));x=mul(x,q[i]);
			out[i+1+len]=mul(w[i],p[i]);
		}
		append(z.b,g.adjoint(w,z.d,1-z.l));trim(out);return out;
	}
	void push(int u)
	{
		auto& z=tr[u];int d=z.d;
		if(!d)return;
		for(int k=0;k<2;k++)
		{
			int v=u*2+k;
			if(!tr[v].d)continue;
			poly a=g.eval(k?tr[u*2].p:tr[u*2+1].q,d);
			for(int i=0;i<d;i++)a[i]=mul(a[i],z.b[i]);
			append(tr[v].b,g.adjoint(a,tr[v].d));
		}
		poly().swap(z.b);
	}
	poly solve(int u)
	{
		auto& z=tr[u];if(!z.d)return {};
		if(z.r-z.l==1)
		{
			poly w(z.d);
			for(int j=0;j<z.d;j++)w[j]=mul(z.b[j],z.q[j]);
			poly r=g.cauchy(w,z.d,z.l+1);r.insert(r.begin(),0);
			return st[z.l].side.apply(r,st[z.l].hs,ans);
		}
		push(u);
		poly a=solve(u*2),b=feed(u*2+1,move(a)),c=solve(u*2+1);
		append(b,c);trim(b);return b;
	}
	poly run()
	{
		poly a=feed(1,{0,1}),b=solve(1);append(a,b);trim(a);return a;
	}
};
struct Solver
{
	int n,ans=0;
	Grid grid;
	vector<vector<int>> e,ch;
	vector<int> sz,h,son;
	vector<poly> f;
	Solver(int n):n(n),grid(4*n+10),e(n+1),ch(n+1),sz(n+1,1),h(n+1,1),son(n+1),f(n+1){}
	int run()
	{
		vector<int> ord={1},fa(n+1);
		for(int i=0;i<n;i++)
		{
			int u=ord[i];
			for(int v:e[u])if(v!=fa[u])fa[v]=u,ch[u].push_back(v),ord.push_back(v);
		}
		for(int i=n-1;i>=0;i--)
		{
			int u=ord[i];
			for(int v:ch[u])
			{
				sz[u]+=sz[v];
				if(!son[u]||h[v]>h[son[u]])son[u]=v;
			}
			if(son[u])h[u]=h[son[u]]+1;
		}
		for(int i=n-1;i>=0;i--)
		{
			int head=ord[i];if(fa[head]&&son[fa[head]]==head)continue;
			vector<int> path;
			for(int u=head;son[u];u=son[u])path.push_back(u);
			if(path.empty()){f[head]={0,1};continue;}
			reverse(path.begin(),path.end());vector<Step> st;
			for(int u:path)
			{
				vector<int> side;
				for(int v:ch[u])if(v!=son[u])side.push_back(v);
				Step s;s.hs=sz[son[u]];s.side.init(side,f,sz);st.push_back(move(s));
			}
			Chain c(move(st),grid,ans);f[head]=c.run();f[head].resize(h[head]+1);
		}
		return ans;
	}
};
int main()
{
	ios::sync_with_stdio(false);cin.tie(nullptr);
	int c,T;cin>>c>>T;
	while(T--)
	{
		int n;cin>>n;Solver s(n);
		for(int i=1,u,v;i<n;i++)cin>>u>>v,s.e[u].push_back(v),s.e[v].push_back(u);
		cout<<s.run()<<'\n';
	}
}
```

## 附录 B：实测优化版完整代码

```cpp
#include<bits/stdc++.h>
using namespace std;
using ll=long long;
using poly=vector<int>;
constexpr int mod=998244353;
int add(int x,int y){return x+y>=mod?x+y-mod:x+y;}
int mul(int x,int y){return (ll)x*y%mod;}
int det_scalar(const int* a,const int* b,int* c,int n,int s)
{
	unsigned long long sum=0;
	for(int i=0;i<n;i++)
	{
		int z=mul(a[i],b[i]);sum+=z;
		int w=s==1?z:s==2?add(z,z):mul(s,z);
		c[i]=a[i]>=w?a[i]-w:a[i]-w+mod;
	}
	return sum%mod;
}
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include<immintrin.h>
__attribute__((target("avx2")))
inline __m128i mul4(__m128i a,__m128i b)
{
	__m256d x=_mm256_cvtepi32_pd(a),y=_mm256_cvtepi32_pd(b);
	__m128i q=_mm256_cvttpd_epi32(_mm256_mul_pd(_mm256_mul_pd(x,y),_mm256_set1_pd(1.0/mod)));
	__m256i p=_mm256_set1_epi64x(mod);
	__m256i z=_mm256_sub_epi64(_mm256_mul_epu32(_mm256_cvtepu32_epi64(a),_mm256_cvtepu32_epi64(b)),_mm256_mul_epu32(_mm256_cvtepu32_epi64(q),p));
	z=_mm256_add_epi64(z,_mm256_and_si256(_mm256_cmpgt_epi64(_mm256_setzero_si256(),z),p));
	z=_mm256_sub_epi64(z,_mm256_and_si256(_mm256_cmpgt_epi64(z,_mm256_set1_epi64x(mod-1)),p));
	return _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(z,_mm256_setr_epi32(0,2,4,6,0,0,0,0)));
}
__attribute__((target("avx2")))
int det_avx2(const int* a,const int* b,int* c,int n,int s)
{
	__m256i acc=_mm256_setzero_si256();__m128i p=_mm_set1_epi32(mod),ss=_mm_set1_epi32(s);int i=0;
	for(;i+4<=n;i+=4)
	{
		__m128i x=_mm_loadu_si128((const __m128i*)(a+i)),y=_mm_loadu_si128((const __m128i*)(b+i));
		__m128i z=mul4(x,y),w=z;
		if(s==2)
		{
			w=_mm_add_epi32(z,z);
			w=_mm_sub_epi32(w,_mm_and_si128(_mm_cmpgt_epi32(w,_mm_set1_epi32(mod-1)),p));
		}
		else if(s!=1)w=mul4(z,ss);
		__m128i r=_mm_sub_epi32(x,w);
		r=_mm_add_epi32(r,_mm_and_si128(_mm_cmpgt_epi32(_mm_setzero_si128(),r),p));
		_mm_storeu_si128((__m128i*)(c+i),r);acc=_mm256_add_epi64(acc,_mm256_cvtepu32_epi64(z));
	}
	alignas(32) unsigned long long tmp[4];_mm256_store_si256((__m256i*)tmp,acc);
	unsigned long long sum=tmp[0]+tmp[1]+tmp[2]+tmp[3];
	for(;i<n;i++)
	{
		int z=mul(a[i],b[i]);sum+=z;int w=s==1?z:s==2?add(z,z):mul(s,z);
		c[i]=a[i]>=w?a[i]-w:a[i]-w+mod;
	}
	return sum%mod;
}
#endif
int det_update(const int* a,const int* b,int* c,int n,int s)
{
#if defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
	static bool avx=__builtin_cpu_supports("avx2");
	if(avx&&n>=32)return det_avx2(a,b,c,n,s);
#endif
	return det_scalar(a,b,c,n,s);
}
void ntt(poly& a,bool inv)
{
	int n=a.size();static poly rt={0,1};
	if((int)rt.size()<n)
	{
		int k=__builtin_ctz((unsigned)rt.size());rt.resize(n);
		for(;1<<k<n;k++)
		{
			int z=1,b=3,e=(mod-1)>>(k+1);
			for(;e;e>>=1,b=mul(b,b))if(e&1)z=mul(z,b);
			for(int i=1<<(k-1);i<1<<k;i++)rt[i*2]=rt[i],rt[i*2+1]=mul(rt[i],z);
		}
	}
	for(int i=1,j=0;i<n;i++)
	{
		int k=n>>1;for(;j&k;k>>=1)j^=k;j^=k;
		if(i<j)swap(a[i],a[j]);
	}
	for(int k=1;k<n;k<<=1)for(int i=0;i<n;i+=k*2)for(int j=0;j<k;j++)
	{
		int x=a[i+j],y=mul(a[i+j+k],rt[k+j]);
		a[i+j]=add(x,y);a[i+j+k]=x>=y?x-y:x-y+mod;
	}
	if(inv)
	{
		reverse(a.begin()+1,a.end());int z=mod-(mod-1)/n;
		for(int& x:a)x=mul(x,z);
	}
}
poly multiply(const poly& a,const poly& b)
{
	int n=a.size(),m=b.size();poly c(n+m-1);
	if(min(n,m)<=8||(ll)n*m<=2048)
	{
		if(n==1&&a[0]==1)return b;
		if(m==1&&b[0]==1)return a;
		for(int k=0;k<n+m-1;k++)
		{
			int l=max(0,k-m+1),r=min(n-1,k);
			for(int i=l;i<=r;i+=16)
			{
				unsigned long long sum=0;
				for(int j=i;j<=min(r,i+15);j++)sum+=(unsigned long long)a[j]*b[k-j];
				c[k]=add(c[k],sum%mod);
			}
		}
		return c;
	}
	int z=1;while(z<n+m-1)z<<=1;
	poly x=a,y=b;x.resize(z);y.resize(z);ntt(x,0);ntt(y,0);
	for(int i=0;i<z;i++)x[i]=mul(x[i],y[i]);
	ntt(x,1);x.resize(n+m-1);return x;
}
int power(int x,int e)
{
	int y=1;for(;e;e>>=1,x=mul(x,x))if(e&1)y=mul(y,x);return y;
}
poly inverse_series(const poly& a,int n)
{
	poly b={power(a[0],mod-2)};
	while((int)b.size()<n)
	{
		int m=min(n,(int)b.size()*2);
		poly x(a.begin(),a.begin()+min(m,(int)a.size()));
		x=multiply(x,b);x.resize(m);
		for(int& v:x)if(v)v=mod-v;
		x[0]=add(x[0],2);b=multiply(b,x);b.resize(m);
	}
	return b;
}
struct Dist{int lo=1;poly a={1};};
struct Node{Dist f;int l=-1,r=-1,sz=0;};
struct Solver
{
	int n,ans=0;
	vector<vector<int>> e,ch;
	vector<int> sz,h,inv,fac,ifac,son;
	vector<Dist> f;
	vector<Node> tr;
	Solver(int n):n(n),e(n+1),ch(n+1),sz(n+1,1),h(n+1,1),inv(n+1),fac(2*n+1,1),ifac(2*n+1,1),son(n+1),f(n+1)
	{
		inv[1]=1;for(int i=2;i<=n;i++)inv[i]=mod-mul(mod/i,inv[mod%i]);
		for(int i=1;i<=2*n;i++)fac[i]=mul(fac[i-1],i);
		ifac[2*n]=power(fac[2*n],mod-2);
		for(int i=2*n;i;i--)ifac[i-1]=mul(ifac[i],i);
	}
	int build(const vector<int>& v,int l,int r)
	{
		int u=tr.size();tr.emplace_back();
		if(r-l==1){tr[u].f=move(f[v[l]]);tr[u].sz=sz[v[l]];}
		else
		{
			int m=(l+r)/2,x=build(v,l,m),y=build(v,m,r);
			tr[u].l=x;tr[u].r=y;tr[u].f.lo=tr[x].f.lo+tr[y].f.lo;
			tr[u].f.a=multiply(tr[x].f.a,tr[y].f.a);
		}
		return u;
	}
	poly trans(const poly& a,const poly& b,int n)
	{
		poly r(n);int m=b.size();
		if((ll)n*m<=2048||m==1)
		{
			for(int i=0;i<n;i++)for(int l=0;l<m;l+=16)
			{
				unsigned long long sum=0;
				for(int j=l;j<min(m,l+16);j++)sum+=(unsigned long long)a[i+j]*b[j];
				r[i]=add(r[i],sum%mod);
			}
			return r;
		}
		poly c=multiply(a,poly(b.rbegin(),b.rend()));
		for(int i=0;i<n;i++)r[i]=c[i+m-1];
		return r;
	}
	void push(int u,const poly& r,Dist& out,int hs)
	{
		auto& z=tr[u];
		if(z.l<0)
		{
			int p=0;
			for(int i=0;i<(int)r.size();i++)
			{
				int k=z.f.lo+i,x=mul(mul(k,z.f.a[i]),r[i]);p=add(p,x);
				out.a[k+1-out.lo]=add(out.a[k+1-out.lo],x);
			}
			ans=add(ans,mul((hs-z.sz+mod)%mod,p));return;
		}
		poly a=trans(r,tr[z.r].f.a,tr[z.l].f.a.size());
		poly b=trans(r,tr[z.l].f.a,tr[z.r].f.a.size());
		push(z.l,a,out,hs);push(z.r,b,out,hs);
	}
	void batch(const vector<int>& path,int l,int r,const vector<int>& key)
	{
		int q=r-l,S=accumulate(key.begin(),key.end(),0),mx=key.back();
		Dist F=move(f[son[path[l]]]);int H=F.lo+(int)F.a.size()-1;
		poly a(F.a.size()),kernel(H+q);
		for(int i=1;i<H+q;i++)kernel[i]=mul(fac[i-1],ifac[i+S]);
		for(int i=0;i<(int)a.size();i++)
		{
			int k=F.lo+i;a[i]=mul(F.a[i],mul(fac[k+S-1],ifac[k-1]));
		}
		poly c=multiply(poly(a.rbegin(),a.rend()),kernel),A(q);
		for(int t=0;t<q;t++)A[t]=c[H+t];
		poly w(mx+1);
		for(int k:key)w[k]=add(w[k],mul(k,mul(fac[k+S],ifac[k])));
		kernel.resize(mx+q);
		for(int i=1;i<mx+q;i++)kernel[i]=mul(fac[i-1],ifac[i+S]);
		c=multiply(poly(w.rbegin(),w.rend()),kernel);poly B(q);B[0]=1;
		for(int j=1;j<q;j++)B[j]=c[mx+j]?mod-c[mx+j]:0;
		poly R=multiply(A,inverse_series(B,q));R.resize(q);
		int lo=min(F.lo+q,key[0]+1),hi=H+q;
		Dist out;out.lo=lo;out.a.assign(hi-lo+1,0);
		for(int i=0;i<(int)a.size();i++)
		{
			int k=F.lo+i+q;out.a[k-lo]=mul(a[i],mul(fac[k-1],ifac[k+S-1]));
		}
		c=multiply(w,poly(R.rbegin(),R.rend()));
		for(int k=key[0]+1;k<=mx+q;k++)
			out.a[k-lo]=add(out.a[k-lo],mul(c[k-1],mul(fac[k-1],ifac[k+S-1])));
		for(int t=0;t<q;t++)
		{
			int u=path[l+t],hs=sz[son[u]],sum=0;ll v=0;
			for(int x:ch[u])if(x!=son[u])sum+=sz[x],v+=(ll)(hs-sz[x])*f[x].lo,f[x].a.clear();
			int z=v%mod;if(z<0)z+=mod;ans=add(ans,add(sum,mul(z,R[t])));
		}
		f[path[r-1]]=move(out);
	}

	void step(int u)
	{
		int heavy=son[u];
		if(ch[u].size()==1){f[u]=move(f[heavy]);f[u].lo++;return;}

			vector<int> side;int lo=f[heavy].lo+1;
			for(int v:ch[u])if(v!=heavy)side.push_back(v),lo=min(lo,f[v].lo+1),ans=add(ans,sz[v]);
			bool det=true;
			for(int v:side)if(f[v].a.size()!=1)det=false;
			if(det)
			{
				int s=0;ll weighted=0;
				for(int v:side)s+=f[v].lo,weighted+=(ll)sz[v]*f[v].lo;
				Dist out;out.lo=lo;out.a.assign(h[u]-lo+1,0);
				const auto& F=f[heavy];
				int r=det_update(F.a.data(),inv.data()+F.lo+s,out.a.data()+F.lo+1-lo,F.a.size(),s);
				int c=((ll)sz[heavy]*s-weighted)%mod;
				if(c<0)c+=mod;
				ans=add(ans,mul(c,r));
				for(int v:side)
				{
					int k=f[v].lo;
					out.a[k+1-lo]=add(out.a[k+1-lo],mul(k,r));
					f[v].a.clear();
				}
				f[u]=move(out);poly().swap(f[heavy].a);return;
			}
			tr.clear();tr.reserve(side.size()*2);int root=build(side,0,side.size());
			const auto& A=tr[root].f;const auto& F=f[heavy];int H=F.a.size(),D=A.a.size(),b=F.lo+A.lo;
			poly r(D);
			Dist out;out.lo=lo;out.a.assign(h[u]-lo+1,0);
			if((ll)H*D<=80000||min(H,D)<=4)
			{
				vector<unsigned long long> sums(D);
				for(int i=0;i<H;i++)
				{
					int total=0;
					for(int l=0;l<D;l+=16)
					{
						unsigned long long sum=0;
						for(int j=l;j<min(D,l+16);j++)
						{
							int z=mul(F.a[i],inv[b+i+j]);sums[j]+=z;
							sum+=(unsigned long long)z*A.a[j];
						}
						total=add(total,sum%mod);
					}
					out.a[F.lo+i+1-lo]=mul(F.lo+i,total);
				}
				for(int j=0;j<D;j++)r[j]=sums[j]%mod;
			}
			else
			{
				poly k(H+D-1);for(int i=0;i<H+D-1;i++)k[i]=inv[b+i];
				poly x=multiply(poly(F.a.rbegin(),F.a.rend()),k),y=multiply(poly(A.a.rbegin(),A.a.rend()),k);
				for(int j=0;j<D;j++)r[j]=x[H-1+j];
				for(int i=0;i<H;i++)out.a[F.lo+i+1-lo]=mul(mul(F.lo+i,F.a[i]),y[D-1+i]);
			}
			push(root,r,out,sz[heavy]);f[u]=move(out);poly().swap(f[heavy].a);
	}
	int run()
	{
		vector<int> ord={1},fa(n+1);
		for(int i=0;i<n;i++)for(int v:e[ord[i]])if(v!=fa[ord[i]])fa[v]=ord[i],ch[ord[i]].push_back(v),ord.push_back(v);
		for(int i=n-1;i>=0;i--)
		{
			int u=ord[i];for(int v:ch[u]){sz[u]+=sz[v];if(!son[u]||h[v]>h[son[u]])son[u]=v;}
			if(son[u])h[u]=h[son[u]]+1;
		}
		for(int i=n-1;i>=0;i--)
		{
			int u=ord[i];if(fa[u]&&son[fa[u]]==u)continue;
			vector<int> path;for(int x=u;son[x];x=son[x])path.push_back(x);
			reverse(path.begin(),path.end());int q=path.size();vector<vector<int>> keys(q);
			for(int t=0;t<q;t++)
			{
				int x=path[t];bool ok=true;
				for(int v:ch[x])if(v!=son[x])
				{
					if(f[v].a.size()!=1){ok=false;break;}keys[t].push_back(f[v].lo);
				}
				if(!ok)keys[t].clear();else sort(keys[t].begin(),keys[t].end());
			}
			for(int l=0;l<q;)
			{
				int r=l+1;
				if(!keys[l].empty())while(r<q&&keys[r]==keys[l])r++;
				if(!keys[l].empty()&&r-l>=128&&(ll)(r-l)*(f[son[path[l]]].lo+(int)f[son[path[l]]].a.size()-1)+(ll)(r-l)*(r-l)/2>=200000)batch(path,l,r,keys[l]);else for(int t=l;t<r;t++)step(path[t]);
				l=r;
			}
		}
		return ans;
	}

};
int main()
{
	ios::sync_with_stdio(false);cin.tie(nullptr);int c,T;cin>>c>>T;
	while(T--){int n;cin>>n;Solver s(n);for(int i=1,u,v;i<n;i++)cin>>u>>v,s.e[u].push_back(v),s.e[v].push_back(u);cout<<s.run()<<'\n';}
	return 0;
}
```
