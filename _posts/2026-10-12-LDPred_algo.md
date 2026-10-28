---
layout: blog_post
title: 'LDPred Algorithm Note'
date: 2026-10-12
tags:
  - Learning
  - Algorithm
---

Learning notes about LDPred algorithm.




# 1. Principle Analysis

---
## Model & Overview

**Linear Model**

Suppose the GWAS summary statistics are calculated by standardized phenotype $y \in \mathbb{R}^{N}$ and genotype $X \in \mathbb{R}^{N \times M}$

$$
\tilde{\beta}
= \frac{X^\top y}{N}.
$$

The linear model assumption implies $y$ and $X$ have the following relation:
$$
y = X \beta + \varepsilon, \ \ \ \varepsilon \sim \mathcal{N}\left(0,(1-h^2)I\right).
$$

Hence the GWAS summary statistics satisfies the following form under linear model assumption:

$$
\begin{aligned}
\tilde{\beta}
&= \frac{X^\top y}{N} \\
&= \frac{X^\top X}{N}\beta + \frac{X^\top \varepsilon}{N} \\
&= D\beta + \frac{X^\top \varepsilon}{N}
\end{aligned}
$$

where $D = \frac{X^\top X}{N}$ is the LD matrix.


Under a marginal view,

$$
\boxed{
\tilde{\beta}_j = D_{jj}\beta_j + \sum_{k \ne j}D_{jk} \beta_k + \frac{1}{N} X_j^\top \varepsilon
}
\tag{1}
$$


However, the value of $X$ and $y$ are often unavaliable, we can only observe GWAS summary statistics $\tilde{\beta}$. The rough idea of LDPred is:
- First replacing $X$ by a reachable reference panel $X_r$. This assumes both $X$ and $X_r$ are samples from the same population so that $X^TX$ and $X_r^TX_r$ share the same structure.
- Then using equation (1) to infer $\beta$ given $\tilde{\beta}, D_r$.

**In the following text, we use $D$ to represent the LD matrix calculated by reference panel.**




The design of $(1-h^2)I$ is because phenotypic variance is normalized to $1$ for each individual. If the genotype component explains variance $h^2$, then we have

$$
\operatorname{Var}(y_n) = \operatorname{Var}((X\beta)_n) + \operatorname{Var}(\varepsilon_n) = h^2 + (1 - h^2) = 1.
$$


**Bayesian Framework**

Given a prior distribution of $f(\beta)$, equation $(1)$ allows us to derive $f(\tilde{\beta} \mid \beta)$. Bayesian method suggests


$$
f(\beta \mid \tilde{\beta}) \propto
\underbrace{f(\beta)}_{\substack{\text{Initial belief on SNP effect} \\ \text{ (Sparse / Dense)}}}
\cdot
\underbrace{f(\tilde{\beta} \mid \beta)}_{\substack{\text{Likelihood of observed data} \\ \text{ (Linked / Unlinked)}}}
$$

Each components of posterior distribution is influenced by 2 factors, to illustrate this algorithm, we analized three scenarios:
- Dense + Unlinked
- Sparse + Unlinked
- Sparse + Linked


The final goal is using posterior expectation $E(\beta \mid \tilde{\beta})$ as estimated effect weights, the general approach to this is:



$$
\begin{array}{c}
\boxed{\text{Set Dense/Sparse + Unlinked/Linked Assumption}} \\[4pt]
\downarrow \\[4pt]
\boxed{\text{Determine Prior } f(\beta) \text{ and Likelihood } f(\tilde{\beta} \mid \beta)} \\[4pt]
\downarrow \\[4pt]
\boxed{\text{Find Posterior } f(\beta \mid \tilde{\beta})} \\[4pt]
\downarrow \\[4pt]
\boxed{\text{Is } \mathbb{E}(\beta \mid \tilde{\beta}) \text{ analytic?}} \\[8pt]
\begin{array}{cc}
\text{Yes} & \text{No} \\[2pt]
\downarrow & \downarrow \\[2pt]
\boxed{\text{Direct computation}} & \boxed{\text{Sample } f(\beta \mid \tilde{\beta}) \text{ and average}}
\end{array}
\end{array}
$$




**Sparse/Dense Assumption**

- Dense implies every SNP matters

$$
\beta_j
\overset{\text{i.i.d.}}{\sim}
\mathcal{N}\left(0, \frac{h^2}{M}\right)
$$

- Sparsity assumes each SNP may or may not be causal with the same probability $p$

$$
\begin{aligned}
\beta_j
&\overset{\text{i.i.d.}}{\sim}
\begin{cases}
\mathcal{N}\left(0, \dfrac{h^2}{Mp}\right), & \text{with probability } p \\[6pt]
0, & \text{with probability } 1-p
\end{cases} \\[10pt]
&\hspace{2em}= (1-p)\,\mathbf{1}_{(\beta_j = 0)} + p\,\mathcal{N}\left(0, \dfrac{h^2}{Mp}\right)
\end{aligned}
$$


**Linked/Unlinked Assumption**

Recall that

$$
\tilde{\beta}_j
=
\underbrace{D_{jj}\beta_j}_{\text{True signal}}
+
\underbrace{\sum_{k \ne j}D_{jk} \beta_k}_{\text{LD leakage}}
+
\underbrace{\frac{1}{N} X_j^\top \varepsilon}_{\text{Sampling noise}}
$$

- Under the unlinked assumption, LD leakage does not exist. For the given standardized sample $X$, $D_{jj}=1$ holds, and $\sum_{k \ne j}D_{jk}\beta_k$ can also be regarded as random sampling noise, i.e.

$$
\tilde{\beta}_j=\beta_j + \text{noise}_j,  \ \ \ \text{noise}_j \sim \mathcal{N}\left(0, \dfrac{1}{N}\right)
$$




- Under the linked assumption, LD leakage is assumed to exist only between adjacent SNPs. The paper further partitions SNPs into blocks by a given window size and **here we assume that all SNPs within the same block share the same neighborhood**. Suppose $\beta_j$ is attributed into $l$-th block $\mathcal{B}_l$, then

    $$
    \tilde{\beta}_j = \beta_j + \sum_{k \in \mathcal{B}_l,k\ne j} D_{jk} \beta_k + \frac{1}{N} X_j^\top \eta^{(l)},  \ \ \ \eta^{(l)} \sim \mathcal{N}\left(0, \dfrac{1- h_l^2}{N}\right)
    $$

    The error term $\eta^{(l)}$ differs across blocks, since it includes LD leakage from outside the block. And $h_l^2$ can be interpreted as heritability explained by the effect weights inside the block.

    In matrix form, For the $l$-th LD block $\mathcal{B}_l$, denote the true effects in this region by $\beta^{(l)}$, the marginal estimates by $\tilde{\beta}^{(l)}$, the genotype SNPs $X^{(l)}$, and the LD matrix by $D^{(l)}$.


    $$
    \tilde{\beta}^{(l)} = D^{(l)}\beta^{(l)} + \frac{(X^{(l)})^{\top}\eta^{(l)}}{N}
    $$

---
## (i) Dense + Unlinked Scenario

This represents the most idealized and simplest setting.

**Dense prior:**

$$
\beta_j
\sim
\mathcal{N}\left(0,\frac{h^2}{M}\right)
$$

**Unlinked likelihood:**

$$
\tilde{\beta}_j\mid\beta_j
\sim
\mathcal{N}\left(\beta_j,\frac{1}{N}\right)
$$

**Posterior:**

The above yields a classic conjugate-normal case, in which the posterior admits a closed-form solution

$$
\beta_j\mid\tilde{\beta}_j
\sim
\mathcal{N}\left(
\frac{1}{1+\frac{M}{Nh^2}}\tilde{\beta}_j,
\frac{1}{N+\frac{M}{h^2}}
\right)
$$

with posterior mean

$$
E(\beta_j\mid\tilde{\beta}_j)
=
\frac{1}{1+\frac{M}{Nh^2}}\tilde{\beta}_j
$$

The implementation is simple: apply the same shrinkage factor independently to each SNP.


---
## (ii) Sparse + Unlinked Scenario


**Sparse prior:**

$$
\beta_j
\sim
(1-p)\mathbf{1}_{(\beta_j = 0)}
+
p\mathcal{N}\left(0,\frac{h^2}{Mp}\right).
$$

**Unlinked likelihood:**

$$
\tilde{\beta}_j\mid\beta_j
\sim
\mathcal{N}\left(\beta_j,\frac{1}{N}\right).
$$

**Posterior:**

The posterior can be directly constructed by

$$
f(\beta_j | \tilde{\beta_j})
=
\frac{
f(\tilde{\beta}_i \mid \beta_i) f(\beta_i)
}{
\int_{-\infty}^{\infty} f(\tilde{\beta}_i \mid \beta_i) f(\beta_i) \, d\beta_i
}
$$


However, the paper provides an alternative way,

$$
\begin{aligned}
f(\beta_j \mid \tilde{\beta}_j)
&= P(\beta_j = 0 \mid \tilde{\beta}_j) \cdot f(\beta_j \mid \beta_j = 0, \tilde{\beta}_j) \\
&\quad + P(\beta_j \ne 0 \mid \tilde{\beta}_j) \cdot f(\beta_j \mid \beta_j \ne 0, \tilde {\beta}_j) \\
\end{aligned}
\tag{2}
$$

Next we derive each term of above equation.

- $f(\beta_j \mid \beta_j=0, \tilde{\beta}_j)$

$$
\begin{aligned}
f(\beta_j \mid \beta_j=0, \tilde{\beta}_j)
&= \frac{
        \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_j = 0)}
        \cdot \color{#15803d}{f(\beta_j | \beta_j=0)}
        \cdot \color{black}{P(\beta_j=0)}}
    {
        \color{#1d4ed8}{f(\tilde{\beta}_j | \beta_j=0) }\cdot \color{black}{P(\beta_j=0)}
    } \\
&= \frac{
        \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j = 0)}
        \cdot \color{#15803d}{f(\beta_j | \beta_j=0)}
    }{
        \color{#1d4ed8}{f(\tilde{\beta}_j | \beta_j=0)}
    } \\
&= \color{#15803d}{f(\beta_j | \beta_j=0)} \\
&= \mathbf{1}_{(\beta_j = 0)}
\end{aligned}
$$


- $f(\beta_j \mid \beta_j\ne0, \tilde{\beta}_j)$

$$
\begin{aligned}
f(\beta_j \mid \beta_j\ne0, \tilde{\beta}_j)
&= \frac{
        \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_j \ne 0)}
        \cdot \color{#15803d}{f(\beta_j | \beta_j\ne0)}
        \cdot \color{black}{P(\beta_j\ne0)}}
    {
        \color{#1d4ed8}{f(\tilde{\beta}_j | \beta_j\ne0) }\cdot \color{black}{P(\beta_j\ne0)}
    } \\
&= \frac{
        \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_j \ne 0)}
        \cdot \color{#15803d}{f(\beta_j | \beta_j\ne0)}
    }{
        \color{#1d4ed8}{f(\tilde{\beta}_j | \beta_j\ne0)}
    } \\
&= \color{#1d4ed8}{\text{const}}
        \cdot \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_j \ne 0)}
        \cdot \color{#15803d}{f(\beta_j | \beta_j\ne0)}
\end{aligned}
$$

&emsp; &emsp; &emsp; &emsp;The rest 2 density functions are known normal pdf:

$$
(\tilde{\beta}_j \mid \beta_j, \beta_j \ne 0) \sim \mathcal{N}\left(\beta_j,\frac{1}{N}\right)
$$

$$
(\beta_j \mid \beta_j \ne 0) \sim \mathcal{N}\left(0,\frac{h^2}{Mp}\right)
$$


&emsp; &emsp; &emsp; &emsp;This is also a conjugate-normal case highly similar to senario 1, we can directly get

$$
(\beta_j \mid \beta_j\ne0, \tilde{\beta}_j)
\sim
\mathcal{N}\left(
\frac{1}{1+\frac{Mp}{Nh^2}}\tilde{\beta}_j,
\frac{1}{N+\frac{Mp}{h^2}}
\right)
$$


- $P(\beta_j \ne 0 \mid \tilde{\beta}_j) $

    &emsp; &emsp; We simply denote this term by $\bar{p}_j$ , and recall the assumptions of current cenario:

    $$
    \tilde{\beta}_j=\beta_j + \text{noise}_j,  \ \ \ \text{noise}_j \sim \mathcal{N}\left(0, \dfrac{1}{N}\right) = \mathcal{N}\left(0, \sigma_e^2\right)
    $$

    $$(\beta_j | \beta_j \ne 0) \sim  \mathcal{N}\left(0, \dfrac{h^2}{Mp}\right)= \mathcal{N}\left(0, \sigma_b^2\right)$$

    &emsp; &emsp; Then we can derive that

    $$
    \begin{aligned}
    (\tilde{\beta}_j \mid \beta_j \ne 0) &\sim \mathcal{N}\left(0,\ \sigma_e^2 + \sigma_b^2\right) \\
    (\tilde{\beta}_j \mid \beta_j = 0) &\sim \mathcal{N}\left(0,\ \sigma_e^2\right)
    \end{aligned}
    $$

    &emsp; &emsp; Back to the goal of writing $\bar{p}_j$

    $$
    \begin{aligned}
    \bar{p}_j
    =
    P(\beta_j\ne0\mid\tilde{\beta}_j)
    &=
        \frac{
        P(\beta_j \ne 0) \cdot P(\tilde{\beta_j}|\beta_j \ne 0)
        }{
        P(\beta_j \ne 0) \cdot P(\tilde{\beta_j}|\beta_j \ne 0)
        +
        P(\beta_j = 0) \cdot P(\tilde{\beta_j}|\beta_j = 0)
        } \\[18pt]
    &=
        \frac{
        p\phi(\tilde{\beta}_j;0,\sigma_b^2+\sigma_e^2)
        }{
        p\phi(\tilde{\beta}_j;0,\sigma_b^2+\sigma_e^2)
        +
        (1-p)\phi(\tilde{\beta}_j;0,\sigma_e^2)
        } \\
    \end{aligned}
    $$

- $P(\beta_j = 0 \mid \tilde{\beta}_j) = 1 - \bar{p}_j$




Plugging the above four results into the posterior in $(2)$, we get

$$
\begin{aligned}
f(\beta_j \mid \tilde{\beta}_j)
&= (1 - \bar{p}_j) \cdot f(\beta_j \mid \beta_j = 0, \tilde{\beta}_j)
 + \bar{p}_j \cdot f(\beta_j \mid \beta_j \ne 0, \tilde {\beta}_j) \\
\end{aligned}
$$



Therefore, the posterior expectation is:

$$
\begin{aligned}
E[\beta_j \mid \tilde{\beta}_j]
&= \int \beta_j \left[
(1 - \bar{p}_j) \, f(\beta_j \mid \beta_j = 0, \tilde{\beta}_j)
+
\bar{p}_j \, f(\beta_j \mid \beta_j \ne 0, \tilde{\beta}_j)
\right] d\beta_j \\
&= (1 - \bar{p}_j)\int \beta_j \cdot \mathbf{1}_{\{\beta_j = 0\}} \, d\beta_j
    +
    \bar{p}_j \int \beta_j \cdot \phi(\beta_j; \mu, \sigma^2) \, d\beta_j \\
&= 0 + \frac{\bar{p}_j}{1+\frac{Mp}{Nh^2}}\tilde{\beta}_j \\
&= \frac{\bar{p}_j}{1+\frac{Mp}{Nh^2}}\tilde{\beta}_j
\end{aligned}
$$



Compared with first scenario, the only additional key quantity is $\bar{p}_j$, the posterior causal probability.

---
## (iii) Sparse + Linked Scenario

This is the core LDpred scenario.
**Sparse Prior:**

$$
\beta_j
\overset{\text{i.i.d.}}{\sim}
(1-p)\mathbf{1}_{(\beta_j = 0)}
+
p\mathcal{N}\left(0,\frac{h^2}{Mp}\right).
$$

$$
f(\beta) = \prod_{j=1}^{M} \left( p \sqrt{\frac{Mp}{2\pi h^2}} \exp\left\{ -\frac{Mp}{2h^2} \beta_j^2 \right\} + (1-p)\mathbf{1}_{(\beta_j = 0)} \right)
$$

**Linked Likelihood:**

We simplified the paper by assuming that SNP linkage only appears within a local LD block, for the $l$-th block $\mathcal{B}_l$,

$$
\begin{aligned}
\tilde{\beta}^{(l)} = D^{(l)}\beta^{(l)} + \frac{(X^{(l)})^\top\eta^{(l)}}{N}.
\end{aligned}
$$


Under this view, likelihood is:


$$
\tilde{\beta}^{(l)} |  \beta^{(l)} \sim \mathcal{N}\left(D^{(l)}\beta^{(l)}, \frac{1-h_l^2}{N}D^{(l)}\right) .
$$

$$
\tilde{\beta}_j^{(l)} \mid \beta^{(l)}
\sim
\mathcal{N}\left(
 \beta_j^{(l)}+\sum_{k\in\mathcal{B}_l,k\ne j}D_{jk}^{(l)}\beta_k^{(l)},\
\frac{1-h_l^2}{N}
\right)
$$

For notational cleanliness in what follows, **we drop the block-related superscript notation \((\cdot)^{(l)}\) from here on**. We hereby declare that all subsequent discussion takes place within a single block. (Alternatively, one may "imagine" that only one block is defined, i.e., the full LD matrix.)


**Posterior:**

$$
f(\beta \mid \tilde{\beta}) =  \frac{ f(\tilde{\beta} \mid \beta) f(\beta)}{\int f(\tilde{\beta} \mid \beta) f(\beta) \, d\beta}
$$

This posterior is not a common distribution, directly calculating

$$
E(\beta \mid \tilde{\beta}) = \int \frac{\beta_i f(\tilde{\beta} \mid \beta) f(\beta)}
{\int f(\tilde{\beta} \mid \beta) f(\beta) \, d\beta} \, d\beta,
\tag{3}
$$

is extremely complicated, since this expectation consists of two $M$-dimensional integrands.


**MCMC Sampling:**

A natural alternative is estimating $(3)$ by sampling from $f(\beta \mid \tilde{\beta})$ and take the average. Gibbs sampler, as a MCMC Sampling technique, suggests:

$$
\begin{aligned}
&\textbf{1. Set } \beta = \beta_{\text{init}} \\[4pt]
&\textbf{2. For } j = 1, 2, \dots, M: \\[4pt]
&\qquad \text{Sample } \beta_j^* \sim f(\beta_j \mid \beta_{-j}, \tilde{\beta}) \\[4pt]
&\qquad \text{Update }  \beta_j \leftarrow  \beta_j^*\\[4pt]
&\textbf{3. Return } \beta
\end{aligned}
$$

as one Gibbs iteration. After several burn-in rounds, sampled values would follow the posterior distribution.

So all we need to do is sampling from



$$\begin{aligned}
f(\beta_j \mid \beta_{-j}, \tilde{\beta})= \frac{f(\tilde{\beta}, \beta)}{f(\tilde{\beta}, \beta_{-j})}
&= \frac{f(\tilde{\beta} \mid \beta) f(\beta)}{f(\tilde{\beta} \mid \beta_{-j}) f(\beta_{-j})}
\overset{\beta_1, ... \beta_M \scriptstyle\text{ are i.i.d.}}{=} \frac{f(\tilde{\beta} \mid \beta) f(\beta_i)}{f(\tilde{\beta} \mid \beta_{-j})} \\
&= \frac{f(\tilde{\beta} \mid \beta) f(\beta_i)}{\int f(\tilde{\beta} \mid \beta) f(\beta_i) \, d\beta_i}.
\end{aligned}$$

However, sampling from this distribution is still not trivial. So the paper arose **an assumption for splification**:

$$
f(\beta_j \mid \beta_{-j}, \color{#15803d}{\tilde{\beta}} \color{black}{) \approx f(\beta_j \mid \beta_{-j},} \color{#15803d}{\tilde{\beta}_j} \color{black}{)}
\tag{\text{Assumption 1}}
$$


Therefore, the key point is sampling from

$$
\begin{aligned}
f(\beta_j \mid \beta_{-j}, \tilde{\beta}_j)
&= P(\beta_j = 0 \mid  \beta_{-j}, \tilde{\beta}_j) \cdot f(\beta_j \mid  \beta_{-j}, \tilde{\beta}_j,  \beta_j = 0 )\\
&\quad + P(\beta_j \ne 0 \mid  \beta_{-j}, \tilde{\beta}_j) \cdot f(\beta_j \mid  \beta_{-j}, \tilde {\beta}_j, \beta_j \ne 0 )\\
\end{aligned}
$$


Very similar to the second scenario, we derive each term of above equation.

- $f(\beta_j \mid \beta_{-j}, \tilde{\beta}_j, \beta_j = 0)$

$$
\begin{aligned}
f(\beta_j \mid \beta_{-j}, \tilde{\beta}_j, \beta_j = 0)
&= \frac{
        \color{#b91c1c}{f(\beta_{-j}, \tilde{\beta}_j \mid \beta_j, \beta_j = 0)}
        \cdot \color{#15803d}{f(\beta_j \mid \beta_j = 0)}
        \cdot \color{black}{P(\beta_j = 0)}
    }{
        \color{#1d4ed8}{f(\beta_{-j}, \tilde{\beta}_j \mid \beta_j = 0)}
        \cdot \color{black}{P(\beta_j = 0)}
    } \\
&= \frac{
        \color{#b91c1c}{f(\beta_{-j}, \tilde{\beta}_j \mid \beta_j = 0)}
        \cdot \color{#15803d}{f(\beta_j \mid \beta_j = 0)}
    }{
        \color{#1d4ed8}{f(\beta_{-j}, \tilde{\beta}_j \mid \beta_j = 0)}
    } \\
&= \color{#15803d}{f(\beta_j \mid \beta_j = 0)} \\
&= \mathbf{1}_{\{\beta_j = 0\}}.
\end{aligned}
$$



- $f(\beta_j \mid  \beta_{-j}, \tilde {\beta}_j, \beta_j \ne 0 )$

$$
\begin{aligned}
f(\beta_j \mid \beta_{-j}, \tilde{\beta}_j, \beta_j \ne 0)
&= \frac{
        \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_{-j}, \beta_j \ne 0)}
        \cdot \color{#15803d}{f(\beta_j \mid \beta_{-j}, \beta_j \ne 0)}
        \cdot \color{black}{P(\beta_j \ne 0)}
    }{
        \color{#1d4ed8}{f(\tilde{\beta}_j, \beta_{-j} \mid \beta_j \ne 0)}
        \cdot \color{black}{P(\beta_j \ne 0)}
    } \\
&= \frac{
        \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_{-j}, \beta_j \ne 0)}
        \cdot \color{#15803d}{f(\beta_j \mid \beta_{-j}, \beta_j \ne 0)}
    }{
        \color{#1d4ed8}{f(\tilde{\beta}_j, \beta_{-j} \mid \beta_j \ne 0)}
    } \\
&= \color{#1d4ed8}{\text{const}}
        \cdot \color{#b91c1c}{f(\tilde{\beta}_j \mid \beta_j, \beta_{-j}, \beta_j \ne 0)}
        \cdot \color{#15803d}{f(\beta_j \mid \beta_{-j}, \beta_j \ne 0)}.
\end{aligned}
$$

&emsp; &emsp; &emsp; &emsp;The rest 2 densitise are also known normal pdf

$$
\begin{aligned}
(\beta_j \mid \beta_{-j}, \beta_j \ne 0)
&= (\beta_j \mid \beta_j \ne 0)
\sim \mathcal{N}\left(0, \frac{h^2}{Mp}\right), \\[6pt]
(\tilde{\beta}_j \mid \beta_j, \beta_{-j}, \beta_j \ne 0)
&\sim \mathcal{N}\left(\beta_j + \sum_{k \ne j} D_{jk} \beta_k, \frac{1 - h^2}{N}\right).
\end{aligned}
$$

&emsp; &emsp; &emsp; &emsp;This is still a conjugate-normal case, we can directly get

$$
(\beta_j \mid  \beta_{-j}, \tilde {\beta}_j, \beta_j \ne 0 )
\sim
\mathcal{N}\left(
\frac{\tilde{\beta}_j - \sum_{k \ne j}D_{jk}\beta_k}{1 + \dfrac{Mp}{Nh^2}} ,\
\frac{1}{1 + \dfrac{Mp}{Nh^2}} \frac{1}{N}
\right)
$$



-  $P(\beta_j \ne 0 \mid  \beta_{-j}, \tilde{\beta}_j)$

&emsp; &emsp; &emsp; &emsp;Define $a_j = \sum_{k \ne j}D_{jk}\beta_k$, then the current scenario implies:

$$
\tilde{\beta}_j = \beta_j + a_j + \frac{1}{N}X_j^\top\varepsilon.
$$

$$
\beta_j
\overset{\text{i.i.d.}}{\sim}
(1-p)\mathbf{1}_{(\beta_j = 0)}
+
p\mathcal{N}\left(0,\frac{h^2}{Mp}\right).
$$



&emsp; &emsp; &emsp; &emsp;Let $\sigma_e^2 = \frac{1 - h^2}{N}$ and  $\sigma_b^2 = \frac{h^2}{Mp}$, then we can derive that:

$$
(\tilde{\beta}_j \mid \beta_{-j}, \beta_j=0)
\sim
\mathcal{N}(a_j,\sigma_e^2),
$$

$$
(\tilde{\beta}_j \mid \beta_{-j}, \beta_j\ne0)
\sim
\mathcal{N}(a_j,\sigma_e^2+\sigma_b^2),
$$



&emsp; &emsp; &emsp; &emsp;The target is

$$
\begin{aligned}
P(\beta_j \ne 0 \mid \beta_{-j}, \tilde{\beta}_j)
&=
\frac{
p \cdot P(\tilde{\beta}_j \mid \beta_{-j}, \beta_j \ne 0)
}{
p \cdot P(\tilde{\beta}_j \mid \beta_{-j}, \beta_j \ne 0)
+
(1-p) \cdot P(\tilde{\beta}_j \mid \beta_{-j}, \beta_j = 0)
} \\[16pt]
&=
\frac{
p\,\phi(\tilde{\beta}_j;a_j,\sigma_e^2+\sigma_b^2)
}{
p\,\phi(\tilde{\beta}_j;a_j,\sigma_e^2+\sigma_b^2)
+
(1-p)\,\phi(\tilde{\beta}_j;a_j,\sigma_e^2)
} \\[16pt]
&=
\bar{\pi}_j
\end{aligned}
$$

- $P(\beta_j = 0 \mid  \beta_{-j}, \tilde{\beta}_j) = 1-\bar{\pi}_j$

To sum up, the whole picture of steps above is:


$$
\boxed{
\begin{aligned}
& \text{Sample from } f(\beta_j \mid \beta_{-j}, \tilde{\beta}) \\[6pt]
& \downarrow \ \scriptstyle{\text{Assumption 1}} \\[6pt]
& \text{Sample from } f(\beta_j \mid \beta_{-j}, \tilde{\beta}_j) \\[6pt]
& \downarrow \ \scriptstyle{\text{}} \\[6pt]
&\text{Sample from } (1 - \bar{\pi}_j)\,\mathbf{1}_{\{\beta_j = 0\}}
+
\bar{\pi}_j\, f(\beta_j \mid \beta_{-j}, \tilde{\beta}_j, \beta_j \ne 0) \\[6pt]
& =
\begin{cases}
0, & \text{w.p. } 1 - \bar{\pi}_j, \\[6pt]
\mathcal{N}\left(
\dfrac{\tilde{\beta}_j - a_j}{1 + \dfrac{Mp}{Nh^2}},\
\dfrac{1}{1 + \dfrac{Mp}{Nh^2}} \dfrac{1}{N}
\right), & \text{w.p. } \bar{\pi}_j.
\end{cases}
\end{aligned}
}
$$

Sweeping this over $M$ times is called one Gibbs sampler iteration, after $T$ Gibbs sampler iterations, we would have $\{\beta^{[1]}, \beta^{[2]}, \dots, \beta^{[T]}\}.$ Then posterior mean can be numerically approximated by:

$$
\frac{1}{T} \sum_{t=1}^T \beta^{[t]}.
$$




---
# 2. Complexity Analysis


## (i) LD Calculation Cost

LDPred materializes $D^{(l)} = \frac{(X^{(l)})^T X^{(l)}}{N}$, in each block, suppose each block has the same size $K \times K$, the calculation of one block is:
$$
O(NK^2).
$$

And there are $\frac{M}{K}$ blocks in total, so the overall cost for LD calculation is:
$$
O(NMK)
$$

## (ii) Gibbs Sampler Cost

Here we only analyze the core method in **1.3 Sparse + Linked Scenario**, namely the time complexity of the approximate Gibbs sampler.

Notation:

- $M$: number of genome-wide SNPs.
- $T$: number of Gibbs sampler iteration, meaning that every SNP is updated once.
- $L$: number of LD neighborhood blocks.
- $\mathcal{B}_l$: $l$-th local LD neighborhood block.
- $K$: block size, or number of SNPs within a LD block, we assume all blocks are of the same size.


For each Gibbs iteration, LDpred iterates over SNPs. Considering the update for a given SNP $j$, the steps required are:

| Step | Operation for SNP $j$ | Formula / object | Time Complexity |
|---|---|---|---|
| 1 | Compute the LD contribution from other SNPs | $a_j = \sum_{k \in \mathcal{B}_l,\ k\ne j}D_{jk}\beta_k$ | $O(K)$ |
| 2 | Compute the posterior causal probability | $\bar{\pi}_j = P(\beta_j\ne0\mid \beta_{-j},\tilde{\beta}_j)$ | $O(1)$ |
| 3 | Sample the causal / non-causal state | $u \sim \mathcal{U}(0,1)$, then compare with $\bar{\pi}_j$ | $O(1)$ |
| 4.a | If non-causal, set effect to zero | $\beta_j \leftarrow 0$ | $O(1)$ |
| 4.b | If causal, sample from the Gaussian posterior | $\mu_j = \dfrac{\tilde{\beta}_j-a_j}{1 + \dfrac{Mp}{Nh^2}}$ <br> sample from $\mathcal{N}\left( \mu_j , v_j \right)$ | $O(1)$ |
| 6 | Save the current sample | update $\beta_j$ | $O(1)$ |

Therefore, for a single SNP $j$, the sampling complexity is:

$$
O(K).
$$


One Gibbs iteration sweeps all SNPs,  costing:

$$
O(MK).
$$

For $T$ Gibbs iterations, the total cost of LDPred is:

$$
O(TMK).
$$

Suming up with LD calculation, the overall complexity of `LDPred` is:
$$
O[(T+N)MK]
$$



# References
- Vilhjálmsson B, Yang J, Finucane H et al. Modeling Linkage Disequilibrium Increases Accuracy of Polygenic Risk Scores, *The American Journal of Human Genetics*, 97, 576-592
