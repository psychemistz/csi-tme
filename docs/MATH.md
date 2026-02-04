# From HSIC Fundamentals to SPLISOSM: A Complete Mathematical Framework for Spatial Isoform Analysis

## Executive Summary

This document provides a ground-up explanation of the Hilbert-Schmidt Independence Criterion (HSIC) and its application in SPLISOSM for detecting spatially variable isoform patterns. We begin with the fundamental statistical problem of testing independence between variables, progress through kernel methods and reproducing kernel Hilbert spaces, and culminate in the complete SPLISOSM framework with worked numerical examples throughout.

**Application Focus**: This framework is designed for analyzing **cytokines and secreted proteins** in spatial transcriptomics, where alternative splicing serves as a regulatory layer complementing classical proteolytic activation mechanisms.

---

## Part I: The Fundamental Problem of Independence Testing

### 1.1 Why Testing Independence Matters

Consider a spatial transcriptomics experiment where we measure isoform expression across tissue locations. The core biological question is: **Do isoform ratios depend on spatial location?** If isoform A is preferentially used in the tumor core while isoform B dominates at the invasive margin, there exists a statistical dependence between isoform usage and spatial position.

Formally, we have two random variables:
- **X**: Isoform composition (e.g., proportion of each isoform)
- **Y**: Spatial coordinates (e.g., x,y position on tissue slide)

The null hypothesis is independence: P(X,Y) = P(X)P(Y), meaning knowing the spatial location tells us nothing about expected isoform ratios.

### 1.2 Classical Approaches and Their Limitations

**Pearson correlation** tests linear dependence between scalar variables:

$$r = \frac{\sum_{i=1}^{n}(x_i - \bar{x})(y_i - \bar{y})}{\sqrt{\sum(x_i-\bar{x})^2}\sqrt{\sum(y_i-\bar{y})^2}}$$

*Limitation*: Pearson correlation equals zero for many dependent relationships. Consider X ~ Uniform(-1,1) and Y = X². These are clearly dependent (Y is completely determined by X), yet Corr(X,Y) = 0 because the relationship is nonlinear.

**Spearman correlation** handles monotonic nonlinearity but still fails for non-monotonic relationships.

**Chi-squared tests** require discretization, losing information and requiring arbitrary bin choices.

**The fundamental challenge**: We need a statistic that captures *any* form of dependence, including nonlinear and multivariate relationships, without making parametric assumptions about the dependence structure.

---

## Part II: Kernels and Feature Spaces — A Complete Introduction

Kernels are the mathematical foundation underlying HSIC and all kernel-based statistical methods. Understanding kernels deeply is essential for appreciating why HSIC can detect arbitrary forms of dependence and why SPLISOSM's choice of spatial kernel matters for statistical power.

### 2.1 The Fundamental Problem: Limitations of Linear Methods

Many statistical methods assume linear relationships between variables. Linear methods fail when relationships are nonlinear. The natural solution is to transform the data into a space where nonlinear relationships become linear. If we could find a transformation φ such that φ(X) and φ(Y) are linearly related, we could apply standard linear methods to the transformed data. This is the core idea behind kernel methods.

### 2.2 Feature Maps and Feature Spaces

A **feature map** is a function φ: X → H that transforms data from the original input space X into a (potentially higher-dimensional) feature space H. The key insight is that nonlinear patterns in X may become linear patterns in H.

**Concrete Example: Polynomial Features**

Consider data points on a line that cannot be separated by a single threshold:
```
Class A: x = -2, -1, 1, 2
Class B: x = -0.5, 0, 0.5
```

These are interleaved and linearly inseparable. Now apply the transformation φ(x) = (x, x²):
```
Class A: φ(x) = (-2, 4), (-1, 1), (1, 1), (2, 4)
Class B: φ(x) = (-0.5, 0.25), (0, 0), (0.5, 0.25)
```

In this two-dimensional feature space, a horizontal line at x₂ = 0.5 perfectly separates the classes.

**The Dimensionality Challenge**

For complex nonlinear relationships, we may need very high-dimensional (or even infinite-dimensional) feature spaces. Consider all polynomial features up to degree d for input in ℝⁿ:

$$\phi(x) = (1, x_1, x_2, \ldots, x_n, x_1^2, x_1x_2, \ldots, x_n^d)$$

The dimension of this feature space is O(nᵈ), which quickly becomes intractable. For d = 10 and n = 100, the feature space has over 10²⁰ dimensions.

### 2.3 The Kernel Trick: Computing Inner Products Implicitly

The **kernel trick** is the observation that many algorithms only require inner products between data points, not the explicit feature representations. If we can compute the inner product ⟨φ(x), φ(x')⟩ directly without computing φ(x) and φ(x') individually, we avoid the dimensionality explosion.

**Definition (Kernel Function)**: A kernel is a function k: X × X → ℝ that computes the inner product between feature representations:

$$k(x, x') = \langle \phi(x), \phi(x') \rangle_{\mathcal{H}}$$

### 2.4 Common Kernel Functions

**Linear Kernel**:
$$k(x, x') = x^T x'$$

The feature map is the identity: φ(x) = x.

**Polynomial Kernel**:
$$k(x, x') = (x^T x' + c)^d$$

The feature space contains all monomials up to degree d.

**Gaussian (Radial Basis Function / RBF) Kernel**:
$$k(x, x') = \exp\left(-\frac{\|x - x'\|^2}{2\sigma^2}\right)$$

This is the most important kernel for HSIC. It corresponds to an **infinite-dimensional** feature space. The bandwidth parameter σ controls the "reach" of similarity.

**Properties of the Gaussian Kernel**:
- k(x, x) = 1 for all x (self-similarity is maximal)
- k(x, x') → 0 as ||x - x'|| → ∞ (distant points are dissimilar)
- k(x, x') > 0 for all x, x' (all points have some positive similarity)
- The kernel is **universal** (see Section 2.7)

**Laplacian Kernel**:
$$k(x, x') = \exp\left(-\frac{\|x - x'\|}{\sigma}\right)$$

Similar to Gaussian but uses L1 norm instead of squared L2 norm. Less smooth, which can be advantageous for detecting sharp boundaries.

### 2.5 Positive Definite Kernels and Mercer's Theorem

**Definition (Positive Definite Kernel)**: A function k: X × X → ℝ is a positive definite kernel if for any finite set of points {x₁, ..., xₙ} ⊂ X, the **Gram matrix** K with entries K_{ij} = k(x_i, x_j) is positive semi-definite:

$$\sum_{i,j} c_i c_j k(x_i, x_j) \geq 0 \quad \text{for all } c_1, \ldots, c_n \in \mathbb{R}$$

**Mercer's Theorem**: Every positive definite kernel corresponds to an inner product in some feature space.

### 2.6 Reproducing Kernel Hilbert Spaces (RKHS)

For every positive definite kernel k, there exists a unique Hilbert space H_k called the **reproducing kernel Hilbert space** (RKHS) with the following properties:

1. **Feature map definition**: φ(x) = k(x, ·) ∈ H_k
2. **Reproducing property**: f(x) = ⟨f, k(x, ·)⟩_H for all f ∈ H_k

Setting f = k(x', ·) yields:
$$k(x', x) = \langle k(x', \cdot), k(x, \cdot) \rangle_{\mathcal{H}}$$

This confirms that the kernel computes inner products between feature representations.

**The Representer Theorem**: For optimization problems over an RKHS:
$$\min_{f \in \mathcal{H}} \sum_{i=1}^n L(y_i, f(x_i)) + \lambda \|f\|_{\mathcal{H}}^2$$

the minimizer has the form:
$$f^*(x) = \sum_{i=1}^n \alpha_i k(x_i, x)$$

This makes kernel methods computationally tractable even in infinite-dimensional spaces.

### 2.7 Universal Kernels and Characteristic Kernels

**Definition (Universal Kernel)**: A kernel k on a compact domain X is universal if its RKHS H is dense in C(X), the space of continuous functions on X.

**Definition (Characteristic Kernel)**: A kernel k is characteristic if the mean embedding μ_P := E_{x~P}[φ(x)] uniquely determines the probability distribution P.

**Theorem**: The Gaussian kernel is both universal and characteristic on any compact subset of ℝⁿ.

**Why This Matters for HSIC**: When using characteristic kernels, the cross-covariance operator C_xy is zero if and only if X and Y are independent. This is the theoretical foundation for HSIC's ability to detect any form of dependence.

### 2.8 Kernel Matrices and Their Properties

In practice, we work with the **kernel matrix** (or Gram matrix) K ∈ ℝⁿˣⁿ defined by K_{ij} = k(x_i, x_j).

**Properties of Kernel Matrices**:
1. **Symmetry**: K = Kᵀ
2. **Positive semi-definiteness**: All eigenvalues λᵢ ≥ 0
3. **Eigendecomposition**: K = UΛUᵀ

**The Centered Kernel Matrix**:
$$\tilde{K} = HKH \quad \text{where} \quad H = I - \frac{1}{n}\mathbf{1}\mathbf{1}^T$$

---

## Part III: The Hilbert-Schmidt Independence Criterion

This section presents the rigorous theoretical foundations of HSIC as established by Gretton et al. (2007).

### 3.1 The Formal Independence Testing Problem

**Problem (Independence Testing)**: Let P_xy be a Borel probability measure defined on a domain X × Y, and let P_x and P_y be the respective marginal distributions. Given an i.i.d. sample Z := (X, Y) = {(x₁, y₁), ..., (xₘ, yₘ)}, does P_xy factorize as P_x P_y?

**Statistical Hypothesis Framework**:
- **Null hypothesis H₀**: P_xy = P_x P_y (independence)
- **Alternative hypothesis H₁**: P_xy ≠ P_x P_y (dependence)

### 3.2 Cross-Covariance Operators in RKHS

Let F be an RKHS with kernel k and feature map φ(x), and G be an RKHS on Y with kernel l and feature map ψ(y).

**Definition (Cross-Covariance Operator)**: The cross-covariance operator C_xy : G → F is defined such that:

$$\langle f, C_{xy}g \rangle_{\mathcal{F}} = \mathbb{E}_{xy}\left[(f(x) - \mathbb{E}_x[f(x)])(g(y) - \mathbb{E}_y[g(y)])\right]$$

Explicitly:
$$C_{xy} := \mathbb{E}_{xy}\left[(\phi(x) - \mu_x) \otimes (\psi(y) - \mu_y)\right]$$

**Key Theorem**: When F and G are **universal** RKHSs, the largest singular value of C_xy equals zero if and only if x ⊥⊥ y.

### 3.3 HSIC Definition and Population Expression

**Definition (Hilbert-Schmidt Independence Criterion)**:

$$\text{HSIC}(P_{xy}, \mathcal{F}, \mathcal{G}) := \|C_{xy}\|_{HS}^2$$

The population expression:

$$\text{HSIC} = \mathbb{E}_{xx'yy'}[k(x,x')l(y,y')] + \mathbb{E}_{xx'}[k(x,x')]\mathbb{E}_{yy'}[l(y,y')] - 2\mathbb{E}_{xy}\left[\mathbb{E}_{x'}[k(x,x')]\mathbb{E}_{y'}[l(y,y')]\right]$$

**Interpretation**:
1. **First term**: Expected similarity in X-space times similarity in Y-space under the joint distribution
2. **Second term**: Product of expected similarities under marginals (what we'd expect under independence)
3. **Third term**: Correction term for centering

When X and Y are independent, HSIC = 0.

### 3.4 Empirical Estimators

**The Trace Formula (V-statistic)**:

$$\widehat{\text{HSIC}}_b(Z) = \frac{1}{m^2}\text{tr}(KHLH)$$

where K is the m×m kernel matrix for X, L is the m×m kernel matrix for Y, and H = I - (1/m)11^T is the centering matrix.

### 3.5 Asymptotic Distribution Theory

**Theorem (Distribution under H₁ - Dependence)**:
$$m^{1/2}\left(\widehat{\text{HSIC}}_b(Z) - \text{HSIC}(P_{xy})\right) \xrightarrow{D} \mathcal{N}(0, \sigma_u^2)$$

**Theorem (Distribution under H₀ - Independence)**:
$$m \cdot \widehat{\text{HSIC}}_b(Z) \xrightarrow{D} \sum_{l=1}^{\infty} \lambda_l z_l^2$$

where z_l ~ N(0,1) and {λ_l} are eigenvalues of the integral operator.

**Key Insight**: The null distribution depends on the eigenvalues of the centered kernel matrices. This is why SPLISOSM's choice of spatial kernel fundamentally affects test behavior.

### 3.6 Gamma Approximation for P-value Computation

The infinite sum of chi-squared variables can be approximated with a two-parameter Gamma distribution:

$$m \cdot \widehat{\text{HSIC}}_b(Z) \approx \text{Gamma}(\alpha, \beta)$$

where:
$$\alpha = \frac{(\mathbb{E}[\widehat{\text{HSIC}}_b])^2}{\text{Var}(\widehat{\text{HSIC}}_b)}, \quad \beta = \frac{m \cdot \text{Var}(\widehat{\text{HSIC}}_b)}{\mathbb{E}[\widehat{\text{HSIC}}_b]}$$

---

## Part IV: SPLISOSM's Spatial Kernel Innovation and Theoretical Contributions

SPLISOSM introduces two major theoretical innovations: (1) a proof that all low-rank kernel approximations sacrifice statistical power (Theorem 1), and (2) a mathematically sound approach to handle undefined ratios in compositional data (Theorem 2).

### 4.1 From Correlation to Dependence Testing

SPLISOSM's empirical HSIC estimator:

$$\widehat{\text{HSIC}}(X, Y) = \frac{1}{(n-1)^2}\text{tr}(KHLH)$$

The test statistic follows an asymptotic χ² mixture under the null:

$$T_{\text{HSIC}} := \frac{1}{n}\text{tr}(KHLH) \xrightarrow{d}_{n\to\infty} \sum_{i,j=1}^{n} \lambda_i \mu_j z_{ij}^2$$

where λ_i and μ_j are eigenvalues of centered kernels HKH and HLH, and z_ij are i.i.d. standard Gaussian variables.

### 4.2 Theorem 1: Low-Rank Kernel Approximations Sacrifice Statistical Power

**Lemma 1 (Poincaré Separation / Cauchy Interlacing Theorem)**: Let A be an n × n symmetric matrix with eigenvalues λ₁ ≤ λ₂ ≤ ··· ≤ λₙ. For any n × m matrix P with orthonormal columns where m < n, let B = PᵀAP with eigenvalues μ₁ ≤ ··· ≤ μ_m. Then:

$$\lambda_i \leq \mu_i \leq \lambda_{n-m+i} \quad \text{for all } i \in \{1, 2, \ldots, m\}$$

**Theorem 1 (Limited Power of Low-Rank Kernel Tests)**:

Consider a one-dimensional target variable X ∈ ℝⁿ. Let K_sp be an n × n spatial kernel matrix of rank d < n - 1. Then:

1. There exists a non-constant spatial pattern X₀ for which I(X₀; K_sp) has **no detection power**
2. Such undetectable patterns form a subspace of dimension **at least n - d**

**Proof**: The test statistic with linear kernel K_X can be expressed as:

$$T(X, K_{sp}) = X^T (H K_{sp} H) X = X^T \tilde{K}_{sp} X$$

Using eigendecomposition K̃_sp = UΛUᵀ:
$$T(X, K_{sp}) = (U^T X)^T \Lambda (U^T X)$$

For T(X, K_sp) = 0, construct X such that UᵀX only has non-zero elements in positions corresponding to zero eigenvalues. ∎

**Corollary (Orthogonal Decomposition)**: Any spatial pattern X can be decomposed into:
$$X = X_\perp + X_\parallel$$
such that T(X_⊥, K_sp) = 0.

**Implications for SPARK-X**: SPARK-X uses a low-rank spatial kernel of rank ≤ 22. For typical datasets with n ≫ 22, the space of undetectable spatial patterns is **much larger** than the space of detectable patterns.

### 4.3 The Smoother ICAR Spatial Kernel

**Definition 4.1 (Smoother ICAR Kernel)**: Let W be the n × n adjacency matrix of a spatial graph where W_ij = 1 if i and j are mutual k-nearest neighbors. The ICAR kernel is:

$$K = (D - \rho W)^{-1}$$

where D = diag{Σⱼ W_ij} is the degree matrix and ρ ∈ [0, 1) controls spatial autocorrelation.

**Properties**:
1. K is positive definite and defines a multivariate Gaussian distribution
2. K captures topological properties and is invariant to scaling, translation, and rotation
3. K⁻¹ is sparse, enabling efficient eigendecomposition
4. The eigenvectors represent meaningful spatial frequency components

**Proposition 4.1 (Connection to Spectral Graph Theory)**: Let L = D - W be the graph Laplacian. The normalized kernel K̃ and normalized Laplacian L̃ share the same eigenvectors, forming the basis for the Graph Fourier Transform.

For eigenvector v of L̃ with eigenvalue λ:
$$\tilde{K}v = \frac{1}{1-\rho+\rho\lambda}v$$

**Rayleigh Quotient Interpretation**: The smoothness of pattern X:
$$R_{\tilde{L}}(X) = \frac{\sum_{i \sim j}(x_i - x_j)^2 (d_i d_j)^{-1/2}}{\sum_i x_i^2}$$

Higher values indicate higher frequency (rapid spatial changes).

### 4.4 Low-Rank GFT Kernel Approximation

**Definition 4.2 (Low-rank GFT Kernel)**: The rank-r GFT kernel is:

$$K_{sp}^r = \sum_{i=1}^{r} \gamma_i u_i u_i^T := U_r \Gamma_r U_r^T$$

where uᵢ is the i-th lowest-frequency spatial component.

The test statistic becomes:
$$T(X, K_{sp}^r) = \sum_{i=1}^{r} \gamma_i \hat{X}_i^2$$

**Key Insight**: As rank increases, the undetectable space shifts toward high-frequency patterns (typically noise).

### 4.5 Theorem 2: Handling Missing Data in Compositional Analysis

**Challenges**: ~50% of spots in Visium have zero UMI counts for a median gene, leading to undefined isoform ratios.

**Definition 4.3 (Zero-Padded Centered Kernel)**:

$$K_{X_n} = \begin{pmatrix} H_m K_{X_m} H_m & 0 \\ 0 & 0 \end{pmatrix}_{n \times n}$$

This effectively replaces undefined ratios with the mean in feature space.

**Theorem 2 (Mean-Replacement Preserves Test Validity)**:

Let T_m be the test statistic for m spots with defined ratios, T_n for all n spots. Then:

1. **Statistic relationship**: T_m = (m/n)T_n
2. **Null threshold bounds**: Can be computed using only the spectrum of the global kernel L_n

### 4.6 SPLISOSM's Three Test Statistics

**HSIC-GC (Gene Counts)**: Tests whether total gene expression varies spatially.
$$T_{\text{GC}} = \frac{1}{n}\text{tr}(K_g H L_Y H), \quad K_g = gg^T$$

**HSIC-IR (Isoform Ratios)**: Tests whether isoform usage proportions vary spatially—the key test for **spatially variable processing (SVP)**.
$$T_{\text{IR}} = \frac{1}{n}\text{tr}(K_X H L_Y H)$$

**HSIC-IC (Isoform Counts)**: Detects either isoform usage or expression changes.

**Key Design Choice**: SPLISOSM uses a linear kernel without log-ratio transformation for compositional data.

### 4.7 Conditional HSIC Test for Spatial Confounding Control

**Problem**: X (isoform usage) and Y (biological covariate) may share spatial patterns simply due to mutual dependence on spatial coordinates.

**Solution**: Learn regression functions f_X: S → X and f_Y: S → Y, then test independence of residuals.

For non-linear relationships:
$$X - \hat{X} = \lambda(K_S + \lambda I)^{-1}X := R_S X$$

**Conditional Kernel**:
$$K_{X|S} := R_S K_X R_S$$

---

## Part V: Comparison with SPARK-X

### 5.1 Statistical Philosophy Differences

| Aspect | SPARK-X | SPLISOSM |
|--------|---------|----------|
| **Framework** | Variance component test | Kernel independence test |
| **Null hypothesis** | Spatial variance = 0 | X ⊥⊥ Y (independence) |
| **Statistic** | Score statistic T = Σᵢ λᵢχᵢ² | HSIC = (1/n)tr(KHLH) |
| **Compositionality** | Not designed for ratios | Native support |

### 5.2 Kernel Construction Differences

| Aspect | SPARK-X | SPLISOSM |
|--------|---------|----------|
| **Rank** | 22 (fixed, low-rank) | Full-rank (n-1 after centering) |
| **Spatial encoding** | Distance-based projections | Graph-based ICAR |
| **Pattern selection** | 11 predefined transformations | Data-adaptive frequency weighting |
| **Compositional data** | Not supported | ZPC kernel with Theorem 2 guarantees |
| **Tissue topology** | Ignored (Euclidean distance) | Respected (neighborhood graph) |

### 5.3 Empirical Power Comparison

From SPLISOSM simulations (n = 1000 spots, 3 isoforms, spatial gradient):

| Effect Size | SPARK-X (adapted) | SPLISOSM HSIC-IR |
|-------------|-------------------|------------------|
| 0.1 (weak) | 12% | **34%** |
| 0.2 | 31% | **67%** |
| 0.3 | 58% | **89%** |
| 0.5 (strong) | 82% | **98%** |

The power advantage is most pronounced at weak effect sizes—precisely the regime relevant for detecting subtle cytokine isoform gradients.

### 5.4 When to Use Each Method

**Use SPARK-X when**:
- Detecting spatially variable gene expression (not isoforms)
- Dataset has >100,000 spots
- Only interested in large-scale patterns

**Use SPLISOSM when**:
- Analyzing isoform ratios or compositional data
- Detecting localized patterns (e.g., tumor boundaries)
- Statistical power is critical
- Need spatial confounding control

---

## Part VI: P-value Computation via Liu's Method

Under the null hypothesis, HSIC-based test statistics follow:
$$T \xrightarrow{d} \sum_i w_i \chi_i^2(1)$$

**Liu's method** approximates this with a scaled chi-squared:
$$T \approx c \cdot \chi^2(d) + \mu$$

Parameters matched via moments:
- μ = Σ w_i
- σ² = 2Σ w_i²
- κ = 8(Σ w_i³)/(Σ w_i²)^{3/2}
- d = 12/κ²
- c = σ/√(2d)

---

## References

**Core HSIC and Kernel Methods**:
1. Gretton A, Bousquet O, Smola A, Schölkopf B. Measuring statistical dependence with Hilbert-Schmidt norms. *ALT 2005*. Lecture Notes in Computer Science, vol 3734. doi: 10.1007/11564089_7

2. Gretton A, Fukumizu K, Teo CH, Song L, Schölkopf B, Smola AJ. A kernel statistical test of independence. *Advances in Neural Information Processing Systems 20 (NIPS 2007)*. MIT Press; 2008:585-592.

3. Zhang K, Peters J, Janzing D, Schölkopf B. Kernel-based conditional independence test and application in causal discovery. *UAI 2011*:804-813.

4. Liu H, Tang Y, Zhang HH. A new chi-square approximation to the distribution of non-negative definite quadratic forms in non-central normal variables. *Computational Statistics & Data Analysis* 2009;53:853-856.

**Spatial Transcriptomics Methods**:
5. Zhu J, Sun S, Zhou X. SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. *Genome Biology* 2021;22:184. doi: 10.1186/s13059-021-02404-0

6. Sun S, Zhu J, Zhou X. Statistical analysis of spatial expression pattern for spatially resolved transcriptomic studies. *Nature Methods* 2020;17:193-200.

7. Svensson V, Teichmann SA, Stegle O. SpatialDE: identification of spatially variable genes. *Nature Methods* 2018;15:343-346.

**SPLISOSM Publication**:
8. Su J, Qu Y, Schertzer M, Yang H, Jiang J, Lhakhang T, Nelson TM, Park S, Lai Q, Fu X, Choi SW, Knowles DA, Rabadan R. Mapping isoforms and regulatory mechanisms from spatial transcriptomics data with SPLISOSM. *Nature Biotechnology* 2026. doi: 10.1038/s41587-025-02965-6

**Spatial Graph and ICAR Models**:
9. Besag J. Spatial interaction and the statistical analysis of lattice systems. *Journal of the Royal Statistical Society B* 1974;36:192-236.

10. Su J, Schertzer M, Bhavna I, Mukherjee S, Bhatt P, Rabadan R. Smoother: a unified and modular framework for incorporating structural dependency in spatial omics data. *Genome Biology* 2023;24:291.

11. Shuman DI, Narang SK, Frossard P, Ortega A, Vandergheynst P. The emerging field of signal processing on graphs. *IEEE Signal Processing Magazine* 2013;30:83-98. doi: 10.1109/MSP.2012.2235192

12. Chung FRK. *Spectral Graph Theory*. CBMS Regional Conference Series in Mathematics, No. 92. American Mathematical Society; 1997. ISBN: 978-0-8218-0315-8

---

**Software Resources**:
- **SPLISOSM Python Package**: https://github.com/JiayuSuPKU/SPLISOSM
- **Documentation**: https://splisosm.readthedocs.io/
- **Paper Reproducibility Scripts**: https://github.com/JiayuSuPKU/SPLISOSM_paper
- **Processed Data (Zenodo)**: doi: 10.5281/zenodo.16905935

---

*Document Version: 2.0 | Last Updated: February 2026*
*Covers: Mathematical foundations (HSIC/Kernels), SPLISOSM methodology, Theorem 1 (low-rank limitations), Theorem 2 (compositional data), Comparison with SPARK-X*
