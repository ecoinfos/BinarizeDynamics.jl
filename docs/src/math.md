# Mathematical Background

This document articulates the theoretical framework underlying `BinarizeDynamics.jl`. While the predecessor package, `BinarizeComp.jl`, focused on quantifying **Positional Heterogeneity** (1st-order statistics), `BinarizeDynamics.jl` expands this framework into **Structural Dynamics** (2nd-order interactions), quantifying how variations at different positions are mechanistically coupled.

This expansion is achieved by transitioning from scalar metrics (means and variances) to high-dimensional vector operations on binarized pairwise comparisons.

## 1. The Combinatorial Transform ($S \to \mathbf{B}$)

The fundamental operation in this framework is the transformation of categorical sequence data into a binary pairwise interaction space.

### 1.1 The Mapping Function

Let $S$ be a matrix of $n$ aligned sequences of length $L$, where $S_{u,p}$ represents the character at position $p$ in sequence $u$. The set of all unique unordered pairs of sequences is denoted as $\mathcal{P} = \{(u,v) \mid 1 \le u < v \le n\}$, with cardinality $M = \binom{n}{2}$.

For any given position $p$, we define a **Binary Difference Vector** $\mathbf{b}_p \in \{0, 1\}^M$. Each element of this vector, corresponding to a specific pair $k=(u,v)$, is determined by the binary distance metric $d$:

$$b_{p}^{(k)} = d(S_{u,p}, S_{v,p}) = \begin{cases} 
0 & \text{if } S_{u,p} = S_{v,p} \quad (\text{Consensus state}) \\
1 & \text{if } S_{u,p} \neq S_{v,p} \quad (\text{Divergent state})
\end{cases}$$

### 1.2 The Tensor Structure

The complete transformation results in a 3D **Binarized Tensor** $\mathbf{B}$ of dimensions $L \times M \times 1$.

*   **Dimension 1 ($L$)**: The aligned positions (columns of the original alignment).
*   **Dimension 2 ($M$)**: The pairwise comparisons (rows representing sequence pairs).
*   **Dimension 3 (1)**: The singleton dimension representing the binary state.

### 1.3 Connection to BinarizeComp.jl

This framework is a strict superset of the methodology presented in `BinarizeComp.jl`. In the previous work, the vector $\mathbf{b}_p$ was collapsed into a single scalar summary statistic, the **Mismatch Proportion ($P_d$)**:

$$P_d(p) = \mathbb{E}[\mathbf{b}_p] = \frac{1}{M} \sum_{k=1}^{M} b_{p}^{(k)}$$

Consequently, the **Positional Variance ($\sigma_d^2$)** derived in `BinarizeComp.jl` describes the variance of this vector:

$$\sigma_d^2(p) = \text{Var}(\mathbf{b}_p) = P_d(p)(1 - P_d(p))$$

`BinarizeDynamics.jl` preserves the full topology of the vector $\mathbf{b}_p$ rather than collapsing it. This allows us to calculate not just the magnitude of variation, but the *synchronicity* of variation between positions.

## 2. Quantifying Structural Coupling (The Phi Coefficient)

To measure the structural coupling between two positions $i$ and $j$, we analyze the statistical dependence between their respective binary vectors $\mathbf{b}_i$ and $\mathbf{b}_j$.

### 2.1 Contingency Table via Bitwise Logic

The relationship between vectors $\mathbf{b}_i$ and $\mathbf{b}_j$ can be summarized in a $2 \times 2$ contingency table. In Julia, these counts are computed efficiently using bitwise operations (`&`, `~`, `count_ones`) on `BitVector` structures, offering significant performance gains over standard integer arithmetic.

Let $N = M$ (the total number of pairs). The counts are defined as:

*   **$n_{11}$ (Co-divergence)**: Pairs that mismatch at both positions $i$ and $j$.
    $$n_{11} = \mathbf{b}_i \cdot \mathbf{b}_j \equiv \text{popcount}(\mathbf{b}_i \ \& \ \mathbf{b}_j)$$

*   **$n_{00}$ (Co-consensus)**: Pairs that match at both positions.
    $$n_{00} = (1-\mathbf{b}_i) \cdot (1-\mathbf{b}_j) \equiv \text{popcount}(\sim\mathbf{b}_i \ \& \ \sim\mathbf{b}_j)$$

*   **$n_{10}$ and $n_{01}$ (Discordance)**: Pairs that mismatch at one position but not the other.
    $$n_{10} = \text{popcount}(\mathbf{b}_i \ \& \ \sim\mathbf{b}_j), \quad n_{01} = \text{popcount}(\sim\mathbf{b}_i \ \& \ \mathbf{b}_j)$$

### 2.2 The Phi Coefficient Metric

We employ the **Phi Coefficient ($\phi$)** as the primary metric for structural coupling. It is the Pearson correlation coefficient for two binary variables:

$$\phi_{ij} = \frac{n_{11}n_{00} - n_{10}n_{01}}{\sqrt{n_{1\cdot}n_{0\cdot}n_{\cdot1}n_{\cdot0}}}$$

**Marginal Definitions:**

The marginal totals correspond directly to the mismatch ($P_d$) and match ($P_m$) counts derived in `BinarizeComp.jl`:

*   $n_{1\cdot} = n_{11} + n_{10} = M \times P_d(i)$
*   $n_{0\cdot} = n_{01} + n_{00} = M \times P_m(i)$

**Geometric Interpretation:**

The denominator represents the geometric mean of the variances of the two positions (scaled by $M^2$):

$$\text{Denominator} \propto \sqrt{\sigma_d^2(i) \times \sigma_d^2(j)}$$

**Dynamics Interpretation:**

*   $\phi_{ij} \approx +1$: **Perfect Positive Coupling**. A mismatch at $i$ guarantees a mismatch at $j$ (and consensus at $i$ guarantees consensus at $j$). The positions evolve as a single unit.
*   $\phi_{ij} \approx -1$: **Perfect Negative Coupling**. A mismatch at $i$ guarantees consensus at $j$. The positions are mutually exclusive in their variation.
*   $\phi_{ij} \approx 0$: **Independence**. The variation at $i$ provides no information about the state of $j$.

## 3. Signal Correction (APC)

In biological and high-dimensional categorical data, global noise often obscures specific functional couplings. This is typically driven by phylogeny or population structure, which causes high-entropy positions to show weak, spurious correlations with many other positions.

To isolate direct structural couplings, we apply **Average Product Correction (APC)**, a standard method in Direct Coupling Analysis (DCA) and protein structure prediction (e.g., AlphaFold's Evoformer logic).

### 3.1 The Background Model

We estimate the background signal shared between positions $i$ and $j$ based on their mean global connectivity:

$$\text{Background}_{ij} = \frac{\bar{\phi}_i \times \bar{\phi}_j}{\bar{\phi}_{total}}$$

Where:
*   $\bar{\phi}_i = \frac{1}{L} \sum_{k \neq i} \phi_{ik}$ is the average correlation of position $i$ with all other positions.
*   $\bar{\phi}_{total} = \frac{1}{L^2 - L} \sum_{i \neq j} \phi_{ij}$ is the global average correlation of the matrix.

### 3.2 The Corrected Metric ($\mathcal{F}$)

The final **Interaction Score $\mathcal{F}_{ij}$** is obtained by subtracting the additive background component:

$$\mathcal{F}_{ij} = \phi_{ij} - \text{Background}_{ij}$$

This corrected metric $\mathcal{F}$ serves as the input for all heatmap visualizations in `BinarizeDynamics.jl`.

## 4. Differential Dynamics (2.5D Analysis)

While 2D analysis reveals the static architecture of a system, **2.5D Analysis** quantifies how that architecture changes between distinct states or groups (e.g., Control vs. Experiment, Taxon A vs. Taxon B).

### 4.1 The Differential Algebra

Let $\mathbf{F}^{(A)}$ and $\mathbf{F}^{(B)}$ be the APC-corrected interaction matrices for Group A and Group B, respectively. The **Differential Map $\Delta$** is defined element-wise:

$$\Delta_{ij} = \mathcal{F}_{ij}^{(B)} - \mathcal{F}_{ij}^{(A)}$$

*   $\Delta_{ij} \gg 0$: Coupling strength has **increased** in Group B (Gain of Structure).
*   $\Delta_{ij} \ll 0$: Coupling strength has **decreased** in Group B (Loss of Structure).

### 4.2 Statistical Validation (Permutation Test)

To distinguish genuine structural shifts from stochastic fluctuations, we employ a non-parametric permutation test.

**Hypothesis ($H_0$)**: The assignment of sequences to Group A or Group B does not influence the structural coupling $\phi_{ij}$.

**Procedure:**

1.  Pool all $n_A + n_B$ sequences.
2.  Shuffle sequence labels randomly to form pseudo-groups $A'$ and $B'$ (preserving original sizes $n_A, n_B$).
3.  Recompute matrices and differential $\Delta'_{ij}$ for $N_{perm}$ iterations (default: 1000).
4.  Calculate $p$-value:

$$p_{ij} = \frac{1}{N_{perm}} \sum_{k=1}^{N_{perm}} \mathbb{I}(|\Delta'_{ij}| \ge |\Delta_{observed}|)$$

Only couplings where $p_{ij} < \alpha$ (e.g., 0.05) are considered statistically significant structural shifts.
