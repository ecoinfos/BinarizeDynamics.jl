# Scientific Concepts

BinarizeDynamics.jl is built on a specific method for analyzing sequence relationships called **Pairwise Binarization**. This section explains the underlying logic and metrics.

## 1. Pairwise Binarization

Standard sequence analysis often deals with characters (A, C, G, T, etc.) directly. However, when analyzing the **dynamics** or differences between populations, it works better to convert the data into a binary "difference" map.

For a set of $N$ sequences of length $L$, we consider all unordered pairs of sequences $(i, j)$ where $1 \le i < j \le N$. There are $N(N-1)/2$ such pairs.

For each position $k$, we define a binary variable $X_{pair, k}$:

```math
X_{(i,j), k} = \begin{cases} 
1 & \text{if } Sequence_i[k] \neq Sequence_j[k] \\
0 & \text{if } Sequence_i[k] = Sequence_j[k]
\end{cases}
```

This transforms our data from `(N sequences x L positions)` to `(TotalPairs x L positions)`. The resulting matrix $X$ is stored efficiently as a `BitMatrix`.

## 2. Structural Analysis (Phi Coefficient)

To understand how positions in the sequence are structurally coupled, we calculate the correlation between the columns of the binary matrix $X$.

The **Phi coefficient** is the Pearson correlation coefficient applied to binary variables. For two positions $u$ and $v$::

```math
\phi_{uv} = \frac{n_{11}n_{00} - n_{10}n_{01}}{\sqrt{n_{1\cdot}n_{0\cdot}n_{\cdot 1}n_{\cdot 0}}}
```

Where:
* $n_{11}$ is the count of pairs where both positions show a difference ($X_{u}=1$ AND $X_{v}=1$).
* $n_{00}$ is the count of pairs where neither position shows a difference.
* The denominator is the geometric mean of the marginal counts.

A high positive $\phi$ indicates that when a mutation (difference) occurs at position $u$, it also tends to occur at position $v$ across the population.

## 3. Average Product Correction (APC)

In many biological datasets (like phylogenies), there is a strong background signal due to shared ancestry (phylogenetic bias). This global signal can obscure specific pairwise interactions.

**Average Product Correction (APC)** is a method to remove this background.

```math
APC_{uv} = \frac{\bar{\phi}_u \times \bar{\phi}_v}{\bar{\phi}_{total}}
```

Where $\bar{\phi}_u$ is the average correlation of position $u$ with all other positions.

The corrected metric is:
```math
\phi^{corrected}_{uv} = \phi_{uv} - APC_{uv}
```

BinarizeDynamics.jl computes this automatically when `analyze_structure(..., method=:phi, apc=true)` is called.
