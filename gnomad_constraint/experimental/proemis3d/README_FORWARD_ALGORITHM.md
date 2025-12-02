<head>
  <style>
    .katex-display > .katex {
        text-align: left !important;
        margin: 0em 0em 1.5em 1.5em;
    }
  </style>
</head>

# Example of the forward algorithm

## Mock data (10 residues)

| Residue<br>(pos) | Expected variants<br>(exp) | Observed variants<br>(obs) |
| :--: | :--: | :--: |
|   1  | 0.9  |   0  |
|   2  | 1.0  |   0  |
|   3  | 1.1  |   0  |
|   4  | 0.8  |   0  |
|   5  | 1.2  |   1  |
|   6  | 1.0  |   1  |
|   7  | 0.9  |   1  |
|   8  | 1.1  |   1  |
|   9  | 1.0  |   1  |
|  10  | 0.8  |   1  |

## Compute the null NLL

### Fit the null parameter:
$$
\hat\gamma = \frac{\sum_i \text{obs}_i}{\sum_i \text{exp}_i}
$$

### Per-residue Poisson mean under the null:
$$
\lambda_i = \hat\gamma \cdot \text{exp}_i
$$

### Per-residue negative log-likelihood:
$$
\text{NLL}_i = \lambda_i - \text{obs}_i\log \lambda_i + \log(\text{obs}_i!) = \lambda_i - \text{obs}_i\log\lambda_i + \mathrm{lgamma}(\text{obs}_i{+}1)
$$

### Total null NLL:
$$
\text{NLL}_{\text{null}} = \sum_i \text{NLL}_i
$$


**per-position terms**: shows each $k_i=\text{obs}_i$, $e_i=\text{exp}_i$, $\lambda_i$, the per-site term, and a running cumulative sum.

### Summary of the Null-model fit:

$$
\sum \text{obs} = 6 \\
\sum \text{exp} = 9.8
$$
$$
\hat\gamma_{\text{null}} = \frac{6}{9.8} = 0.612245 \\
\text{NLL}_\text{null} \approx 8.99461
$$


## Round 1: evaluate the addition of candidate regions

The workflow is shown in detail for only the **best candidate region**, $R = [1,2,3,4]$, in full detail.

### Compute the observed and expected values for the candidate region:

$$
\text{obs}: [0,0,0,0] \\
\text{exp}: [0.9,1.0,1.1,0.8]
$$
$$
\sum\text{obs}_R = 0+0+0+0 = 0 \\
\sum\text{exp}_R = 0.9+1.0+1.1+0.8 = 3.8
$$


### Compute the observed and expected values for the remainining set of residues ${R^c} = [5,6,7,8,9,10]$:

$$
\text{obs}: [1,1,1,1,1,1] \\
\text{exp}: [1.2,1.0,0.9,1.1,1.0,0.8]
$$
$$
\sum\text{obs}*{R^c} = 6 \\
\sum\text{exp}*{R^c} = 1.2+1.0+0.9+1.1+1.0+0.8 = 6.0
$$

### MLEs for Poisson rates (one parameter per part)

$$
\hat\gamma_R =\frac{0}{3.8}=0.0
$$
$$
\hat\gamma_{R^c}=\frac{6}{6.0}=1.0
$$

### Per-position lambdas and NLL terms

#### Poisson NLL term for position $i$:
$$
\text{NLL}_i=\lambda_i - k_i\log\lambda_i + \gamma(k_i{+}1)
$$
$$
\lambda_i=\hat\gamma\cdot \text{exp}_i
$$

#### Region $[1,2,3,4]$ (using $\hat\gamma_R=0$):

All $\lambda_i=0$, all $k_i=0$ → each term $=0$

$$
\text{1: } \lambda=0, k=0 → term =0 \\
\text{2: } \lambda=0, k=0 → term =0 \\
\text{3: } \lambda=0, k=0 → term =0 \\
\text{4: } \lambda=0, k=0 → term =0 \\
$$
$$
\Rightarrow \text{NLL}_R = 0.000000
$$

#### Remaining [5,6,7,8,9,10] (using $\hat\gamma_{R^c}=1)$:

Here $\lambda_i = \text{exp}_i$ and $k_i=1$, so each term $=\lambda - \log\lambda$ (since $\gamma(2)=0$):

$$
\text{5: } term = 1.2 - \ln(1.2) = 1.017678444 \\
\text{6: } term = 1.0 - \ln(1.0) = 1.000000000 \\
\text{7: } term = 0.9 - \ln(0.9) = 1.005360516 \\
\text{8: } term = 1.1 - \ln(1.1) = 1.004689820 \\
\text{9: } term = 1.0 - 0 = 1.000000000 \\
\text{10: } term = 0.8 - \ln(0.8) = 1.023143551
$$

#### Sum:
$$
\text{NLL}_{R^c} = 1.017678444+1.0+1.005360516+1.004689820+1.0+1.023143551
= \mathbf{6.050872331}
$$

### Candidate total NLL and LRT

$$
\text{NLL}*\text{cand}=\text{NLL}*R+\text{NLL}*{R^c}=0+6.050872331=\mathbf{6.050872331} \\
\Delta \text{NLL}=\text{NLL}*\text{null}-\text{NLL}*\text{cand}=8.994610 - 6.050872331=\mathbf{2.943737669} \\
D = 2\Delta\text{NLL} = \mathbf{5.887475338}
$$
With **df = 1**,
$$
p*{\text{LRT}} = \Pr{\chi^2_1 \ge D} \approx \mathbf{0.015249}
$$

### Decision

$\alpha=0.001$

Compare $p=0.015249$ to $0.001$ → **not significant** at this stringent level → **do not accept** [1,2,3,4] in Round 1 (under this α).
