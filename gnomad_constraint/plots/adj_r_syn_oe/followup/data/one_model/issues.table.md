| Problem (metric) | Current | Ridge rr | × Roulette | × Roulette + ridge | Plateau: weighted | Plateau: saturating | Plateau: logistic | Plateau: per substitution | 1M: groups + coverage | 1M: + MQ | 1M: + DNM regional + MQ | 1M: + gnomAD regional + MQ | 1M: full |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Synonymous level is off 1 (total obs/exp) | 1.015 | no | worse | worse | yes | yes | yes | yes | yes | yes | yes | yes | yes |
| Contexts are offset from each other (weighted SD of the 32 contexts' obs/exp) | 6.7% | no | partly | partly | no | no | worse | partly | partly | partly | partly | partly | partly |
| Low coverage is miscalibrated (obs/exp at AN 20–49% / 80–89%) | 1.34 / 0.95 | no | no | no | no | no | no | no | yes | yes | yes | yes | yes |
| Regional variation: per-gene noise (noise, core / no QC flag (1 = Poisson)) | 2.91 / 3.54 | partly | partly | partly | no | no | no | no | no | no | partly | partly | partly |
| Maternal hotspot regions run high (control genes' obs/exp, highest region (8p)) | 1.32 | partly | partly | yes | no | no | no | no | no | no | partly | partly | partly |
| Least constrained LOEUF decile runs high (obs/exp in decile 10 minus deciles 1–6) | +0.110 | partly | partly | partly | no | no | no | no | no | no | partly | partly | partly |
| Outliers with no QC flag (genes with \|z\| > 8) | 86 | partly | yes | yes | no | no | no | no | no | no | partly | partly | yes |
| Low-mapping-quality genes look depleted (flagged genes with z < −8) | 132 | no | no | no | no | no | no | no | no | yes, overshoots | yes, overshoots | yes, overshoots | yes, overshoots |
| Histone genes run high (outliers) | 24 | no | no | no | no | no | no | no | no | no | no | no | no |
| *Side effect:* pLoF expected vs synonymous | 0 | +0.7% | +3.6% | +4.4% | n/t | +0.9% | -1.9% | +3.6% | +3.0% | +3.0% | +3.1% | +4.3% | +7.1% |
| *Outcome:* NDD AUPRC change | 0 | +0.003 | -0.000 | +0.002 | n/t | +0.001 | -0.000 | -0.000 | +0.001 | +0.012 | +0.014 | +0.017 | +0.014 |
