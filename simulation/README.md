# SingleCellAPA_evaluation

<p align="center">
<img src="https://github.com/LuChenLab/SCAPE/blob/main/image/simulation.jpg" alt="Simulation" width="600" />
</p>


## Package list

Here is the all package i used to perform performance testing, each package could be install following given steps.



|   Name    | Language |     Identification of  pA sites      |      PolyA_USE       | Invidual cell  matrix | Version |                          website                          |
| :-------: | :------: | :----------------------------------: | :------------------: | :-------------------: | :-----: | :-------------------------------------------------------: |
|  MAAPER   |    R     |     Modeling based on known PASs     |          No          |          No           |  1.1.1  |           https://github.com/Vivianstats/MAAPER           |
|   scAPA   | Shell, R | Homer and mclust::Mclust (R package) |          No          |          Yes          |  0.1.0  |             https://github.com/ElkonLab/scAPA             |
| scAPAtrap |    R     | derfinder::regionMatrix (R package)  | Yes; correct pA site |          Yes          |  0.1.0  |            https://github.com/BMILAB/scAPAtrap            |
| SCAPTURE  |  Shell   |                Homer                 |          No          |          Yes          |    1    |            https://github.com/YangLab/SCAPTURE            |
| scDaPars  |  python  |               DaPars2                |          No          |          No           |  0.1.0  |          https://github.com/YiPeng-Gao/scDaPars           |
|  Sierra   |    R     |    fit a Guassian with NLS or MLE    |          No          |          Yes          | 0.99.27 |              https://github.com/VCCRI/Sierra              |
|   SCAPE   |  python  |   Modeling based on insert size    |          Yes          |          Yes          |  1.0.0  |            https://github.com/LuChenLab/SCAPE             |
| polyApipe |  python  |   Soft clipped reads with poly(A)    |         Yes          |          Yes          |  0.1.0  | https://github.com/MonashBioinformaticsPlatform/polyApipe |



## How to use it



**Requirement**

1. Install all packages were listed above.
2. R (4.0.3)
3. python (3.7.3)
4. snakemake (6.6.1)



**Run**

```shell
snakemake -s simu_run.py
```



## Publication

Zhou *et al*. SCAPE: A mixture model revealing single-cell polyadenylation diversity and cellular dynamic during cell differentiation and reprogramming. *Nucleic Acids Research*.
