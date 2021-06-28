# MOMF
Deconvolution analysis with the bulk RNA-seq data and single-cell RNA-seq data

## Installation
```R
### install devtools packages (devtools package)
install.packages("devtools")

### install MOMF package
devtools::install_github("sqsun/MOMF")
```

## Example data
The example data (`toy_example.Rdata`) is a simulated data.<br>
* `sc_counts`: scRNA-seq gene expression matrix (#cells x #genes); 
* `sc_cell_type`: cell types for scRNA-seq data (#cells x 1);
* `bulk_counts`: bulk RNA-seq gene expression matrix (#individuals x #genes). <br>

## Example Code
Two main functions `momf.fit` and `momf.computeRef` are used to do deconvoluation analysis.
```R
### load MOMF package
> library(MOMF)

### load example data
> load("toy_example.RData")

### compute the cell type specific expression level as reference
> priorU <- momf.computeRef(sc_counts, sc_cell_type)

### create the gene list for MOMF 
> GList <- list(X1 = t(sc_counts), X2 = t(bulk_counts))

### run MOMF
> momf_res <- momf.fit(DataX = GList, DataPriorU=priorU, method="KL", rho=2, num_iter=100)

### output the cell type proportions
> cell_prop <- momf_res$cell.prop
> heatmap(cell_prop)
```
## Citation
Xifang Sun, Shiquan Sun, and Sheng Yang. *An efficient and flexible method for deconvoluting bulk RNAseq data
with single-cell RNAseq data*, 2019, DOI: 10.5281/zenodo.3373980. 

## Supports
Please reach out Xifang Sun or Sheng Yang (email: xfangsun@126.com or yangsheng@njmu.edu.cn) if you have any questions.
