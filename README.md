# MOMF
Deconvolution analysis with the bulk RNA-seq data and single-cell RNA-seq data

## Installation
```R
#install devtools packages (devtools package)
install.packages("devtools")
#install MOMF package
devtools::install_github("sqsun/MOMF")
```

## Example
### Simulation data
The example data (`toy_example.Rdata`) is a simulated data.<br>
* `sc_counts`: scRNA-seq gene expression matrix (#cells x #genes); 
* `bulk_counts`: bulk RNA-seq gene expression matrix (#individuals x #genes). <br>

### Toy example
Two main functions `momf.fit` and `momf.computeRef` are used to do deconvoluation analysis.
```R
### load MOMF package
library(MOMF)

### load example data
> load("toy_example.RData")

### compute the cell type specific expression level as reference
> priorU <- momf.computeRef(sc_counts, sc_cell_type)

### create the gene list for MOMF 
> GList <- list(X1 = t(sc_counts), X2 = t(bulk_counts))

### run MOMF
> momf_res <- momf.fit(DataX = GList, DataPriorU=priorU, method="KL", rho=2, num_iter=3)

### output the cell type proportions
> cell_prop <- momf_res$cell_prop
> heatmap(cell_prop)
```

## Support
If that doesn't work, please reach out Sheng Yang (email:yangsheng@njmu.edu.cn).

