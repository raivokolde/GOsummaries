GOsummaries
===========

An R package that visualizes the GO enrichment results as word clouds and arranges them together with figures of experimental data. This allows us to draw informative summary plots for analyses such as differential expression or clustering, where for each gene list we display its behaviour in the experiment alongside with its GO annotations. The approach is especially interesting for Principal Component Analysis (PCA), where we can annotate the principal axes functionally based on the weights of the genes.

## Installation
The package is available ar CRAN, so the installation can be done as usual. 
```s
install.packages("GOsummaries")`
```

More comprehensive user guide can be found in the [vignette](http://cran.r-project.org/web/packages/GOsummaries/vignettes/GOsummaries-basics.pdf).

## Examples

### PCA 
```s
# Perform PCA on samples
pcr = prcomp(t(tissue_example$exp))

# Create gosummaries object
gs_pca = gosummaries(pcr, annotation = tissue_example$annot)

# Plot
plot(gs_pca, classes = "Tissue")
```
![pca](http://raivokolde.github.com/GOsummaries/images/pca.png)

### K-means clustering
```s
library(GOsummaries)

data(tissue_example)

# Filter genes and perform k-means
sd = apply(tissue_example$exp, 1, sd)
exp2 = tissue_example$exp[sd > 0.75,]
exp2 = exp2 - apply(exp2, 1, mean)
kmr = kmeans(exp2, centers = 6, iter.max = 100)

# Create gosummaries object
gs_kmeans = gosummaries(kmr, exp = exp2, annotation = tissue_example$annot)
plot(gs_kmeans, components = 1:2)
```
![kmeans](http://raivokolde.github.com/GOsummaries/images/kmeans.png)
