library(GOsummaries)
load_all(package("GOsummaries"))

setwd("/Users/kolde/Raivo/Projects/GOsummaries/GOsummaries/visual_test/vtest/")

dirname = date()
system(sprintf("mkdir '%s'", dirname))

setwd(dirname)

data(gs_pca)
data(gs_limma)
data(gs_limma_exp)
data(gs_kmeans)

cat("\nPCA default")
plot(gs_pca, filename = "pca_default.pdf")

cat("\nPCA default, Tissue")
plot(gs_pca, classes = "Tissue", filename = "pca_tissue.pdf")

cat("\nPCA default, Tissue, fontsize = 10")
plot(gs_pca, classes = "Tissue", fontsize = 10, filename = "pca_tissue_fs10.pdf")

cat("\nPCA default, Tissue, fontsize = 15")
plot(gs_pca, classes = "Tissue", fontsize = 15, filename = "pca_tissue_fs15.pdf")

cat("\nPCA default, Tissue, fontsize = 25")
plot(gs_pca, classes = "Tissue", fontsize = 25, filename = "pca_tissue_fs25.pdf")

cat("\nPCA default, Tissue, panel_height = 0")
plot(gs_pca, classes = "Tissue", panel_height = 0, filename = "pca_tissue_ph0.pdf")

cat("\nLimma default")
plot(gs_limma, filename = "limma_default.pdf")

cat("\nLimma default, panel_height = 0")
plot(gs_limma, filename = "limma_default_ph0.pdf", panel_height = 0)

cat("\nLimma default, panel_height = 10")
plot(gs_limma, filename = "limma_default_10.pdf")

cat("\nLimma exp default")
plot(gs_limma_exp, classes = "Tissue", filename = "limma_exp_default.pdf")

cat("\nLimma exp violin")
plot(gs_limma_exp, classes = "Tissue", panel_plot = panel_violin, filename = "limma_exp_violin.pdf")

cat("\nLimma exp violin + box")
plot(gs_limma_exp, classes = "Tissue", panel_plot = panel_violin_box, filename = "limma_exp_violin_box.pdf")

cat("\nLimma exp violin + box, new color scheme")
cust = function(p, par){
	p = p + scale_fill_brewer(par$classes, type = "qual", palette = 1)
	return(p)
}
plot(gs_limma_exp, classes = "Tissue", panel_plot = panel_violin_box, panel_customize = cust, filename = "limma_exp_violin_box_colors.pdf")

cat("\nLimma exp box, new color scheme")
cust = function(p, par){
	p = p + scale_fill_brewer(par$classes, type = "qual", palette = 2)
	return(p)
}
plot(gs_limma_exp, classes = "Tissue", panel_plot = panel_boxplot, customize = cust, filename = "limma_exp_box_colors.pdf")

cat("\nKmeans exp default")
plot(gs_kmeans, classes = "Tissue", filename = "kmeans_exp_default.pdf")

cat("\nKmeans exp violin")
plot(gs_kmeans, classes = "Tissue", panel_plot = panel_violin, filename = "kmeans_exp_violin.pdf")

cat("\nKmeans exp violin + box")
plot(gs_kmeans, classes = "Tissue", panel_plot = panel_violin_box, filename = "kmeans_exp_violin_box.pdf")
















