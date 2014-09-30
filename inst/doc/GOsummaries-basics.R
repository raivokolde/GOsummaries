### R code from vignette source 'GOsummaries-basics.Rnw'

###################################################
### code chunk number 1: style-Sweave
###################################################
BiocStyle::latex()


###################################################
### code chunk number 2: myCodeBlock
###################################################
library(GOsummaries, quietly=TRUE)
library(vegan, quietly=TRUE)


###################################################
### code chunk number 3: example1
###################################################
# Define gene lists
genes1 = c("203485_at", "209469_at", "209470_s_at", "203999_at", 
	"205358_at", "203130_s_at", "210222_s_at", "202508_s_at", "203001_s_at", 
	"207957_s_at", "203540_at", "203000_at", "219619_at","221805_at",
	 "214046_at", "213135_at", "203889_at", "209990_s_at", "210016_at", 
	"202507_s_at","209839_at", "204953_at", "209167_at", "209685_s_at",  
	"211276_at", "202391_at", "205591_at","201313_at")
genes2 = c("201890_at", "202503_s_at", "204170_s_at", "201291_s_at", 
	"202589_at", "218499_at", "209773_s_at", "204026_s_at", "216237_s_at", 
	"202546_at", "218883_s_at", "204285_s_at", "208659_at", "201292_at", 
	"201664_at")
gl = list(List = list(genes1, genes2)) # Two lists per component

# Construct gosummaries objects
gs = gosummaries(gl)

plot(gs, fontsize = 8, filename = "figure2.pdf")


###################################################
### code chunk number 4: Example2
###################################################
data(tissue_example)

# Filter genes and perform k-means
sd = apply(tissue_example$exp, 1, sd)
exp2 = tissue_example$exp[sd > 0.75,]
exp2 = exp2 - apply(exp2, 1, mean)
kmr = kmeans(exp2, centers = 6, iter.max = 100)

# Create gosummaries object
exp2[1:6, 1:5]
head(tissue_example$annot)

gs_kmeans = gosummaries(kmr, components = 1:2, exp = exp2, annotation = tissue_example$annot)
plot(gs_kmeans, fontsize = 8, classes = "Tissue", filename = "figure3.pdf")


###################################################
### code chunk number 5: Example3
###################################################
cust = function(p, par){
	p = p + scale_fill_brewer(par$classes, type = "qual", palette = 2)
	return(p)
}
plot(gs_kmeans, panel_plot = panel_violin, panel_customize = cust, 
classes = "Tissue", components = 1:2, filename = "ex3.pdf")


###################################################
### code chunk number 6: ExampleUserSupplied
###################################################
wcd1 = data.frame(Term = c("KLF1", "KLF2", "POU5F1"), Score = c(0.05, 0.001, 0.0001))
wcd2 = data.frame(Term = c("CD8", "CD248", "CCL5"), Score = c(0.02, 0.005, 0.00001))


###################################################
### code chunk number 7: ExampleUserSupplied2
###################################################
gs = gosummaries(wc_data = list(Results1 = wcd1, Results2 = wcd2))
plot(gs, filename = "figure5.pdf")


###################################################
### code chunk number 8: ExampleUserSupplied3
###################################################
# To get two word clouds per block use neted lists
gs = gosummaries(wc_data = list(Results = list(wcd1, wcd2)))
plot(gs, filename = "figure6.pdf")


###################################################
### code chunk number 9: ExampleMetagenomic
###################################################
data(metagenomic_example)

# Run Principal Coordinate Analysis on Bray-Curtis dissimilarity matrix 
pcoa = cmdscale(vegdist(t(metagenomic_example$otu), "bray"), k = 3)

# By turning off the GO analysis we can show the names of taxa
gs = gosummaries(pcoa, metagenomic_example$otu, metagenomic_example$annot,
                 show_genes = T, gconvert_target = NULL, n_genes = 30)

plot(gs, class = "BodySite", fontsize = 8, file = "figure7.pdf")


###################################################
### code chunk number 10: SessionInfo
###################################################
sessionInfo()


