## Dummydata object
dummydata = function(gl, max){
	res = list( 
		mat = data.frame(x = factor(1:length(gl)), y = unlist(lapply(gl, length))),
		max = max
	)
	
	class(res) = c(class(res), "dummyData")
	
	return(res)
}

is.dummyData = function(x) inherits(x, "dummyData")

add_dummydata.gosummaries = function(gosummaries){
	max = max(unlist(lapply(gosummaries, function(x) lapply(x$Gene_lists, length))))
	for(i in seq_along(gosummaries)){
		gosummaries[[i]]$Data = dummydata(gosummaries[[i]]$Gene_lists, max)
	}
	
	return(gosummaries)
}

##

## gosummaries object constructor and related functions
 
#' Constructor for gosummaries object
#' 
#' Constructor for gosummaries object that contains all the necessary information to d
#' raw the figure, like gene lists and their annotations, expression data and all the 
#' relevant texts.
#' 
#' The object is a list of "components", with each component defined by a gene list or a 
#' pair of gene lists. Each "component" has the slots as follows:
#' \itemize{
#'   \item \bold{Title}: title string of the component. (Default: the names of the gene 
#' lists)
#'   \item \bold{Gene_lists}: list of one or two gene lists with names gl1 (and gl2 if 
#' present).
#'   \item \bold{GPR}: g:Profiler results based on the Gene_lists slot. 
#'   \item \bold{Data}: the related data (expression values, PCA rotation, ...) that is 
#' used to draw the "panel" i.e. theplot above the wordclouds. In principle there is no 
#' limitation what  kind of data is there as far as the function that is provided to 
#' draw that in \code{\link{plot.gosummaries}} can use it.
#'   \item \bold{Percentage}: a text that is drawn on the right top corner of every 
#' component. In case of PCA this is the percentage of variation the component explains, 
#' by default it just depicts the number of genes in the Gene_lists slot.
#' }
#' 
#' The GO enrichment analysis is performed using g:Profiler web toolkit and its 
#' associated R package \code{gProfileR}. This means the computer has to have internet 
#' access to annotate the gene lists. Since g:Profiler can accept a wide range of gene 
#' IDs then user usually does not have to worry about converitng the gene IDs into right 
#' format. To be absolutely sure the tool recognizes the gene IDs one can check if they 
#' will give any results in \url{http://biit.cs.ut.ee/gprofiler/gconvert.cgi}. 
#' 
#' There can be a lot of results for a typical GO enrichment analysis but usually these 
#' tend to be pretty redundant. Since one can fit only a small number of categories into 
#' a word cloud we have to bring down the number of categories to show an reduce the 
#' redundancy. For this we apply an algortihm that selects from groups of related 
#' categories only the ones with the best p-values. The algorithm works as follows. 
#' \itemize{
#'   \item Delete all the results where category size is smaller than \code{min_size}, 
#' larger than \code{max_size}, enrichment p-value is larger than 
#' \code{p_value_threshold} or the category does not belong to one of the 
#' \code{go_branches}.
#'   \item Map all the remaining results back to the GO graph and find connected 
#' components.
#'   \item From each connected component retain the category with smallest p-value. 
#'   \item If more than \code{max_signif} categories are still present return 
#' \code{max_signif} ones with the best p-values
#' } 
#' 
#' The slots of the object can be filled with custom information using a function 
#' \code{\link{add_to_slot.gosummaries}}. 
#' 
#' By default the Data slot is filled with a dataset that contains the number of genes 
#' in the Gene_lists slot. Expression data can be added to the object for example by 
#' using function \code{\link{add_expression.gosummaries}}. It is possible to derive 
#' your own format for the Data slot as well, as long as a panel plotting function for 
#' this data is alaso provided (See \code{\link{panel_boxplot}} for further 
#' information).
#' 
#' There are several constructors of gosummaries object that work on common analysis 
#' result objects, such as \code{\link{gosummaries.kmeans}}, 
#' \code{\link{gosummaries.MArrayLM}} and \code{\link{gosummaries.prcomp}} corresponding 
#' to k-means, limma and PCA results.
#'
#' @param x list of arrays of gene names (or list of lists of arrays of gene names)
#' @param organism the organism that the gene lists correspond to. The format should be 
#' as follows: "hsapiens", "mmusculus", "scerevisiae", etc.
#' @param \dots Description of parameter 
#' @return   A gosummaries type of object
#' 
#' @seealso \code{\link{gosummaries.kmeans}}, 
#' \code{\link{gosummaries.MArrayLM}}, \code{\link{gosummaries.prcomp}}
#' 
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  # Example1
#' 
#' @rdname gosummaries
#' @export
gosummaries = function(x, ...){
	UseMethod("gosummaries", x)
}

#' @param go_branches GO tree branches and pathway databases as denoted in g:Profiler (Possible values: BP, CC, MF, ke, re) 
#' @param p_value_threshold threshold for p-values that are corrected for multiple testing
#' @param min_size minimal size of functional category to be considered
#' @param max_size maximal size of functional category to be considered
#' @param max_signif maximal number of categories returned per query
#' @param ordered_query logical showing if the lists are ordered or not (it determines if the ordered query algorithm is used in g:Profiler)
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  # Example1
#' 
#' @rdname gosummaries
#' @method gosummaries default
#' @S3method gosummaries default
#' @export
gosummaries.default = function(x, organism = "hsapiens", go_branches = c("BP", "ke", "re"), p_value_threshold = 1e-2, min_size = 50, max_size = 1000, max_signif = 40, ordered_query = T, ...){
	components = 1:length(x)
	
	# Find the number of lists per component (assumed to be the same over components)
	if(is.list(x[[1]])){
		k = 2
	}
	else{
		k = 1
	}
	
	# Create the resulting data structure
	res = list()
	
	for(i in components){
		comp = list(
			Title = names(x)[i],
			Gene_lists = NULL,
			GPR = NULL,
			Data = NULL,
			Percentage = NULL,
			Organism = organism
		)
		
		if(k == 1){
			comp$Gene_lists = list(gl1 = x[[i]])
			comp$GPR = list(gpr1 = NULL)
			comp$Percentage = sprintf("n: %d", length(x[[i]]))
		}
		
		if(k == 2){
			comp$Gene_lists = list(gl1 = x[[i]][[1]], gl2 = x[[i]][[2]])
			comp$GPR = list(gpr1 = NULL, gpr2 = NULL)
			comp$Percentage = sprintf("Up: %d\nDown: %d", length(x[[i]][[2]]), length(x[[i]][[1]]))
		}
		
		res[[as.character(components[i])]] = comp
	}
	
	class(res) = "gosummaries"
	
	res = add_dummydata.gosummaries(res)
	res = annotate.gosummaries(res, organism = organism, go_branches = go_branches, p_value_threshold = p_value_threshold, min_size = min_size, max_size = max_size, max_signif = max_signif, ordered_query = ordered_query)
	
	return(res)
} 

 
#' @param gosummaries a gosummaries object
#' @rdname add_to_slot.gosummaries
#' @export
is.gosummaries = function(x) inherits(x, "gosummaries")

#' @param \dots not used
#' 
#' @rdname add_to_slot.gosummaries
#' @method print gosummaries
#' @S3method print gosummaries
#' @export
print.gosummaries = function(x, ...){
	for(a in x){
		cat(sprintf("Component title: %s\n", a$Title))
		cat("\n")
		cat(sprintf("%s\n", "Head of gene lists"))
		for(l in a$Gene_lists){
			print(head(l))
		}
		cat("\n")
		cat(sprintf("%s\n", "Top annotation results"))
		for(l in a$GPR){
			print(head(l))
		}
		cat("\n")
		cat(sprintf("%s\n", "Data"))
		print(head(a$Data))
		
		cat("\n")
		cat(sprintf("%s\n", "Percentage slot:"))
		cat(a$Percentage)
		
		cat("\n===================================================\n\n")
	}
}

 
#' Functions to interact with gosummaries object
#'
#' @param x a gosummaries object
#' @param slot the component slot name to be filled (e.g Title, Percentage, etc.)
#' @param values list of values where each element corresponds to one component
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  #Example1
#' 
#' @rdname add_to_slot.gosummaries
#' @export
add_to_slot.gosummaries = function(x, slot, values){
	if(length(x) != length(values)) stop("Length of gosummaries object and values does not match")
	
	for(i in seq_along(x)){
		x[[i]][[slot]] = values[[i]]
	}
	
	return(x)
}
##

## Annotate gosummaries with g:Profiler
filter_gprofiler_default = function(gpr, go_branches = c("BP", "KEGG", "REAC"), p_value_threshold, min_size, max_size, max_signif){
	gpr = gpr[gpr$domain %in% go_branches,]
	
	if(nrow(gpr) == 0) return(gpr)
	
	gpr$sg = 0
	gpr$sg[which(gpr$term.size > max_size | gpr$relative.depth == 1)]= 1
	gpr$sg2 = cumsum(gpr$sg)
	gpr = gpr[gpr$term.size <= max_size,]
	gpr = ddply(gpr, "sg2", function(x) x[which.min(x$p.value),, drop = F])
	
	gpr = ddply(gpr, "query.number", function(x){
		x[order(x$p.value), ][1:min(max_signif, nrow(x)),]
	})
	
	return(gpr)
}

filter_gprofiler_min_pvalue = function(gpr, go_branches = c("BP", "KEGG", "REAC"), p_value_threshold, min_size, max_size, max_signif){
	gpr = gpr[gpr$domain %in% go_branches,]
	
	gpr = gpr[gpr$term.size > min_size & gpr$term.size < max_size & gpr$p.value < p_value_threshold, ]
	
	gpr = ddply(gpr, "query.number", function(x) {
		rank = rank(x$p.value)
		return(x[rank <= max_signif, ])
	})
	
	return(gpr)
}

filter_gprofiler = function(gpr, go_branches = c("BP", "KEGG", "REAC"), p_value_threshold, min_size, max_size, max_signif, heuristic = "default"){
	res = switch(heuristic, 
		default=filter_gprofiler_default(gpr, go_branches = go_branches, p_value_threshold = p_value_threshold, min_size = min_size, max_size = max_size, max_signif = max_signif), 
		min_pvalue=filter_gprofiler_min_pvalue(gpr, go_branches = go_branches, p_value_threshold = p_value_threshold, min_size = min_size, max_size = max_size, max_signif = max_signif))
	
	return(res)
}

adjust_gprofiler_output = function(gpr){
	for(i in c("query.number", "significant", "p.value", "term.size", "query.size", "overlap.size", "precision", "recall", "subgraph.number", "relative.depth")){
		gpr[[i]] = as.numeric(as.character(gpr[[i]]))
	}
	
	return(gpr)
}

annotate.gosummaries = function(gosummaries, organism, components = 1:length(gosummaries), go_branches, p_value_threshold, min_size, max_size, max_signif, ordered_query, go_pruning_heuristic = "default"){
	
	if(!is.gosummaries(gosummaries)) stop("Function requires a gosummaries type of  object")
	
	components = as.character(components)
	
	#Compile gene lists 
	gl = NULL
	for(i in seq_along(components)){
		for(j in seq_along(gosummaries[[components[i]]]$Gene_lists)){
			l = gosummaries[[components[i]]]$Gene_lists[[j]]
			if(length(l) == 0){
				l = c("uuuuuuu1", "3uuuuuuuuuuu5")
			}
			gl = c(gl, list(l))
		}
	}
	
	# Run g:Profiler analysis 
	gpr = gprofiler(query = gl, organism = organism, ordered_query = ordered_query)
	
	# Clean and filter the results
	gpr$query.number = as.numeric(as.character(gpr$query.number))
	gpr = adjust_gprofiler_output(gpr)
	gpr = filter_gprofiler(gpr, go_branches = go_branches, p_value_threshold = p_value_threshold, min_size = min_size, max_size = max_size, max_signif = max_signif, heuristic = go_pruning_heuristic)
	gpr = gpr[, c("query.number", "term.id", "p.value", "term.name")]
	colnames(gpr) = c("query.number", "Term.id", "P.value", "Term.name")

	k = 1
	for(i in seq_along(components)){
		for(j in seq_along(gosummaries[[components[i]]]$GPR)){
			lname = paste("gpr", j, sep = "")
			gosummaries[[components[i]]]$GPR[[lname]] = gpr[gpr$query.number == k, -1] 
			k = k + 1
		}
	}
	
	return(gosummaries)
}

##

## Fill Data slot in gosummaries object
padz = function(x, n=max(nchar(x))) gsub(" ", "0", formatC(x, width=n))

 
#' Add expression data to gosummaries object
#' 
#' Function to add expression data and its annotations to a gosummaries object.  
#' 
#' The data is added to the object in a "long" format so it would be directly usable by 
#' the ggplot2 based panel drawing functions \code{\link{panel_boxplot}} etc. For each 
#' component it produces a data frame with columns:
#' \itemize{
#'   \item \bold{x} : sample IDs for the x axis, their factor order is the same as on the columns of input matrix 
#'   \item \bold{y} : expression values from the matrix
#'   \item  . . . : sample annotation columns from the annotation table that can be displayed on figure as colors.
#' } 
#'
#' @param gosummaries a gosummaries object
#' @param exp an expression matrix, with row names corresponding to the names in the 
#' Gene_lists slot
#' @param annotation a \code{data.frame} describing the samples, its row names should match 
#' with column names of \code{exp}
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  # Example1
#' 
#' @export
add_expression.gosummaries = function(gosummaries, exp, annotation){
	if(!is.gosummaries(gosummaries)) stop("Function requires a gosummaries type object")
	
	if(!all(colnames(exp) %in% rownames(annotation)) & !is.null(annotation)) 
		stop("Column names of expression matrix and row names of annotation dataframe do not match") 
	
	for(i in seq_along(gosummaries)){
		a = list()
		for(j in seq_along(gosummaries[[i]]$Gene_lists)){
			e = data.frame(exp[gosummaries[[i]]$Gene_lists[[j]], ])
			d = melt(e)
			colnames(d) = c("ID", "y")
			d$x = paste(j, padz(match(d$ID, colnames(e))), d$ID, sep = "")
			a[[j]] = d
		}
		a = do.call("rbind", a)
		if(!is.null(annotation)){
			annotation$ID = make.names(rownames(annotation))
			a = merge(a, annotation)
			a = a[, !(colnames(a) %in% "ID")]
			a$x = factor(a$x, levels = sort(unique(as.character(a$x))))
		}
		if(length(gosummaries[[i]]$Gene_lists) == 1){
			class(a) = c(class(a), "oneListExpData")
		}
		if(length(gosummaries[[i]]$Gene_lists) == 2){
			class(a) = c(class(a), "twoListExpData")
		}
		gosummaries[[i]]$Data = a
	}
	
	return(gosummaries)
}

add_pca.gosummaries = function(gosummaries, pcr, annotation){
	if(!is.gosummaries(gosummaries)) stop("Function requires a gosummaries type object")
	
	if(!all(rownames(pcr$x) %in% rownames(annotation)) & !is.null(annotation)) 
		stop("Column names of expression matrix and row names of annotation dataframe do not match")
	
	for(i in seq_along(gosummaries)){
		a = data.frame(x = pcr$x[, i])
		
		if(!is.null(annotation)){
			a$ID = rownames(pcr$x)
			annotation$ID = rownames(annotation)
			a = merge(a, annotation)
			a = a[, !(colnames(a) %in% "ID")]
		}
		
		class(a) = c(class(a), "pcaData")
		gosummaries[[i]]$Data = a
	}
	
	return(gosummaries)
}

##

## Adjust Wordcloud appearance parameters in gosummaries object
shorten_strings = function(strings, max){
	strings = as.character(strings)
	n = nchar(strings)
	strings[n > max] = paste(substr(strings[n > max], 1, max - 3), "...", sep = "")
	return(strings)
}

adjust_wordcloud_appearance = function(gosummaries, term_length = 35, wordcloud_colors = c("grey70", "grey10")){
	if(all(is.null(ldply(gosummaries, function(x){ldply(x$GPR, function(y) data.frame(y$P.value))})[,2]))){
		return(gosummaries)
	}
	
	
	best_pval = -log10(min(ldply(gosummaries, function(x){ldply(x$GPR, function(y) data.frame(y$P.value))})[,2]))
	
	for(i in seq_along(gosummaries)){
		# Shorten GO names
		for(j in seq_along(gosummaries[[i]]$GPR)){
			gosummaries[[i]]$GPR[[j]]$Term.name = shorten_strings(gosummaries[[i]]$GPR[[j]]$Term.name, term_length)
			gosummaries[[i]]$GPR[[j]]$Colors = colorRampPalette(wordcloud_colors)(100)[ceiling(-log10(gosummaries[[i]]$GPR[[j]]$P.value) / best_pval * 100)]
		}
	}
	
	return(gosummaries)
}
##

## GGplot utility functions (uses some functions from unreleased gtable package, if released should be moved to that)
is.gtable = function(x) {
  inherits(x, "gtable")
}

gtable_trim = function(x) {
  stopifnot(is.gtable(x))
  
  w = range(x$layout$l, x$layout$r)
  h = range(x$layout$t, x$layout$b)
  
  x$widths = x$widths[seq.int(w[1], w[2])]
  x$heights = x$heights[seq.int(h[1], h[2])]
  
  x$layout$l = x$layout$l - w[1] + 1
  x$layout$r = x$layout$r - w[1] + 1
  x$layout$t = x$layout$t - h[1] + 1
  x$layout$b = x$layout$b - h[1] + 1
  
  x
}

gtable_filter = function(x, pattern, fixed = FALSE, trim = TRUE) {
  
  matches = grepl(pattern, x$layout$name, fixed = fixed)
  x$layout = x$layout[matches, , drop = FALSE]
  x$grobs = x$grobs[matches]
  
  if (trim) x = gtable_trim(x)
  
  return(x)
}


## Define zeroGrob
zeroGrob = function(){
	grob(cl = "zeroGrob", name = "NULL")
}
widthDetails.zeroGrob = 
heightDetails.zeroGrob = 
grobWidth.zeroGrob = 
grobHeight.zeroGrob = function(x) unit(0, "cm")
drawDetails.zeroGrob = function(x, recording) {}
is.zero <- function(x) identical(x, zeroGrob())
##

panelize_ggplot2 = function(plot_function, customize_function, par){
	res = function(data, fontsize, legend = F){
		p = ggplot_gtable(ggplot_build(customize_function(plot_function(data, fontsize, par), par)))
		
		if(legend){
			index = grepl("guide-box", p$layout$name, fixed = FALSE)
			if(sum(index) == 0) return(zeroGrob())
			a = p$grobs[[which(index)]]
			width = p$widths[p$layout$l[index]]
			# height = unit(length(getGrob(a, gPath = "key", grep = TRUE, global = T)) / 2 * 1.4 * fontsize, "points") + unit(1, "lines")
			height = a$heights[2] - unit(3, "mm")
			fg = frameGrob()
			fg = packGrob(fg, a, width = width, height = height)
			return(fg)
		}
		else{
			return(gtable_filter(p, "panel"))
		}
	}
	
	return(res)
} 
##

## Plotting utility functions
vplayout = function(x, y){
	return(viewport(layout.pos.row = x, layout.pos.col = y))
}

open_file_con = function(filename, width, height){
	# Get file type
	r = regexpr("\\.[a-zA-Z]*$", filename)
	if(r == -1) stop("Improper filename")
	ending = substr(filename, r + 1, r + attr(r, "match.length"))

	f = switch(ending,
		pdf = function(x, ...) pdf(x, ...),
		png = function(x, ...) png(x, units = "in", res = 300, ...),
		jpeg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
		jpg = function(x, ...) jpeg(x, units = "in", res = 300, ...),
		tiff = function(x, ...) tiff(x, units = "in", res = 300, compression = "lzw", ...),
		bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
		stop("File type should be: pdf, png, bmp, jpg, tiff")
	)
	
	f(filename, width = width, height = height)
}
##

## Raw plotting functions
calc_component_dimensions = function(component, par, legend){
	# Wordcloud height
	nr = max(laply(component$GPR, nrow))
	if(length(component$GPR) > 1){
		wc_height = ifelse(nr > 3, max(nr / 4.5, 3), nr)
		arrows_height = 1.5
	}
	else{
		wc_height = ifelse(nr > 3, max(nr / 8, 3), nr)
		arrows_height = 0.5
	}
	
	# Compile results
	res = list(
		title_height = unit(1.5, "lines"),
		panel_height = unit(par$panel_height, "lines"),
		arrows_height = unit(arrows_height, "lines"),
		wc_height = unit(wc_height, "lines"),
		
		panel_width = unit(par$panel_width, "points"),
		percentage_width = NULL,
		wc_width = unit(par$panel_width / length(component$GPR), "points")
	)
	
	return(res)
}

calc_components_dimensions = function(gosummaries, par, legend){
	component_dims = list()
	component_heights = unit(0, "npc") 
	
	# Calculate percentage slot width
	n = max(laply(gosummaries, function(x) nchar(max(unlist(strsplit(x$Percentage, "\n"))))))
	if(n > 0){
		percentage_width = unit(1.3, "grobwidth", textGrob(paste(rep("k", max(3, n)), collapse = ""), gp = gpar(fontsize = par$fontsize)))
	}
	else{
		percentage_width = unit(0.2, "cm")
	}
	
	for(i in seq_along(gosummaries)){
		dims = calc_component_dimensions(gosummaries[[i]], par, legend)
		dims$percentage_width = percentage_width
		
		component_heights = unit.c(component_heights,  dims$title_height + dims$panel_height + dims$arrows_height + dims$wc_height + unit(1, "lines"))
		
		component_dims[[i]] = dims
	}
	
	# Calculate percentage slot width
	
	component_heights = component_heights[-1]
	
	component_width = component_dims[[1]]$panel_width + percentage_width
	
	return(list(component_dims = component_dims, component_heights = component_heights, component_width = component_width))
}

gen_legend = function(legend_data, par){
	n = length(legend_data$colors)
	# Create Grobs
	title = textGrob(legend_data$title, y = 1, x = 0, just = c(0, 1), gp = gpar(fontface = "bold", cex = 0.8))
	
	rect_height = unit(1.4 * par$fontsize, "points")
	yy = unit(1, "npc") - unit(0.8, "lines") - (0:(n - 1)) * rect_height
	boxes = rectGrob(x = 0, y = yy, height = rect_height, width = rect_height, just = c(0, 1), gp = gpar(col = 0, fill = legend_data$colors))
	
	yyy = yy - rect_height * 0.65
	gl = gList()
	length = c(rep(0, n), convertWidth(grobWidth(title), "in"))
	for(i in 1:n){
		gl[[i]] = textGrob(legend_data$labels[[i]], x = rect_height + unit(3, "points"), y = yyy[i], hjust = 0, gp = gpar(cex = 0.8))
		length[i] = convertWidth(grobWidth(gl[[i]]), "in")
	}
	
	# Calculate size
	height = unit(0.8, "lines") + n * rect_height
	width = rect_height + unit(3, "points") + unit(max(length), "in")
	 
	# Put together a frame
	fg = frameGrob()
	fg = packGrob(fg, rectGrob(width = width, height = height, gp = gpar(col = NA)))
	fg = packGrob(fg, title)
	fg = packGrob(fg, boxes)
	for(i in 1:length(gl)){
		fg = packGrob(fg, gl[[i]])
	}
	
	return(fg)
}

gen_wordcloud_legend = function(gosummaries, par){
	legend_data = list()
	
	legend_data$title = "Enrichment P-value"
	legend_data$colors = colorRampPalette(rev(par$wordcloud_colors))(3)
	
	# Calculate p-value breakpoints
	best_pval = -log10(min(ldply(gosummaries, function(x){ldply(x$GPR, function(y) data.frame(y$P.value))})[,2]))	
	
	breaks = grid.pretty(c(0, best_pval))
	if(length(breaks) %% 2 == 0){
		kesk = mean(c(breaks[length(breaks) / 2], breaks[length(breaks) / 2 + 1]))
	}
	else{
		kesk = breaks[ceiling(length(breaks) / 2)]
	}
	
	legend_data$labels = c(substitute(10 ^ -p, list(p = breaks[length(breaks)])), substitute(10 ^ -p, list(p = kesk)), 1)
	
	return(gen_legend(legend_data, par))
}

plot_wordcloud = function(words, freq, color, algorithm, dimensions){
	if(length(words) > 0){
		return(plotWordcloud(words, freq, colors = color, random.order = F, min.freq = -Inf, rot.per = 0, scale = 0.85, max_min = c(1, 0), algorithm = algorithm, add = F, grob = T, dimensions = dimensions))
	}
	return(zeroGrob())
}

plot_arrow = function(end, par){
	x = switch(end, first = c(0, 0.95), both = c(0.05, 0.95), last = c(0.05, 1))
	res = linesGrob(x = x, y = 0.5, arrow = arrow(ends = end, type = "closed", angle = 15, length = unit(0.1, "inches")), gp = gpar(lwd = 0.3 * par$fontsize, col = "grey40"))
	return(res)
}

plot_component = function(data_component, plot_panel, par, component_dims){
	
	# Create gtable
	heights = with(component_dims, unit.c(title_height, panel_height, arrows_height + wc_height))
	widths = with(component_dims, unit.c(panel_width, percentage_width))
	gtable_component = gtable::gtable(widths, heights)
	
	
	# Add title
	title = textGrob(x = 0, hjust = 0, data_component$Title, gp = gpar(fontface = "bold", fontsize = par$fontsize))
	gtable_component = gtable::gtable_add_grob(gtable_component, title, 1, 1)
	
	# Add plot
	if(par$panel_height != 0){
			p = plot_panel(data_component$Data, par$fontsize)
			b = rectGrob(gp = gpar(lwd = 1.5, col = "grey40"))
			gtable_component = gtable::gtable_add_grob(gtable_component, gList(p, b), 2, 1)
		}
	
	# Add percentage
	p = textGrob(data_component$Percentage, x = 0.1, y = 1, vjust = 1, hjust = 0, gp = gpar(fontsize = par$fontsize * 0.8))
	gtable_component = gtable::gtable_add_grob(gtable_component, p, 2, 2)
	
	# Arrows and wordclouds
	heights = with(component_dims, unit.c(arrows_height, wc_height))
	widths = with(component_dims, rep(wc_width, length(data_component$GPR)))
	gtable_aw = gtable::gtable(widths, heights)
	
	if(length(data_component$GPR) == 1){
		wc = plot_wordcloud(data_component$GPR$gpr1$Term.name, -log10(data_component$GPR$gpr1$P.value), data_component$GPR$gpr1$Colors, algorithm = "leftside", dimensions = with(component_dims, unit.c(wc_width, wc_height)))
		gtable_aw = gtable::gtable_add_grob(gtable_aw, wc, 2, 1)
	}
	
	if(length(data_component$GPR) == 2){
		gtable_aw = gtable::gtable_add_grob(gtable_aw, plot_arrow(end = "first", par), 1, 1)
		gtable_aw = gtable::gtable_add_grob(gtable_aw, plot_arrow(end = "last", par), 1, 2)
		
		wc1 = plot_wordcloud(data_component$GPR$gpr1$Term.name, -log10(data_component$GPR$gpr1$P.value), data_component$GPR$gpr1$Colors, algorithm = "leftside", dimensions = with(component_dims, unit.c(wc_width, wc_height)))
		wc2 = plot_wordcloud(data_component$GPR$gpr2$Term.name, -log10(data_component$GPR$gpr2$P.value), data_component$GPR$gpr2$Colors, algorithm = "rightside", dimensions = with(component_dims, unit.c(wc_width, wc_height)))
		
		gtable_aw = gtable::gtable_add_grob(gtable_aw, wc1, 2, 1)
		gtable_aw = gtable::gtable_add_grob(gtable_aw, wc2, 2, 2)
	}
	
		
	gtable_component = gtable::gtable_add_grob(gtable_component, gtable_aw, 3, 1)
	gtable_component = gtable::gtable_add_padding(gtable_component, unit(0.5, "lines"))
	
	return(gtable_component)

}

plot_motor = function(gosummaries, plot_panel, legend = T, par = list(fontsize = 10, panel_height = 5, panel_width = 405), filename = NA){
	
	# Calculate dimensions for the picture components
	component_dimensions = calc_components_dimensions(gosummaries, par, legend)
	
	# Create legends and their dimensions 
	if(par$panel_height != 0){
		panel_legend = plot_panel(gosummaries[[1]]$Data, par$fontsize, legend = T)
	}
	else{
		panel_legend = zeroGrob()
	}
	
	wordcloud_legend = gen_wordcloud_legend(gosummaries, par)
	
	legend_width = max(grobWidth(panel_legend), grobWidth(wordcloud_legend))
	
	
	
	
	# Open connection to file if filename specified
	if(!is.na(filename)){
		width = convertWidth((component_dimensions$component_width + legend_width)  , "inches", valueOnly = T) * 1.05
		height = convertHeight(sum(component_dimensions$component_heights), "inches", valueOnly = T) * 1.05
		open_file_con(filename, width, height)
	}
	
	
	
	
	
	
	# Define layout parameters
	grid.newpage()
	gp = list(fontsize = par$fontsize)
	pushViewport(viewport(gp = do.call(gpar, gp)))
	
	# Create viewport for padding 
	pushViewport(viewport(width = (component_dimensions$component_width + legend_width)  * 1.05, height = sum(component_dimensions$component_heights) * 1.05))
	
	# Create layout of the picture
	pushViewport(viewport(layout = grid.layout(nrow = length(gosummaries), ncol = 2, widths = unit.c(component_dimensions$component_width, legend_width), heights = component_dimensions$component_heights)))
	
	# Draw components
	for(i in seq_along(gosummaries)){
		pushViewport(vplayout(i, 1))
		plot_component(gosummaries[[i]], plot_panel, par, component_dimensions$component_dims[[i]])
		popViewport()
	}
	
	# Draw legends
	hpl = grobHeight(panel_legend)
	
	pushViewport(vplayout(seq_along(gosummaries), 2))
	pushViewport(viewport(x = unit(0, "npc"), y = unit(1, "npc") - unit(1.5, "lines"), height = hpl, width = grobWidth(panel_legend),  just = c(0, 1)))
	grid.draw(panel_legend)
	popViewport()
	
	pushViewport(viewport(x = unit(0.5, "lines") + unit(2, "mm"), y = unit(1, "npc") - unit(ifelse(as.numeric(convertUnit(hpl, "cm")) > 0, 2.5, 1.5), "lines") - hpl, height = grobHeight(wordcloud_legend), width = grobWidth(wordcloud_legend), just = c(0, 1)))
	grid.draw(wordcloud_legend)
	popViewport()
	popViewport()
	
	# Pop layout viewport
	popViewport()
	
	# Pop padding viewport
	popViewport()
	
	# Close connection if filename specified
	if(!is.na(filename)){
		dev.off()
	}
}
##

## Panel functions for expression data
panel_crossbar = function(data, fontsize = 10, par){
	qq = function(x, ...) {
		res = data.frame(ymin = quantile(x, 0.05), y = median(x), ymax = quantile(x, 0.95))
		
		return(res)
	}
	
	if(!is.null(par$classes)){
		p = qplot(x = x, y = y, geom = "crossbar", stat = "summary", width = 0.7, fun.data = qq, fill = data[, par$classes], data = data) + geom_bar(aes(x = 1, y = 0, fill = data[, par$classes]), stat = "identity")
	}
	else{
		p = qplot(x = x, y = y, geom = "crossbar", stat = "summary", width = 0.7, fun.data = qq, data = data)
	}
	
	if(inherits(data, "twoListExpData")){
		p = p + geom_vline(xintercept = length(unique(data$x))/2 + 0.5)
	}
	
	p = p + theme_bw(base_size = fontsize)
	
	return(p)
}

 
#' Panel drawing functions 
#' 
#' These functions are used to draw the panel portion of every component based on the 
#' Data slots in gosummaries object. These concrete functions assume the data is 
#' presented as does \code{\link{add_expression.gosummaries}}. They provide three 
#' options: boxplot, violin plot (which shows the distrubution more precisely) and both 
#' combined.
#' 
#' These functions specify in principle the general setting for the panels, like which 
#' "geom"-s, how the data is transformed and summarized, etc. To make small adjustments 
#' to the figure such as changing color scheme, write your own customization function 
#' (See \code{\link{customize}} as example).
#' 
#' It is possible to write your own panel plotting function, as long as the parameters  
#' used and the return value are similar to what is specified here. When writing a new 
#' panel function one only has to make sure that it matches the data given in the Data 
#' slot of the gosummaries object.
#' 
#' @param data the data from Data slot of the go sumaries object
#' @param fontsize fontsize in points, mainly used to ensure that the legend fontsizes 
#' match
#' @param par additional parameters for drawing the plot, given as list. These functions 
#' use only \code{classes} slot for figuring out which parameter to use for coloring the 
#' "geom"-s. However, when building a custom function it provides a way to give extra 
#' parameters to the plotting function. 
#' @return  It returns a function that can draw a ggplot2 plot of the data in Data slot 
#' of a gosummaries object. The legend and the actual plots for the panels are extracted 
#' later from the figure produced by this function.
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  # Example1
#' 
#' @rdname panel_boxplot
#' @export
panel_boxplot = function(data, fontsize = 10, par){
	qq = function(x) {
		bs = boxplot.stats(x)$stats 
		data.frame(ymin = bs[1], lower = bs[2], middle = bs[3], upper = bs[4], ymax = bs[5])
	
	}
	if(!is.null(par$classes)){
			p = qplot(x = x, y = y, geom = "boxplot", stat = "summary", fun.data = qq, fill = data[, par$classes], width = 0.7, data = data) 
	}
	else{
			p = qplot(x = x, y = y, geom = "boxplot", stat = "summary", fun.data = qq, width = 0.7, data = data) 
	}
	
	if(inherits(data, "twoListExpData")){
		p = p + geom_vline(xintercept = length(unique(data$x))/2 + 0.5)
	}
	
	p = p + theme_bw(base_size = fontsize)
	
	return(p)
}

#' @rdname panel_boxplot
#' @export
panel_violin = function(data, fontsize = 10, par){
	if(!is.null(par$classes)){
			p = qplot(x, y, geom = "violin", fill = data[, par$classes], colour = I(NA), data = data)
	}
	else{
			p = qplot(x, y, geom = "violin", colour = I(NA), fill = I("gray30"), data = data) 
	}
	
	if(inherits(data, "twoListExpData")){
		p = p + geom_vline(xintercept = length(unique(data$x))/2 + 0.5)
	}
	
	p = p + theme_bw(base_size = fontsize)
	
	return(p)
}

#' @rdname panel_boxplot
#' @export
panel_violin_box = function(data, fontsize = 10, par){
	if(!is.null(par$classes)){
			p = qplot(x, y, geom = "violin", fill = data[, par$classes], colour = I(NA), data = data) + geom_boxplot(fill = NA, outlier.size = 0.25, size = I(0.1), width = 0.7)
	}
	else{
			p = qplot(x, y, geom = "violin", colour = I(NA), data = data) + geom_boxplot(fill = NA, outlier.size = 0.25, size = I(0.1), width = 0.7)
	}
	
	if(inherits(data, "twoListExpData")){
		p = p + geom_vline(xintercept = length(unique(data$x))/2 + 0.5)
	}
	
	p = p + theme_bw(base_size = fontsize)
	
	return(p)
}

##

## Panel functions for pca data
panel_histogram = function(data, fontsize = 10, par){
	if(!is.null(par$classes)){
		p = qplot(x, geom = "bar", fill = data[, par$classes], binwidth = (max(data$x) - min(data$x)) / 20, data = data) 
	}
	else{
		p = qplot(x, geom = "bar", data = data) + geom_rug()
	}
	p = p + theme_bw(base_size = fontsize)
	
	return(p)
}
##

## Panel function for dummyData
panel_dummy = function(data, fontsize = 10, par){
	if(nrow(data$mat) == 1){
		colors = "#336699"
		p = qplot(1, y, geom = "bar", stat = "identity", data = data$mat, fill = x, width = I(0.6), ylim = c(0, data$max), xlim = c(0.5, 1.5)) + scale_fill_manual(values = colors) + theme_bw(base_size = fontsize) + opts(legend.position = "none") + coord_flip()
	} 
	
	if(nrow(data$mat) == 2){
		colors = c("#336699", "#990033")
		data$mat$y[1] = -data$mat$y[1]
		data$mat$x = factor(data$mat$x, labels = c("G1 > G2", "G1 < G2"))
		p = qplot(1, y, geom = "bar", stat = "identity", position = "identity", data = data$mat, fill = x, ylim = c(-data$max, data$max), width = I(0.6), xlim = c(0.5, 1.5)) + scale_fill_manual("Regulation direction", values = colors) + theme_bw(base_size = fontsize)  + coord_flip()
	}
	
	return(p)
}

##

## Customization functions
customize_dummy = function(p, par){
	return(p)
}

 
#' Customization function for panel
#' 
#' This function is supposed to make small changes in the panel function appearance like 
#' changing color scheme for example. It has to match with the output of the 
#' corresponding panel function. Check examples in plot.gosummaries to see how to write 
#' one yourself.
#' 
#' @param p a ggplot2 plot object
#' @param par parameters object like in \code{\link{panel_boxplot}}
#' @return  a ggplot2 plot object with added customizations
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  # Example1
#' 
#' @export
customize = function(p, par){
	p = p + scale_fill_discrete(par$classes)
	return(p)
}
##

## plot.gosummaries

# TODO figure out the classes parameter situation either by creating good examples and explanatory text or changing the parameter

#' Plot the GOsummaries figure
#' 
#' The function to draw a GOsummaries figure based on a \code{gosummaries} object.  The 
#' GOsummaries figure consists of several components each defined by a gene list ora a 
#' pair of them. The GO annotations of them are shown as wordclouds. Optionally one can 
#' draw related (expression) data on panels atop of the wordclouds. 
#' 
#' In most cases the function can decide which type of plot to draw into the panel part. 
#' If there is no data explicitly put into the Data slots of the gosummaries object, it 
#' just draws a horizontal barplot with the numbers of genes. On visualizing the PCA 
#' data it draws histogram of the samples on the principal axes. For clustering and 
#' differential expression it draws the boxplot of expression values.   
#'
#' @param x a gosummaries object
#' @param components names of the gosummaries components to draw. The names are usually 
#' numbers 1:n in character format.
#' @param panel_plot plotting function for panel  
#' @param panel_customize customization function for the panel plot, menat for making 
#' small changes like changing colour scheme
#' @param panel_par list of arguments passed on to \code{panel_plot} function
#' @param panel_height panel height as number of lines, with given \code{fontsize}. If 
#' set to 0 no panel is drawn. 
#' @param fontsize font size used throughout the figure in points
#' @param term_length maximum length of the dispalyed GO categories in characters, 
#' longer names are cropped to this size
#' @param wordcloud_colors two element vector of colors to define color scheme for 
#' displaying the enrichment p-values across the wordclouds. First element defines the 
#' color for category with worst p-value and the second for the word with the best.
#' @param filename file path where to save the picture. Filetype is decided by the 
#' extension in the path. Currently following formats are supported: png, pdf, tiff, 
#' bmp, jpeg. Even if the plot does not fit into the plotting window, the file size is 
#' calculated so that the plot would fit there.
#' @param \dots not used
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  # Example1
#' 
#' @method plot gosummaries
#' 
#' @export
plot.gosummaries = function(x, components = names(x)[1:min(10, length(x))],  panel_plot = NULL, panel_customize = NULL, panel_par = list(), panel_height = 5, fontsize = 10, term_length = 35, wordcloud_colors = c("grey70", "grey10"), filename = NA, ...){
	
	# Check input
	if(!is.gosummaries(x)) stop("Function requires an object of gosummaries type")
	if(any(!(as.character(components) %in% names(x)))) stop("Selected components are not present in data")
	
	# Take out components of interest
	x = x[as.character(components)]
	if(length(x) < 1) stop("No components selected")
	
	# Add wordcloud colors and adjust the string length
	x = adjust_wordcloud_appearance(x, term_length, wordcloud_colors)
	
	# Attach default plotting method if it is not set
	if(is.null(panel_plot)){
		if (is.null(x[[1]]$Data)){
			panel_height = 0
		}
		else if(inherits(x[[1]]$Data, "dummyData")){
			if(panel_height == 5){
				panel_height = 3
			}
			panel_plot = panel_dummy
		}
		else if(inherits(x[[1]]$Data, "pcaData")){
			panel_plot = panel_histogram
		}
		else if(inherits(x[[1]]$Data, "oneListExpData") | inherits(x[[1]]$Data, "twoListExpData")){
			panel_plot = panel_boxplot
		}
		
		else{
			stop("The panel data is in unknown format, please specify matching panel_plot function")
		}
	}
	
	if(is.null(panel_customize)){
		if(inherits(x[[1]]$Data, "dummyData")){
			panel_customize = customize_dummy
		}
		else{
			panel_customize = customize
		}
	}
	
	
	# Set figure parameters
	par = list(
		fontsize = fontsize, 
		panel_height = panel_height, 
		panel_width = 1.35 * fontsize * term_length, 
		wordcloud_colors = wordcloud_colors
	)
	
	# Take the panel plot to proper format 
	plot_panel = panelize_ggplot2(panel_plot, panel_customize, panel_par)
	
	plot_motor(x, plot_panel = plot_panel, legend = T, par = par, filename = filename)
}

##

## Data type specific convenience functions like for prcomp, kmeans, limma, ...
# TODO implement the data type specific convenience functions maybe using a S3 generic function on these different data types
 
#' Prepare gosummaries object based on PCA results 
#' 
#' The PCA results are converted into a gosummaries object, by extracting genes with the 
#' largest positive and negative weights from each component. 
#' 
#' The usual visualisation of PCA results displays the projections of sample expression  
#' on the principal axes. It shows if and how the samples cluster, but not why do they  
#' behave like that. Actually, it is possible to go further and annotate the axes by 
#' studying genes that have the largest influence in the linear combinations that define 
#' the principal components. For example, high expression of genes with large negative 
#' weights pushes the samples projection to the negative side of the principal axis and 
#' large positive weigths to the positive side. If a sample has highly expressed genes 
#' in both groups it stays most probably in the middle. If we annotate functionally the 
#' genes with highest positive and negative weights for each of the principal axes, then 
#' it is possible to say which biological processes drive the separation of samples on 
#' them.   
#' 
#' This function creates a gosummaries object for such analysis. It expects the results 
#' of \code{\link{prcomp}} function. It assumes that the PCA was done on samples and, 
#' thus, the row names of the rotation matrix can be interpreted as gene names. For each 
#' component it annotates \code{n_genes} elements with highest positive and negative 
#' weights.
#'     
#'
#' @param x an object of class \code{prcomp}
#' @param annotation a \code{data.frame} describing the samples, its row names should match 
#' with column names of the projection matrix in x
#' @param components numeric vector of components to include 
#' @param n_genes how many top genes to use for annotation
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}}
#' @return  A gosummaries object.
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Example1
#' 
#' @method gosummaries prcomp
#' @S3method gosummaries prcomp
#' 
#' @export
gosummaries.prcomp = function(x, annotation = NULL, components = 1:6, n_genes = 500, ...){
	
	gl = list()
	for(i in components){
		gl[[sprintf("Principal component %d", i)]] = list(
			gl1 = rownames(x$rotation)[order((x$rotation[, components[i]]))][1:n_genes],
			gl2 = rownames(x$rotation)[order(-(x$rotation[, components[i]]))][1:n_genes]
		) 
	}
	
	cat("Annotating functionally\n")
	gosummaries = gosummaries.default(gl, ...)
	gosummaries = add_to_slot.gosummaries(gosummaries, "Percentage", paste(round(x$sdev[components] / sum(x$sdev) * 100), "%", sep = ""))
	gosummaries = add_pca.gosummaries(gosummaries, x, annotation)
	
	return(gosummaries)
}

 
#' Prepare gosummaries object based on k-means results 
#' 
#' The gosummaries object is created based on the genes in the clusters, it is possible 
#' to add corresponding gene expression data as well.
#' 
#' The k-means clustering of expression matrix naturally defines a set of gene lists 
#' that can be annotated functionally and displayed as a GOsummaries figure. This 
#' functon takes in a \code{kmeans} object and and converts it to a \code{gosummaries} 
#' object that can be plotted. If expression matrix is attached then the panel shows the 
#' expression values for each gene as boxplots, if not then number of genes is displayed
#' 
#' It is advisable to filter some genes out before doing the clustering since the very 
#' large gene lists (more than 2000 genes) might fail the annotation step and are 
#' usually not too specific either.  
#'
#' @param x an object of class \code{kmeans}
#' @param exp an expression matrix, with row names corresponding to the names of the 
#' genes in clusters (Optional)
#' @param annotation a \code{data.frame} describing the samples, its row names should 
#' match with column names of \code{exp} (Optional)
#' @param components numeric vector of clusters to annotate
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}} 
#' @return  A gosummaries object.
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Example1
#' 
#' @method gosummaries kmeans
#' @S3method gosummaries kmeans
#' 
#' @export
gosummaries.kmeans = function(x, exp = NULL, annotation = NULL, components = 1:length(x$size), ...){
	
	gl = list()
	for(i in components){
		gl[[sprintf("Cluster %d", i)]] = names(x$cluster[x$cluster == i])
	}
	
	cat("Annotating functionally\n")
	gosummaries = gosummaries.default(gl, ...)
	
	cat("Adding expression values\n")
	if(!is.null(exp)){
		gosummaries = add_expression.gosummaries(gosummaries, exp, annotation)
	}
	
	return(gosummaries)
}

 
#' Prepare gosummaries object based on limma results
#' 
#' The gosummaries object is created based on the differentially expresed genes, each 
#' contrast defines one component.
#' 
#' The usual differential expression analysis involves making several comparisons 
#' between treatments ehere each one yields an up and down regulated gene list. In a 
#' GOsummaries figure each comparison is displayed as one component with two wordclouds. 
#' If expression matrix is attached then the panel shows the expression values for each 
#' gene as boxplots, if not then number of genes is displayed 
#'
#' @param x an object of class \code{MArrayLM}
#' @param p.value p-value threshold as defined in topTable
#' @param lfc log fold change threshold as defined in topTable
#' @param adjust.method multiple testing adjustment method as defined in topTable
#' @param exp an expression matrix, with row names corresponding to the names of the 
#' genes in clusters (Optional)
#' @param annotation a \code{data.frame} describing the samples, its row names should 
#' match with column names of \code{exp} (Optional)
#' @param components numeric vector of comparisons to annotate
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}}
#' @return A gosummaries object.
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#' # Example1
#' 
#' @method gosummaries MArrayLM
#' @S3method gosummaries MArrayLM 
#' 
#' @export
gosummaries.MArrayLM = function(x, p.value = 0.05, lfc = 1, adjust.method = "fdr", exp = NULL, annotation = NULL, components = 1:ncol(x),  ...){
	
	# Calculate the gene list
	gl = list()
	perc = list()
	for(i in components){
		tt = topTable(x, coef = i, p.value = p.value, lfc = lfc, adjust.method = adjust.method, number = Inf)
		
		gl_up = as.character(tt$ID[tt$logFC > 0])
		gl_down = as.character(tt$ID[tt$logFC < 0])
		
		g1 = paste(rownames(x$contrasts)[x$contrasts[, i] < 0], collapse = ", ")
		g2 = paste(rownames(x$contrasts)[x$contrasts[, i] > 0], collapse = ", ")
		
		title = sprintf("%s VS %s", g1, g2)
		perc[[i]] = sprintf("G1 > G2: %d\nG1 < G2: %d", length(gl_down), length(gl_up))
		
		gl[[title]] = list(
			gl1 = gl_down,
			gl2 = gl_up
		) 
	}
	
	cat("Annotating functionally\n")
	gosummaries = gosummaries.default(gl, ...)
	
	cat("Adding expression values\n")
	if(!is.null(exp)){
		gosummaries = add_expression.gosummaries(gosummaries, exp, annotation)
	}
	gosummaries = add_to_slot.gosummaries(gosummaries, "Percentage", perc)
	
	return(gosummaries)
}



##
# 
# ## Test Examples
# # PCA
# library(AnnotatedPCA)
# library(gProfileR)
# data(lukk_small)
# pcr = prcomp(t(lukk_small$exp))
# 
# gp = gosummaries(pcr, annotation = lukk_small$annot)
# plot(gp)
# 
# # Kmeans
# sd = apply(lukk_small$exp, 1, sd)
# exp2 = lukk_small$exp[sd > 0.75,]
# exp2 = exp2 - apply(exp2, 1, mean)
# kmr2 = kmeans(exp2, centers = 6, iter.max = 100)
# 
# gp = gosummaries(kmr2)
# plot(gp)
# 
# gp = gosummaries(kmr2, exp = exp2)
# plot(gp)
# 
# gp = gosummaries(kmr2, exp = exp2, annotation = lukk_small$annot)
# plot(gp)
# 
# # Limma
# require(limma)
# 
# mm = model.matrix(~ lukk_small$annot$Tissue - 1)
# colnames(mm) = make.names(levels(factor(lukk_small$annot$Tissue)))
# 
# contrast = makeContrasts(brain - cell.line, hematopoietic.system - muscle, cell.line - hematopoietic.system, levels = colnames(mm))
# 
# fit = lmFit(lukk_small$exp, mm)
# fit = contrasts.fit(fit, contrast)
# fit = eBayes(fit)
# 
# gp_dumb = gosummaries(fit)
# plot(gp_dumb)
# 
# gp = gosummaries(fit, exp = lukk_small$exp, annotation = lukk_small$annot)
# plot(gp, panel_par = list(classes = "Tissue"))
# ##
# 
# 
# ## Ciesci
# setwd("~/Raivo/Projects/GOsummaries/CaseStudies")
# 
# exp = read.table("Data/ciesci_exp.tsv", sep = "\t")
# rownames(exp) = str_replace(rownames(exp), "_at", "")
# ann = read.table("Data/ciesci_ann.txt", sep = "\t", head = T, row.names = 1)
# ann = ann[colnames(exp), ]
# ann$Group = factor(ann$Group, c("mca_1", "pam_2", "msc_3", "carb_4"))
# exp = exp[, order(ann$Cancer_stage)]
# ann = ann[order(ann$Cancer_stage),]
# 
# customize_ciesci = function(p, par){
# 	p = p + scale_fill_discrete(par$classes, c = 80, l = 70, h.start = 30)
# 	return(p)
# }
# 
# # PCA
# pcr = prcomp(t(exp))
# 
# gp = gosummaries(pcr, annotation = ann, organism = "mmusculus")
# plot(gp, components = 1:3, panel_par = list(classes = "Group"), panel_customize = customize_ciesci, filename = "~/Desktop/GPplus/pca.pdf")

# # K-means
# sd = apply(exp, 1, sd)
# exp2 = exp[sd > 0.75,]
# exp2 = exp2 - apply(exp2, 1, mean)
# kmr2 = kmeans(exp2, centers = 6, iter.max = 100)
# 
# gp = gosummaries(kmr2, exp = exp2, annotation = ann, organism = "mmusculus")
# # plot(gp, panel_height = 0, filename = "~/Desktop/GPplus/kmeans_0.pdf")
# plot(gp, panel_par = list(classes = "Group"), panel_customize = customize_ciesci, filename = "~/Desktop/GPplus/kmeans.pdf")
# 
# # Limma 
# 
# library(limma)
# mm = model.matrix(~ ann$Group - 1)
# colnames(mm) = levels(ann$Group)
# 
# contrast = makeContrasts(pam_2 - mca_1, msc_3 - pam_2, carb_4 - msc_3, levels = colnames(mm))
# 
# fit = lmFit(exp, mm)
# fit = contrasts.fit(fit, contrast)
# fit = eBayes(fit)
# 
# gp = gosummaries(fit, exp = exp2, annotation = ann, organism = "mmusculus")
# 
# plot(gp, panel_height = 0, filename = "~/Desktop/GPplus/limma_0.pdf")
# plot(gp, panel_par = list(classes = "Group"), panel_customize = customize_ciesci, filename = "~/Desktop/GPplus/limma.pdf")
# 
# gp = add_dummydata.gosummaries(gp)
# plot(gp, panel_height = 3, panel_par = list(classes = "Group"), filename = "~/Desktop/GPplus/limma_bar.pdf")
# 
# 
# 
# 
# 
# ##
# 
# 
# 
