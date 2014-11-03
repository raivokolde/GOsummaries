## Dummydata object
dummyData = function(gl, max){
    res = list( 
        mat = data.frame(x = factor(1:length(gl)), 
                         y = unlist(lapply(gl, length))),
        max = max
    )
    
    class(res) = c(class(res), "dummyData")
    
    return(res)
}

is.dummyData = function(x) inherits(x, "dummyData")

add_dummydata.gosummaries = function(gosummaries){
    max = max(unlist(lapply(gosummaries, function(x){
        lapply(x$Gene_lists, length)
    })))
    
    for(i in seq_along(gosummaries)){
        gosummaries[[i]]$Data = dummyData(gosummaries[[i]]$Gene_lists, max)
    }
    
    return(gosummaries)
}
##

## gosummaries object constructor and related functions
gosummaries_base = function(gl = NULL, wc_data = NULL, score_type = "p-value", wc_algorithm = "middle", wordcloud_legend_title = NULL){
    # Check the input parameters for consistency
    if(is.null(gl) & is.null(wc_data)){
        stop("Either gene lists or the word cloud data have to be specified")
    }
    
    if(!(score_type %in% c("p-value", "count"))){
        stop("score_type has to be either: p-value or count")
    }
    
    if(!(wc_algorithm %in% c("middle", "top"))){
        stop("wc_algorithm has to be either: middle or top")
    }
    
    if(!is.null(gl) & !is.null(wc_data)){
        # Find length of the components
        if(length(gl) != length(wc_data)){
            stop("The gene lists and word cloud data have to have the same length")
        }
        
        components = 1:length(gl)
        
        # Find the number of lists per component (based on first component)
        gl_k = ifelse(is.list(gl[[1]]), 2, 1) 
        
        wc_condition = is.list(wc_data[[1]]) & !is.data.frame(wc_data[[1]])
        wc_data_k = ifelse(wc_condition, 2, 1)
        
        if(gl_k != wc_data_k){
            stop("The gene lists and the word cloud data have to have a similar structure")
        }
        else{
            k = gl_k
        }
        
        # Check if the names in the gene list and Go results are the same
        if(!identical(names(gl), names(wc_data))){
            stop("The gene lists and the word cloud data have to have same names")
        }
        else{
            names = names(gl)
        }
    }
    
    if(!is.null(gl) & is.null(wc_data)){
        components = 1:length(gl)
        k = ifelse(is.list(gl[[1]]), 2, 1) 
        names = names(gl)
    }
    
    if(is.null(gl) & !is.null(wc_data)){
        components = 1:length(wc_data)
        k = ifelse(is.list(wc_data[[1]]) & !is.data.frame(wc_data[[1]]), 2, 1)
        names = names(wc_data)
    }
    
    # Create the resulting data structure
    res = list()
    
    for(i in components){
        comp = list(
            Title = names[i],
            Gene_lists = NULL,
            WCD = switch(k, 
                list(wcd1 = wc_data[[i]]), 
                list(wcd1 = wc_data[[i]][[1]], wcd2 = wc_data[[i]][[2]])
            ),
            Data = NULL,
            Percentage = " "
        )
        
        if(!is.null(gl)){
            comp_gl = gl[[i]]
            if(k == 1){
                comp$Gene_lists = list(gl1 = comp_gl)
                comp$Percentage = sprintf("n: %d", length(comp_gl))
            }
            else{
                comp$Gene_lists = list(gl1 = comp_gl[[1]], gl2 = comp_gl[[2]])
                comp$Percentage = sprintf("G1: %d\nG2: %d",
                                           length(comp_gl[[1]]), 
                                          length(comp_gl[[2]]))
            }
        }
        
        res[[i]] = comp
    }
    
    # Set word cloud legend title
    if(is.null(wordcloud_legend_title)){
        if(score_type == "p-value"){
            wordcloud_legend_title = "P-value"
        }
        if(score_type == "count"){
            wordcloud_legend_title = "Count"
        }
    }
    
    res = structure(res, 
        class = "gosummaries", 
        score_type = score_type, 
        wordcloud_legend_title = wordcloud_legend_title, 
        wc_algorithm = wc_algorithm
    )
    
    return(res)
}

# gosummaries_base()
# gosummaries_base(list(Tere = "A"))
# gosummaries_base(list(Tere = list("A", "B")))
# gosummaries_base(NULL, list(Tere = data.frame(Term = "ime", Score = "muna")))
# gosummaries_base(list(Tere = "A"), list(Tere = data.frame(Term = "ime", Score = "muna")))
# gosummaries_base(list(Tere = "A"), list(Tore = data.frame(Term = "ime", Score = "muna")))
# gosummaries_base(list(Tere = list("A", "B")), list(Tere = data.frame(Term = "ime", Score = "muna")))

 
#' Constructor for gosummaries object
#' 
#' Constructor for gosummaries object that contains all the necessary 
#' information to draw the figure, like gene lists and their annotations, 
#' expression data and all the relevant texts.
#' 
#' The object is a list of "components", with each component defined by a gene 
#' list or a 
#' pair of gene lists. Each "component" has the slots as follows:
#' \itemize{
#'   \item \bold{Title}: title string of the component. (Default: the names of 
#' the gene lists)
#'   \item \bold{Gene_lists}: list of one or two gene lists
#'   \item \bold{WCD}: g:Profiler results based on the Gene_lists slot or user 
#' entered table. 
#'   \item \bold{Data}: the related data (expression values, PCA rotation, ...) 
#' that is used to draw the "panel" i.e. theplot above the wordclouds. In 
#' principle there is no limitation what  kind of data is there as far as the 
#' function that is provided to draw that in \code{\link{plot.gosummaries}} can 
#' use it.
#'   \item \bold{Percentage}: a text that is drawn on the right top corner of 
#' every component. In case of PCA this is the percentage of variation the 
#' component explains, by default it just depicts the number of genes in the 
#' Gene_lists slot.
#' }
#' 
#' Some visual parameters are stored in the attributes of \code{gosummaries} 
#' object: 
#' \code{score_type} tells how to handle the scores associated to wordclouds, 
#' \code{wc_algorithm} specifies the wordcloud layout algorithm and 
#' \code{wordcloud_legend_title} specifies the title of the wordcloud. One can 
#' change them using the \code{attr} function.
#' 
#' The word clouds are specified as \code{data.frame}s with two columns: "Term" 
#' and "Score". If one wants to use custom data for wordclouds, instead of the 
#' default GO enrichment results, then this is possible to specify parameter 
#' \code{wc_data}. The input structure is similar to the gene list input, only 
#' instead of gene lists one has the two column 
#' \code{data.frame}s.
#' 
#' The GO enrichment analysis is performed using g:Profiler web toolkit and its
#' associated R package \code{gProfileR}. This means the computer has to have 
#' internet access to annotate the gene lists. Since g:Profiler can accept a 
#' wide range of gene IDs then user usually does not have to worry about 
#' converitng the gene IDs into right format. To be absolutely sure the tool 
#' recognizes the gene IDs one can check if they will give any results in 
#' \url{http://biit.cs.ut.ee/gprofiler/gconvert.cgi}. 
#' 
#' There can be a lot of results for a typical GO enrichment analysis but 
#' usually these tend to be pretty redundant. Since one can fit only a small 
#' number of categories into a word cloud we have to bring down the number of 
#' categories to show an reduce the redundancy. For this we use hierarchical 
#' filtering option \"moderate\" in g:Profiler. In g:Profiler the categories 
#' are grouped together when they share one or more enriched parents. The 
#' \"moderate\" option selects the most significant category from each of such 
#' groups. (See more at http://biit.cs.ut.ee/gprofiler/)   
#' 
#' The slots of the object can be filled with custom information using a 
#' function \code{\link{add_to_slot.gosummaries}}. 
#' 
#' By default the Data slot is filled with a dataset that contains the number 
#' of genes in the Gene_lists slot. Expression data can be added to the object 
#' for example by using function \code{\link{add_expression.gosummaries}}. It 
#' is possible to derive your own format for the Data slot as well, as long as 
#' a panel plotting function for this data is alaso provided (See 
#' \code{\link{panel_boxplot}} for further information).
#' 
#' There are several constructors of gosummaries object that work on common 
#' analysis result objects, such as \code{\link{gosummaries.kmeans}}, 
#' \code{\link{gosummaries.MArrayLM}} and \code{\link{gosummaries.prcomp}} 
#' corresponding to k-means, limma and PCA results. 
#'
#' @param x list of arrays of gene names (or list of lists of arrays of gene 
#' names)
#' @param \dots additional parameters for gprofiler function 
#' @return   A gosummaries type of object
#' 
#' @seealso \code{\link{gosummaries.kmeans}}, 
#' \code{\link{gosummaries.MArrayLM}}, \code{\link{gosummaries.prcomp}}
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' # Define gene lists 
#' genes1 = c("203485_at", "209469_at", "209470_s_at", "203999_at", 
#' "205358_at", "203130_s_at", "210222_s_at", "202508_s_at", "203001_s_at", 
#' "207957_s_at", "203540_at", "203000_at", "219619_at", "221805_at", 
#' "214046_at", "213135_at", "203889_at", "209990_s_at", "210016_at", 
#' "202507_s_at", "209839_at", "204953_at", "209167_at", "209685_s_at",  
#' "211276_at", "202391_at", "205591_at", 
#' "201313_at")
#' genes2 = c("201890_at", "202503_s_at", "204170_s_at", "201291_s_at", 
#' "202589_at", "218499_at", "209773_s_at", "204026_s_at", "216237_s_at", 
#' "202546_at", "218883_s_at", "204285_s_at", "208659_at", "201292_at", 
#' "201664_at")
#' 
#' 
#' gl1 = list(List1 = genes1,  List2 = genes2) # One list per component
#' gl2 = list(List = list(genes1, genes2)) # Two lists per component
#' 
#' # Construct gosummaries objects
#' gs1 = gosummaries(gl1)
#' gs2 = gosummaries(gl2)
#' 
#' plot(gs1, fontsize = 8)
#' plot(gs2, fontsize = 8)
#' 
#' # Changing slot contents using using addToSlot.gosummaries 
#' gs1 = add_to_slot.gosummaries(gs1, "Title", list("Neurons", "Cell lines"))
#' 
#' # Adding expression data
#' data(tissue_example)
#' 
#' gs1 = add_expression.gosummaries(gs1, exp = tissue_example$exp, annotation = 
#' tissue_example$annot)
#' gs2 = add_expression.gosummaries(gs2, exp = tissue_example$exp, annotation = 
#' tissue_example$annot)
#'
#' plot(gs1, panel_par = list(classes = "Tissue"), fontsize = 8)
#' plot(gs2, panel_par = list(classes = "Tissue"), fontsize = 8)
#' }
#' 
#' # Using custom annotations for word clouds
#' wcd1 = data.frame(Term = c("KLF1", "KLF2", "POU5F1"), Score = c(0.05, 0.001, 
#' 0.0001))
#' wcd2 = data.frame(Term = c("CD8", "CD248", "CCL5"), Score = c(0.02, 0.005, 
#' 0.00001))
#' 
#' gs = gosummaries(wc_data = list(Results1 = wcd1, Results2 = wcd2))
#' plot(gs)
#' 
#' gs = gosummaries(wc_data = list(Results = list(wcd1, wcd2)))
#' plot(gs)
#' 
#' # Adjust wordcloud legend title
#' gs = gosummaries(wc_data = list(Results = list(wcd1, wcd2)), 
#' wordcloud_legend_title = "Significance score")
#' plot(gs)
#' 
#' @rdname gosummaries
#' @export
gosummaries = function(x = NULL, ...){
    UseMethod("gosummaries", x)
} 

#' @param wc_data precalculated GO enrichment results (see Details)
#' @param organism the organism that the gene lists correspond to. The format 
#' should be as follows: "hsapiens", "mmusculus", "scerevisiae", etc.
#' @param go_branches GO tree branches and pathway databases as denoted in 
#' g:Profiler (Possible values: BP, CC, MF, ke, re) 
#' @param max_p_value threshold for p-values that have been corrected for 
#' multiple testing
#' @param min_set_size minimal size of functional category to be considered
#' @param max_set_size maximal size of functional category to be considered
#' @param max_signif maximal number of categories returned per query
#' @param ordered_query logical showing if the lists are ordered or not (it 
#' determines if the ordered query algorithm is used in g:Profiler)
#' @param hier_filtering a type of hierarchical filtering used when reducing 
#' the number of g:Profiler results (see \code{\link{gprofiler}} for further 
#' information) 
#' @param score_type indicates the type of scores in \code{wc_data}. Possible 
#' values: "p-value" and "count"
#' @param wc_algorithm the type of wordcloud algorithm used. Possible values 
#' are "top" that puts first word to the top corner and "middle" that puts 
#' first word to the middle. 
#' @param wordcloud_legend_title title of the word cloud legend, should reflect 
#' the nature of the score
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' 
#' @rdname gosummaries
#' @method gosummaries default
#' @export
gosummaries.default = function(x = NULL, wc_data = NULL, organism = "hsapiens", go_branches = c("BP", "ke", "re"), max_p_value = 1e-2, min_set_size = 50, max_set_size = 1000, max_signif = 40, ordered_query = TRUE, hier_filtering = "moderate", score_type = "p-value", wc_algorithm = "middle", wordcloud_legend_title = NULL, ...){
    
    # Create basic structure
    res = gosummaries_base(gl = x, 
        wc_data = wc_data, 
        score_type = score_type, 
        wc_algorithm = wc_algorithm, 
        wordcloud_legend_title = wordcloud_legend_title
    )
    
    # Add data and annotations
    if(!is.null(x)){
        res = add_dummydata.gosummaries(res)
    }
    
    if(is.null(wc_data)){
        res = annotate.gosummaries(res, organism = organism, 
                                   go_branches = go_branches, 
                                   max_p_value = max_p_value, 
                                   min_set_size = min_set_size, 
                                   max_set_size = max_set_size, 
                                   max_signif = max_signif, 
                                   ordered_query = ordered_query, 
                                   hier_filtering = hier_filtering, ...)
    
        attr(res, "wordcloud_legend_title") = "Enrichment p-value"
    }
    
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
        for(l in a$WCD){
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

#' @param i index
#' 
#' @rdname add_to_slot.gosummaries
#' @method [ gosummaries
#' @export
"[.gosummaries" = function(x, i, ...) {
    attrs = attributes(x)
    out = unclass(x)
    out = out[i]
    attributes(out) = attrs
    
    return(out)
}


#' Functions for working with gosummaries object
#' 
#' Functions for working with gosummaries object 
#' 
#' Method [ ensures that subsetting gosummaries object will not lose the 
#' attributes.
#' 
#' \code{add_to_slot.gosummaries} allows to add values to specific slots in the 
#' gosummaries object
#'
#' @param x a gosummaries object
#' @param slot the component slot name to be filled (e.g Title, Percentage, 
#' etc.)
#' @param values list of values where each element corresponds to one component
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' data(gs_kmeans)
#' 
#' # Add new title to the components
#' gs_kmeans_new = add_to_slot.gosummaries(gs_kmeans, "Title", 
#' as.list(paste("K-means cluster", 1:length(gs_kmeans))))
#' 
#' print(gs_kmeans_new)
#' 
#' @rdname add_to_slot.gosummaries
#' @export
add_to_slot.gosummaries = function(x, slot, values){
    if(length(x) != length(values)){
        stop("Length of gosummaries object and values does not match")
    } 
    
    for(i in seq_along(x)){
        x[[i]][[slot]] = values[[i]]
    }
    
    return(x)
}
##

## Annotate gosummaries with g:Profiler
filter_gprofiler = function(gpr, go_branches, min_set_size, max_signif){
    gpr = gpr[gpr$domain %in% go_branches & gpr$term.size > min_set_size, ]
    
    gpr = plyr::ddply(gpr, "query.number", function(x) {
        x = x[order(x$p.value),]
        rank = rank(x$p.value)
        return(x[rank <= max_signif, ])
    })
    
    return(gpr)
}

annotate.gosummaries = function(gosummaries, organism, components = 1:length(gosummaries), go_branches, min_set_size, max_p_value, max_set_size, max_signif, ordered_query, hier_filtering, ...){
    
    if(!is.gosummaries(gosummaries)){
        stop("Function requires a gosummaries type of  object")
    } 
    
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
    user_agent = sprintf("gProfileR/%s; GOsummaries/%s", 
                            packageVersion("gProfileR"), 
                            packageVersion("GOsummaries"))
    gProfileR::set_user_agent(ua = user_agent, append = FALSE)
    
    gpr = gProfileR::gprofiler(query = gl, organism = organism, 
                               ordered_query = ordered_query, 
                               max_set_size = max_set_size, 
                               hier_filtering = hier_filtering, 
                               max_p_value = max_p_value, ...)
    
    # Clean and filter the results
    gpr$query.number = as.numeric(as.character(gpr$query.number))
    gpr = filter_gprofiler(gpr, go_branches = go_branches, 
                           min_set_size = min_set_size, max_signif = max_signif)
    gpr = gpr[, c("query.number", "p.value", "term.name")]
    colnames(gpr) = c("query.number", "Score", "Term")
    
    k = 1
    
    for(i in seq_along(components)){
        for(j in seq_along(gosummaries[[i]]$WCD)){
            lname = paste("wcd", j, sep = "")
            results = gpr[gpr$query.number == k, c("Term", "Score")] 
            gosummaries[[i]]$WCD[[lname]] = results
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
#' The data is added to the object in a "long" format so it would be directly 
#' usable by the ggplot2 based panel drawing functions 
#' \code{\link{panel_boxplot}} etc. For each component it produces a data frame 
#' with columns:
#' \itemize{
#'   \item \bold{x} : sample IDs for the x axis, their factor order is the same 
#' as on the columns of input matrix 
#'   \item \bold{y} : expression values from the matrix
#'   \item  . . . : sample annotation columns from the annotation table that 
#' can be displayed on figure as colors.
#' } 
#'
#' @param gosummaries a gosummaries object
#' @param exp an expression matrix, with row names corresponding to the names 
#' in the Gene_lists slot
#' @param annotation a \code{data.frame} describing the samples, its row names 
#' should match with column names of \code{exp}
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' data(gs_limma)
#' data(tissue_example)
#' 
#' # Add just expression without annotations
#' gs_limma_exp1 = add_expression.gosummaries(gs_limma, exp = 
#' tissue_example$exp)
#' 
#' print(gs_limma_exp1)
#' 
#' # Add expression with annotations
#' gs_limma_exp2 = add_expression.gosummaries(gs_limma, exp = 
#' tissue_example$exp, annotation = tissue_example$annot)
#' 
#' print(gs_limma_exp1)
#' }
#' @export
add_expression.gosummaries = function(gosummaries, exp, annotation = NULL){
    # Check arguments
    if(!is.gosummaries(gosummaries)){
        stop("Function requires a gosummaries type object")
    } 
    
    if(!all(colnames(exp) %in% rownames(annotation)) & !is.null(annotation)){
        stop("Column names of expression matrix and row names of annotation data.frame do not match")
    } 
    
    if(is.null(gosummaries[[1]]$Gene_lists)){
        stop("No gene list data provided")
    }
    
    # Do the work
    for(i in seq_along(gosummaries)){
        a = list()
        for(j in seq_along(gosummaries[[i]]$Gene_lists)){
            e = data.frame(exp[gosummaries[[i]]$Gene_lists[[j]], ])
            d = reshape2::melt(e)
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
    if(!is.gosummaries(gosummaries)){ 
        stop("Function requires a gosummaries type object")
    }
    
    if(inherits(pcr, "prcomp")){
        x = pcr$x
    }
    else{
        x = pcr
    }
    
    if(!all(rownames(x) %in% rownames(annotation)) & !is.null(annotation)) {
        stop("Column names of expression matrix and row names of annotation dataframe do not match")
    }
    
    for(i in seq_along(gosummaries)){
        a = data.frame(x = x[, i])
        
        if(!is.null(annotation)){
            a$ID = rownames(x)
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
    
    ind = n > max
    strings[ind] = paste(substr(strings[ind], 1, max - 3), "...", sep = "")
    
    return(strings)
}

convert_scores = function(gosummaries){
    for(i in seq_along(gosummaries)){
        for(j in seq_along(gosummaries[[i]]$WCD)){
            scores = gosummaries[[i]]$WCD[[j]]$Score
            
            if(attr(gosummaries, "score_type") == "p-value"){
                if(any(scores == 0)){
                    scores[scores == 0] = min(scores[scores != 0])
                }
                scores = -log10(scores)
            }
            if(attr(gosummaries, "score_type") == "count"){
                scores = scores
            }
            
            gosummaries[[i]]$WCD[[j]]$Score = scores
        }
    }
    
    return(gosummaries)
}

adjust_wordcloud_appearance = function(gosummaries, term_length = 35, wordcloud_colors = c("grey70", "grey10")){
    
    # Convert different types of scores to common values
    gosummaries = convert_scores(gosummaries)
    
    # Find best score
    scores = plyr::llply(gosummaries, function(x){
        comp_scores = plyr::llply(x$WCD, function(y){
            y$Score
        })
        
        unlist(comp_scores)
    })
    
    scores = unlist(scores)
    
    if(length(scores) == 0){
        return(gosummaries)
    }
     
    best_score = max(scores)
    
    # Adjust Term names and scores
    for(i in seq_along(gosummaries)){
        for(j in seq_along(gosummaries[[i]]$WCD)){
            # Shorten GO names
            original_terms = gosummaries[[i]]$WCD[[j]]$Term
            shortened_terms = shorten_strings(original_terms, term_length)
            gosummaries[[i]]$WCD[[j]]$Term = shortened_terms
            
            # Calculate colors
            comp_scores = gosummaries[[i]]$WCD[[j]]$Score
            palette = colorRampPalette(wordcloud_colors)(100)
            palette_ind = ceiling(comp_scores / best_score * 100)
            gosummaries[[i]]$WCD[[j]]$Colors = palette[palette_ind]
        }
    }
    
    return(gosummaries)
}
##

## Plotting utility functions

## Define zeroGrob
zeroGrob = function(){
    grob(cl = "zeroGrob", name = "NULL")
}
widthDetails.zeroGrob = 
heightDetails.zeroGrob = 
grobWidth.zeroGrob = 
grobHeight.zeroGrob = function(x) unit(0, "cm")
drawDetails.zeroGrob = function(x, recording) {}
is.zero = function(x) identical(x, zeroGrob())
##

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
        tiff = function(x, ...) tiff(x, units = "in", res = 300, 
                                        compression = "lzw", ...),
        bmp = function(x, ...) bmp(x, units = "in", res = 300, ...),
        stop("File type should be: pdf, png, bmp, jpg, tiff")
    )
    
    f(filename, width = width, height = height)
}

panelize_ggplot2 = function(plot_function, customize_function, par){
    res = function(data, fontsize, legend = FALSE){
        p = plot_function(data, fontsize, par)
        p = customize_function(p, par)
        p = ggplot2::ggplot_build(p)
        p = ggplot2::ggplot_gtable(p) 
        
        if(legend){
            if(any(grepl("guide-box", p$layout$name))){
                gb = gtable::gtable_filter(p, "guide-box")
                gb_grob = gb$grob[[1]]$grob[[1]]
                nc = ncol(gb_grob)
                nr = nrow(gb_grob)
                res = gb_grob[2:(nr - 1), 2:(nc - 1)]
                return(res)
            }
            else{
                res = gtable::gtable(widths = unit(0, "cm"), 
                                     heights = unit(0, "cm"))
                return(res)
            }
        }
        else{
            return(gtable::gtable_filter(p, "panel"))
        }
    }
    
    return(res)
} 
##

## Raw plotting functions
calc_component_dimensions = function(component, par){
    # Wordcloud height
    nr = max(plyr::laply(component$WCD, nrow))
    if(length(component$WCD) > 1){
        wc_height = ifelse(nr > 3, max(nr / 6.5, 3), nr)
        # arrows_height = 1.5
        arrows_height = 0.5
    }
    else{
        wc_height = ifelse(nr > 3, max(nr / 10, 3), nr)
        arrows_height = 0.5
    }
    
    # Percentage slot width
    gp = gpar(fontsize = par$fontsize, cex = 0.8)
    percentage_rows = strsplit(component$Percentage, "\n")[[1]]
    percentage_grobs = lapply(as.list(percentage_rows), textGrob, gp = gp)
    percentage_widths = lapply(percentage_grobs, grobWidth)
    percentage_widths_units = do.call(unit.c, percentage_widths)
    
    # Compile results
    lines_in_points = par$fontsize * 1.445
    res = list(
        title_height = unit(1.5 * lines_in_points, "points"),
        panel_height = unit(par$panel_height * lines_in_points, "points"),
        arrows_height = unit(arrows_height * lines_in_points, "points"),
        wc_height = unit(wc_height * lines_in_points, "points"),
    
        panel_width = unit(par$panel_width * lines_in_points, "points"),
        percentage_width = max(percentage_widths_units) * 1.25,
        wc_width = unit(par$panel_width * lines_in_points / 
                        length(component$WCD), "points")
    )
    
    return(res)
}

calc_components_dimensions = function(gosummaries, par){
    component_dims = list()
    
    for(i in seq_along(gosummaries)){
        component_dims[[i]] = calc_component_dimensions(gosummaries[[i]], par)
    }
    
    return(component_dims)
}

gen_legend = function(legend_data, par){
    n = length(legend_data$colors)
    # Create Grobs
    title = textGrob(legend_data$title, y = 1, x = 0, just = c(0, 1),
                     gp = gpar(fontsize = par$fontsize, fontface = "bold",
                    cex = 0.8))
    
    # rect_height = unit(1.7 * par$fontsize, "points")
    rect_height = unit(6.096, "mm")
    yy = unit(1, "npc") - unit(0.8 * par$fontsize * 1.1, "points") - 
         (0:(n - 1)) * rect_height
    boxes = rectGrob(x = 0, y = yy, height = rect_height, width = rect_height,
                      just = c(0, 1), 
                     gp = gpar(col = 0, fill = legend_data$colors))
    
    yyy = yy - rect_height * 0.5
    gl = gList()
    length = c(rep(0, n), convertWidth(grobWidth(title), "in"))
    for(i in 1:n){
        gl[[i]] = textGrob(legend_data$labels[[i]], 
                           x = rect_height + unit(3, "points"), y = yyy[i], 
                           hjust = 0, 
                           gp = gpar(cex = 0.8, fontsize = par$fontsize))
        length[i] = convertWidth(grobWidth(gl[[i]]), "in")
    }
    
    # Calculate size
    height = unit(0.8 * par$fontsize * 1.445, "points") + n * rect_height
    width = rect_height + unit(3, "points") + unit(max(length), "in")
     
    # Put together a frame
    fg = frameGrob()
    fg = packGrob(fg, rectGrob(width = width, height = height, 
                  gp = gpar(col = NA)))
    fg = packGrob(fg, title)
    fg = packGrob(fg, boxes)
    for(i in 1:length(gl)){
        fg = packGrob(fg, gl[[i]])
    }
    
    return(fg)
}

gen_wordcloud_legend = function(gosummaries, par){
    legend_data = list()
    
    legend_data$title = par$wordcloud_legend_title
    legend_data$colors = colorRampPalette(rev(par$wordcloud_colors))(3)
    
    # Calculate p-value breakpoints
    scores = plyr::ldply(gosummaries, function(x){
        plyr::ldply(x$WCD, function(y) data.frame(y$Score))
    })[,2]
    
    if(length(scores) == 0){
        empty_grob = rectGrob(width = unit(0.0001, "cm"), 
                              height = unit(0.0001, "cm"))
        return(empty_grob)
    }
    
    best_score = max(scores)    
    
    breaks = grid.pretty(c(0, best_score))
    if(length(breaks) %% 2 == 0){
        average_low = breaks[length(breaks) / 2]
        average_high = breaks[length(breaks) / 2 + 1]
        average = mean(c(average_low, average_high))
    }
    else{
        average = breaks[ceiling(length(breaks) / 2)]
    }
    
    if(attr(gosummaries, "score_type") == "p-value"){
        legend_data$labels = c(
            substitute(10 ^ -p, list(p = breaks[length(breaks)])), 
            substitute(10 ^ -p, list(p = average)), 
            1
        )
    }
    if(attr(gosummaries, "score_type") == "count"){
        legend_data$labels = c(breaks[length(breaks)], average, 0)
    }
    
    return(gen_legend(legend_data, par))
}

plot_wordcloud = function(words, freq, color, algorithm, dimensions){
    if(length(words) > 0){
        return(plotWordcloud(words, freq, colors = color, random.order = FALSE,
                              min.freq = -Inf, rot.per = 0, scale = 0.85,
                             max_min = c(1, 0), algorithm = algorithm, 
                             add = FALSE, grob = TRUE, 
                             dimensions = dimensions))
    }
    
    return(zeroGrob())
}

plot_arrow = function(end, par){
    x = switch(end, first = c(0, 0.95), both = c(0.05, 0.95), last = c(0.05, 1))
    res = linesGrob(x = x, y = 0.5, arrow = arrow(ends = end, type = "closed",
                     angle = 15, length = unit(0.1, "inches")), 
                    gp = gpar(lwd = 0.3 * par$fontsize, col = "grey40"))
    return(res)
}

plot_component = function(data_component, plot_panel, par, component_dims){
    
    # Create gtable
    heights = with(component_dims, unit.c(title_height, panel_height, 
                                          arrows_height + wc_height))
    widths = with(component_dims, unit.c(panel_width, percentage_width))
    
    gtable_component = gtable::gtable(widths, heights)
    
    # Add title
    title = textGrob(x = 0, 
        hjust = 0, 
        label = data_component$Title, 
        gp = gpar(fontface = "bold", fontsize = par$fontsize)
    )
    
    gtable_component = gtable::gtable_add_grob(gtable_component, title, 1, 1,
                                                 clip = "off")
    
    # Add plot
    if(par$panel_height != 0){
        p = plot_panel(data_component$Data, par$fontsize)
    }
    else{
        p = zeroGrob()
    }
    b = rectGrob(gp = gpar(lwd = 1.5, col = "grey40", fill = NA))
    
    gtable_component = gtable::gtable_add_grob(gtable_component, 
                                               gTree(children = gList(p, b)), 
                                               2, 1, clip = "off")
    
    # Add percentage
    p = textGrob(data_component$Percentage, x = 0.1, y = 1, vjust = 1, 
                 hjust = 0, gp = gpar(fontsize = par$fontsize, cex = 0.8))
    
    gtable_component = gtable::gtable_add_grob(gtable_component, p, 2, 2, 
                                               clip = "off")
    
    # Arrows and wordclouds
    heights = with(component_dims, unit.c(arrows_height, wc_height))
    widths = with(component_dims, rep(wc_width, length(data_component$WCD)))
    
    gtable_aw = gtable::gtable(widths, heights)
    
    if(length(data_component$WCD) == 1){
        wc = plot_wordcloud(words = data_component$WCD$wcd1$Term,
            freq = data_component$WCD$wcd1$Score, 
            color = data_component$WCD$wcd1$Colors, 
            algorithm = switch(par$wc_algorithm, top = "leftside_top", 
                               middle = "leftside"), 
            dimensions = with(component_dims, unit.c(wc_width, wc_height))
        )
        
        gtable_aw = gtable::gtable_add_grob(gtable_aw, wc, 2, 1, 
                                            name = "wordcloud")
    }
    
    if(length(data_component$WCD) == 2){
        wc1 = plot_wordcloud(words = data_component$WCD$wcd1$Term,
            freq = data_component$WCD$wcd1$Score, 
            color = data_component$WCD$wcd1$Colors, 
            algorithm = switch(par$wc_algorithm, top = "leftside_top", 
                                middle = "leftside"), 
            dimensions = with(component_dims, unit.c(wc_width, wc_height))
        )
        wc2 = plot_wordcloud(words = data_component$WCD$wcd2$Term,
            freq = data_component$WCD$wcd2$Score, 
            color = data_component$WCD$wcd2$Colors, 
            algorithm = switch(par$wc_algorithm, top = "rightside_top", 
                               middle = "rightside"), 
            dimensions = with(component_dims, unit.c(wc_width, wc_height))
        )
        
        gtable_aw = gtable::gtable_add_grob(gtable_aw, wc1, 2, 1,
                                             name = "wordcloud-left")
        gtable_aw = gtable::gtable_add_grob(gtable_aw, wc2, 2, 2, 
                                            name = "wordcloud-right")
    }
    
        
    gtable_component = gtable::gtable_add_grob(gtable_component, gtable_aw, 3,
                                               1, name = "arrows-wordcloud")
    gtable_component = gtable::gtable_add_padding(gtable_component, 
                        unit(c(0, 0, 0.5 * par$fontsize * 1.445, 0), "points"))
    
    return(gtable_component)

}

plot_motor = function(gosummaries, plot_panel, par = list(fontsize = 10, panel_height = 5, panel_width = 405), filename = NA){
    
    # Calculate dimensions for the picture components
    component_dimensions = calc_components_dimensions(gosummaries, par)
    
    # Create component grobs
    components = list()
    par$wc_algorithm = attr(gosummaries, "wc_algorithm")
    for(i in seq_along(gosummaries)){
        components[[i]] = plot_component(data_component = gosummaries[[i]], 
            plot_panel = plot_panel, 
            par = par, 
            component_dims = component_dimensions[[i]]
        )
    }
    
    # Create legends 
    if(par$panel_height != 0){
        panel_legend = plot_panel(gosummaries[[1]]$Data, par$fontsize, 
                                  legend = TRUE)
        height = gtable::gtable_height(panel_legend)
        height_in_cm = convertHeight(height, "cm", valueOnly = TRUE)
        if(height_in_cm != 0){
            panel_legend = gtable::gtable_add_padding(panel_legend, 
                        unit(c(0, 0, 0.5 * par$fontsize * 1.445, 0), "points"))
        }
    }
    else{
        panel_legend = gtable::gtable(widths = unit(0, "cm"), 
                                      heights = unit(0, "cm"))
    }
    
    if(par$wordcloud_colors[1] == par$wordcloud_colors[2]){
        wordcloud_legend = gtable::gtable(widths = unit(0.001, "cm"), 
                                          heights = unit(0.001, "cm"))
    }
    else{
        wordcloud_legend = gen_wordcloud_legend(gosummaries, par)
    }
    
    # Calculate legend dimensions to create gtable for it
    pl_height = gtable::gtable_height(panel_legend)
    pl_width = gtable::gtable_width(panel_legend)
    wl_height = grobHeight(wordcloud_legend)
    wl_width = grobWidth(wordcloud_legend)
    
    legend_width = max(pl_width, wl_width)
    legend_height = unit.c(pl_height, wl_height)
    vp = viewport(
        y = unit(1, "npc") - unit(1.5 * par$fontsize * 1.445, "points"), 
        height = sum(legend_height), just = c(0.5, 1))
    
    gtable_legend = gtable::gtable(width = legend_width, 
                                   height = legend_height, vp = vp)
    
    panel_legend = editGrob(panel_legend, vp = viewport(x = unit(0, "npc"),
                             width = pl_width, just = c(0, 0.5)))
    gtable_legend = gtable::gtable_add_grob(gtable_legend, 
                                            grobs = panel_legend, t = 1, l = 1,
                                            name = "panel-legend", clip = "off")
    
    wordcloud_legend = editGrob(wordcloud_legend, 
                                vp = viewport(x = unit(0, "npc"), 
                                            width = wl_width, just = c(0, 0.5)))
    gtable_legend = gtable::gtable_add_grob(gtable_legend, wordcloud_legend, 
        t = 2, 
        l = 1,
        name = "wordcloud-legend", 
        clip = "off"
    )
    
    
    # Create gtable layout for the whole figure
    component_widths = do.call(unit.c, lapply(components, gtable::gtable_width))
    widths = unit.c(max(component_widths), gtable::gtable_width(gtable_legend))
    heights = do.call(unit.c, lapply(components, gtable::gtable_height))
    
    gtable_full = gtable::gtable(widths, heights)
    
    # Add components
    for(i in 1:length(components)){
        component_width = gtable::gtable_width(components[[i]])
        components[[i]] = editGrob(components[[i]], 
                                   vp = viewport(x = 0, 
                                                  width = component_width, 
                                                  just = c(0, 0.5)))
        component_name = paste("Component", i, sep = "-")
        gtable_full = gtable::gtable_add_grob(gtable_full, components[[i]], i, 
                                              1, name = component_name)
    }
    
    # Add legend
    gtable_full = gtable::gtable_add_grob(gtable_full, gtable_legend, 1, 2, 
                                          length(components), name = "legend")
    
    # Add padding
    padding_size = unit(0.5 * par$fontsize * 1.445, "points")
    gtable_full = gtable::gtable_add_padding(gtable_full, padding_size)
    
    # Open connection to file if filename specified
    if(!is.na(filename)){
        width = convertWidth(gtable::gtable_width(gtable_full), "inches",
                              valueOnly = TRUE) 
        height = convertHeight(gtable::gtable_height(gtable_full), "inches",
                                valueOnly = TRUE) 
        open_file_con(filename, width, height)
    }
    
    # Draw
    grid.draw(gtable_full)
    
    # Close connection if filename specified
    if(!is.na(filename)){
        dev.off()
    }
    
    return(gtable_full)
}
##

## Panel functions for expression data
panel_crossbar = function(data, fontsize = 10, par){
    qq = function(x, ...) {
        res = data.frame(ymin = quantile(x, 0.05), y = median(x), 
                         ymax = quantile(x, 0.95))
        
        return(res)
    }
    
    if(!is.null(par$classes)){
        p = ggplot2::qplot(x = data$x, y = data$y, geom = "crossbar", 
                           stat = "summary", width = 0.7, fun.data = qq, 
                           fill = data[, par$classes]) + 
            ggplot2::geom_bar(aes(x = 1, y = 0, fill = data[, par$classes]), 
                              stat = "identity")
    }
    else{
        p = ggplot2::qplot(x = data$x, y = data$y, geom = "crossbar", 
                           stat = "summary", width = 0.7, fun.data = qq)
    }
    
    if(inherits(data, "twoListExpData")){
        p = p + ggplot2::geom_vline(xintercept = length(unique(data$x))/2 + 0.5)
    }
    
    p = p + ggplot2::theme_bw(base_size = fontsize)
    
    return(p)
}

 
#' Panel drawing functions 
#' 
#' These functions are used to draw the panel portion of every component based 
#' on the Data slots in gosummaries object. These concrete functions assume the 
#' data is presented as done by \code{\link{add_expression.gosummaries}}. They 
#' provide three  options: boxplot, violin plot (which shows the distrubution 
#' more precisely) and both combined. Additionally it is possible to combine 
#' the values from one each functional group into one box/violin plot using 
#' the corresponding functions ending with "_classes".
#' 
#' These functions specify in principle the general setting for the panels, 
#' like which "geom"-s, how the data is transformed and summarized, etc. To 
#' make small adjustments to the figure such as changing color scheme, write 
#' your own customization function (See \code{\link{customize}} as example).
#' 
#' It is possible to write your own panel plotting function, as long as the 
#' parameters used and the return value are similar to what is specified here. 
#' When writing a new panel function one only has to make sure that it matches 
#' the data given in the Data slot of the gosummaries object.
#' 
#' @param data the data from Data slot of the gosummaries object
#' @param fontsize fontsize in points, mainly used to ensure that the legend 
#' fontsizes match
#' @param par additional parameters for drawing the plot, given as list. These 
#' functions use only \code{classes} slot for figuring out which parameter to 
#' use for coloring the "geom"-s. However, when building a custom function it 
#' provides a way to give extra parameters to the plotting function. 
#' @return  It returns a function that can draw a ggplot2 plot of the data in 
#' Data slot of a gosummaries object. The legend and the actual plots for the 
#' panels are extracted later from the figure produced by this function.
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' 
#' @examples
#'  
#' \dontrun{
#' data(gs_kmeans)
#' 
#' # Draw default version with plot_boxplot
#' plot(gs_kmeans, components = 1:3, classes = "Tissue")
#' 
#' # Define alternative where one boxplot summarises all values in one class
#' plot_classes = function(data, fontsize, par){
#'     qplot(x = data[, par$classes], y = data$y, geom = "boxplot", 
#'           fill = data[, par$classes]) + theme_bw()
#' }
#' 
#' plot(gs_kmeans, components = 1:3, panel_plot = plot_classes, classes = "Tissue")
#' 
#' 
#' # Flip the boxplots to make them more comparable
#' plot_classes = function(data, fontsize, par){
#'     qplot(x = data[, par$classes], y = data$y, geom = "boxplot", 
#'         fill = data[, par$classes]) + coord_flip() + theme_bw()
#' }
#' 
#' plot(gs_kmeans, components = 1:3, panel_plot = plot_classes, classes = "Tissue")
#' }
#' 
#' @rdname panel_boxplot
#' @export
panel_boxplot = function(data, fontsize = 10, par){
    qq = function(x) {
        bs = boxplot.stats(x)$stats 
        data.frame(ymin = bs[1], lower = bs[2], middle = bs[3], upper = bs[4], 
                   ymax = bs[5])
    
    }
    if(!is.null(par$classes)){
            p = ggplot2::qplot(x = data$x, y = data$y, geom = "boxplot", 
                               stat = "summary", fun.data = qq, 
                               fill = data[, par$classes], width = 0.7) 
    }
    else{
            p = ggplot2::qplot(x = data$x, y = data$y, geom = "boxplot", 
                               stat = "summary", fun.data = qq, width = 0.7) 
    }
    
    if(inherits(data, "twoListExpData")){
        n_samples = length(unique(data$x))/2
        p = p + ggplot2::geom_vline(xintercept = n_samples + 0.5)
    }
    
    p = p + ggplot2::theme_bw(base_size = fontsize)
    
    return(p)
}

#' @rdname panel_boxplot
#' @export
panel_violin = function(data, fontsize = 10, par){
    if(!is.null(par$classes)){
        p = ggplot2::qplot(x = data$x, y = data$y, geom = "violin", 
                           fill = data[, par$classes], colour = I(NA))
    }
    else{
        p = ggplot2::qplot(x = data$x, y = data$y, geom = "violin", 
                           colour = I(NA), fill = I("gray30")) 
    }
    
    if(inherits(data, "twoListExpData")){
        n_samples = length(unique(data$x))/2
        p = p + ggplot2::geom_vline(xintercept = n_samples + 0.5)
    }
    
    p = p + ggplot2::theme_bw(base_size = fontsize)
    
    return(p)
}

#' @rdname panel_boxplot
#' @export
panel_violin_box = function(data, fontsize = 10, par){
    if(!is.null(par$classes)){
            p = ggplot2::qplot(x = data$x, y = data$y, geom = "violin", 
                               fill = data[, par$classes], colour = I(NA)) + 
                ggplot2::geom_boxplot(fill = NA, outlier.size = 0.25, 
                                      size = I(0.1), width = 0.7)
    }
    else{
            p = ggplot2::qplot(x = data$x, y = data$y, geom = "violin", 
                               colour = I(NA)) + 
                ggplot2::geom_boxplot(fill = NA, outlier.size = 0.25, 
                                      size = I(0.1), width = 0.7)
    }
    
    if(inherits(data, "twoListExpData")){
        n_samples = length(unique(data$x))/2
        p = p + ggplot2::geom_vline(xintercept = n_samples + 0.5)
    }
    
    p = p + ggplot2::theme_bw(base_size = fontsize)
    
    return(p)
}

combine_classes = function(data, par){
    if(inherits(data, "twoListExpData")){
        data[, par$classes] = factor(data[, par$classes])
        new_x = paste(substr(data[, "x"], 1, 1), data[, par$classes])
        n_levels = length(levels(data[, par$classes]))
        levels = paste(rep(1:2, times = c(n_levels, n_levels)),
                        levels(data[, par$classes]))
        
        data$x = factor(new_x, levels = levels)
    }
    else{
        data$x = data[, par$classes]
    }
    
    return(data)
}

#' @rdname panel_boxplot
#' @export
panel_boxplot_classes = function(data, fontsize = 10, par){
    data = combine_classes(data, par)
    
    p = panel_boxplot(data, fontsize, par)
    
    return(p)
}

#' @rdname panel_boxplot
#' @export
panel_violin_classes = function(data, fontsize = 10, par){
    data = combine_classes(data, par)
    
    p = panel_violin(data, fontsize, par)
    
    return(p)
}


#' @rdname panel_boxplot
#' @export
panel_violin_box_classes = function(data, fontsize = 10, par){
    data = combine_classes(data, par)
    
    p = panel_violin_box(data, fontsize, par)
    
    return(p)
}


##

## Panel functions for pca data
panel_histogram = function(data, fontsize = 10, par){
    if(!is.null(par$classes)){
        p = ggplot2::qplot(data$x, geom = "bar", fill = data[, par$classes], 
                           binwidth = (max(data$x) - min(data$x)) / 20, 
                           colour = I("grey70")) 
    }
    else{
        p = ggplot2::qplot(data$x, geom = "bar", 
                           binwidth = (max(data$x) - min(data$x)) / 20, 
                           data = data)
    }
    p = p + ggplot2::theme_bw(base_size = fontsize)
    
    return(p)
}
##

## Panel function for dummyData
panel_dummy = function(data, fontsize = 10, par){
    if(nrow(data$mat) == 1){
        colors = "#336699"
        p = ggplot2::qplot(1, data$mat$y, geom = "bar", stat = "identity", 
                           fill = data$mat$x, width = I(0.6), 
                           ylim = c(0, data$max), xlim = c(0.5, 1.5)) + 
            ggplot2::scale_fill_manual(values = colors) + 
            ggplot2::theme_bw(base_size = fontsize) + 
            ggplot2::theme(legend.position = "none") + 
            ggplot2::coord_flip()
    } 
    
    if(nrow(data$mat) == 2){
        colors = c("#336699", "#990033")
        data$mat$y[1] = -data$mat$y[1]
        data$mat$x = factor(data$mat$x, labels = c("G1 > G2", "G1 < G2"))
        p = ggplot2::qplot(1, data$mat$y, geom = "bar", stat = "identity", 
                           position = "identity", fill = data$mat$x, 
                           ylim = c(-data$max, data$max), width = I(0.6), 
                           xlim = c(0.5, 1.5)) + 
            ggplot2::scale_fill_manual("Regulation direction", values = colors)+
            ggplot2::theme_bw(base_size = fontsize) + 
            ggplot2::theme(legend.position = "none") + 
            ggplot2::coord_flip()
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
#' This function is supposed to make small changes in the panel function 
#' appearance like changing color scheme for example. It has to match with the 
#' output of the corresponding panel function. Check examples in 
#' plot.gosummaries to see how to write one yourself.
#' 
#' @param p a ggplot2 plot object
#' @param par parameters object like in \code{\link{panel_boxplot}}
#' @return  a ggplot2 plot object with added customizations
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' data(gs_limma_exp)
#' 
#' cust = function(p, par){
#'     p = p + scale_fill_brewer(par$classes, type = "qual", palette = 1)
#'     return(p)
#' }
#' 
#' plot(gs_limma_exp, classes = "Tissue", panel_plot = panel_boxplot, 
#'         panel_customize = cust, fontsize = 8) 
#' } 
#' 
#' @export
customize = function(p, par){
    p = p + ggplot2::scale_fill_discrete(par$classes)
    return(p)
}
##

## plot.gosummaries


#' Plot the GOsummaries figure
#' 
#' The function to draw a GOsummaries figure based on a \code{gosummaries} 
#' object.  The GOsummaries figure consists of several components each defined 
#' by a gene list ora a pair of them. The GO annotations of them are shown as 
#' wordclouds. Optionally one can draw related (expression) data on panels atop 
#' of the wordclouds. 
#' 
#' In most cases the function can decide which type of plot to draw into the 
#' panel part. If there is no data explicitly put into the Data slots of the 
#' gosummaries object, it just draws a horizontal barplot with the numbers of 
#' genes. On visualizing the PCA data it draws histogram of the samples on the 
#' principal axes. For clustering and differential expression it draws the 
#' boxplot of expression values.   
#'
#' @param x a gosummaries object
#' @param components index for the components to draw. 
#' @param classes name of the variable from annotation data.frame that defines 
#' the colors in the plot
#' @param panel_plot plotting function for panel  
#' @param panel_customize customization function for the panel plot, menat for 
#' making small changes like changing colour scheme
#' @param panel_par list of arguments passed on to \code{panel_plot} function
#' @param panel_height panel height as number of lines, with given 
#' \code{fontsize}. If set to 0 no panel is drawn. 
#' @param panel_width panel width in lines of text
#' @param fontsize font size used throughout the figure in points
#' @param term_length maximum length of the dispalyed GO categories in 
#' characters, longer names are cropped to this size
#' @param wordcloud_colors two element vector of colors to define color scheme 
#' for displaying the enrichment p-values across the wordclouds. First element 
#' defines the color for category with worst p-value and the second for the 
#' word with the best. Set the same value for both if you want to remove the 
#' color scale and the legend. 
#' @param wordcloud_legend_title title of wordcloud legend 
#' @param filename file path where to save the picture. Filetype is decided by 
#' the extension in the path. Currently following formats are supported: png, 
#' pdf, tiff, bmp, jpeg. Even if the plot does not fit into the plotting 
#' window, the file size is calculated so that the plot would fit there.
#' @param \dots not used
#' 
#' @return The \code{\link{gtable}} object containing the figure
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' data(gs_limma)
#' 
#' # Default plot
#' plot(gs_limma, fontsize = 8)
#' 
#' # Omitting the panel area 
#' plot(gs_limma, panel_height = 0, fontsize = 8)
#' 
#' # Selecting only certain components
#' plot(gs_limma, components = c(1, 3), fontsize = 8)
#' 
#' # Cutting the longer terms shorter (see right wordcloud on first component)
#' plot(gs_limma, term_length = 20, fontsize = 8) 
#' 
#' # Change wordcloud colors
#' plot(gs_limma, term_length = 20, wordcloud_colors = c("#C6DBEF", "#08306B"), 
#'      fontsize = 8)
#' 
#' # Adjust panel plot type (see panel_boxplot help for options)
#' data(gs_kmeans)
#' 
#' plot(gs_kmeans, panel_plot = panel_violin, classes = "Tissue", components = 
#'      1:2, fontsize = 8)
#' plot(gs_kmeans, panel_plot = panel_violin_box, classes = "Tissue", 
#'      components = 1:2, fontsize = 8)
#' 
#' # Adjust colorscheme for plot (see customize help for more information) 
#' cust = function(p, par){
#'   p = p + scale_fill_brewer(par$classes, type = "qual", palette = 2)
#'   return(p)
#' }
#' plot(gs_kmeans, panel_plot = panel_violin, panel_customize = cust, 
#'      classes = "Tissue", components = 1:2, fontsize = 8)
#' }
#' @method plot gosummaries
#' 
#' @export
plot.gosummaries = function(x, components = 1:min(10, length(x)), classes = NA, panel_plot = NULL, panel_customize = NULL, panel_par = list(), panel_height = 5, panel_width = 30, fontsize = 10, term_length = 35, wordcloud_colors = c("grey70", "grey10"), wordcloud_legend_title = NULL, filename = NA, ...){
    
    # Check input
    if(!is.gosummaries(x)){
        stop("Function requires an object of gosummaries type")
    } 
    if(any(!(components %in% 1:length(x)))){
        stop("Selected components are not present in data")
    } 
    
    # Add classes to panel_par
    if(!is.na(classes)){
        if(!(classes %in% colnames(x[[1]]$Data))){
            stop("Classes variable has to be present in the data.frame in the component Data slot")
        }
        else{
            panel_par[["classes"]] = classes
        }
    }
    
    # Take out components of interest
    x = x[components]
    if(length(x) < 1) stop("No components selected")
    
    # Add wordcloud colors and adjust the string length
    x = adjust_wordcloud_appearance(x, term_length, wordcloud_colors)
    
    # Adjust wordcloud legend title 
    if(is.null(wordcloud_legend_title)){
        wordcloud_legend_title = attr(x, "wordcloud_legend_title")
    } 
    
    
    # Attach default plotting method if it is not set
    first_comp_data = x[[1]]$Data
    if(is.null(panel_plot)){
        if (is.null(first_comp_data)){
            panel_height = 0
        }
        else if(inherits(first_comp_data, "dummyData")){
            if(panel_height == 5){
                panel_height = 3
            }
            panel_plot = panel_dummy
        }
        else if(inherits(first_comp_data, "pcaData")){
            panel_plot = panel_histogram
        }
        else if(inherits(first_comp_data, "oneListExpData") | 
                inherits(first_comp_data, "twoListExpData")){
            panel_plot = panel_boxplot
        }
        
        else{
            stop("The panel data is in unknown format, please specify matching panel_plot function")
        }
    }
    
    if(is.null(panel_customize)){
        if(inherits(first_comp_data, "dummyData")){
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
        panel_width = panel_width, 
        wordcloud_colors = wordcloud_colors,
        wordcloud_legend_title = wordcloud_legend_title
    )
    
    # Take the panel plot to proper format 
    plot_panel = panelize_ggplot2(panel_plot, panel_customize, panel_par)
    
    invisible(plot_motor(x, plot_panel = plot_panel, par = par, 
                         filename = filename))
}

##

## Data type specific convenience functions like for prcomp, kmeans, limma, ...
convert_gene_ids = function(unique_ids, gconvert_target, organism){
    if(!is.null(gconvert_target)){
        cat(sprintf("%s\n", "Convert IDs"))
        gcr = gProfileR::gconvert(unique_ids, target = gconvert_target, 
                                  organism = organism)
        gcr = ddply(gcr, "alias", function(x){
            if(nrow(x) == 1){
                return(x)
            }
            else{
                return(x[which.min(nchar(as.character(x$target))), ])
            }
        })
        i2g = gcr$target
        names(i2g) = gcr$alias
    }
    else{
        i2g = unique_ids
        names(i2g) = toupper(unique_ids)
    }
    
    return(i2g)
}

filter_wc_data = function(wc_data, i2g, n_genes){
    wc_data$Term = i2g[toupper(as.character(wc_data$Term))]
    wc_data = wc_data[!is.na(wc_data$Term),]
    wc_data = wc_data[!duplicated(wc_data$Term), ]
    if(nrow(wc_data) != 0){
        wc_data = wc_data[1:min(nrow(wc_data), n_genes), ]
    }
    
    return(wc_data)
}

filter_pca_wc_data = function(wc_data){
    index = wc_data$Score / max(wc_data$Score)
    
    return(subset(wc_data, index > 0.1))
}

# Spearman correlation analysis on PCA components
pspearman = function(rho, n, lower.tail = TRUE) {
    q = (n^3 - n) * (1 - rho) / 6
    den = (n * (n^2 - 1)) / 6 + 1 
    r = 1 - q/den
    p = pt(r / sqrt((1 - r ^ 2) / (n - 2)), df = n - 2, 
           lower.tail = !lower.tail)
    
    return(p)
}

spearman_mds = function(pc, expression, n_genes){
    n = ncol(expression)
    cc = cor(t(expression), pc, method = "spearman")[,1]
    res = data.frame(Term = names(cc), Correlation = cc)
    res$Score = pmin(pspearman(res$Correlation, n, lower.tail = FALSE), 
                     pspearman(res$Correlation, n, lower.tail = TRUE))
    res = res[order(abs(res$Correlation), decreasing = TRUE), ]
    
    n_up = sum(res$Correlation > 0)
    n_down = sum(res$Correlation < 0)
    
    res = list(
        res[res$Correlation < 0, c("Term", "Score")][1:min(n_down, 2*n_genes),],
        res[res$Correlation > 0, c("Term", "Score")][1:min(n_up, 2 * n_genes), ]
    )
    
    return(res)
}

 
#' Prepare gosummaries object based on Multi Dimensional Scaling (MDS) results
#' 
#' The Multi Dimensional Scaling (MDS) results are converted into a gosummaries 
#' object, by finding genes that have most significant Spearman correlations 
#' with each component.  
#' 
#' This visualisation of MDS results is very similar to the one performed by
#'  \code{\link{gosummaries.prcomp}}. Difference from PCA is that, in general,
#'  we do not have the loadings for individual genes that could be used to 
#' associate genes with components. However, it is possible to find genes that
#'  are most correlated with each component. This function uses Spearman 
#' correlation coefficient to find most correlated features. The significance 
#' of the correlation values is decided using he approximation with 
#' t-distribution.
#' 
#' The function can also display genes instead of their GO annotations, while 
#' the sizes of the gene names correspond to the Spearman correlation p-values. 
#' The corresponding parameters are described in more detail in 
#' \code{\link{gosummaries.MArrayLM}}. This feature is important in 
#' applications, like metabolomics and metagenomics, where the features are not 
#' genes and it is not possible to run GO enrichment analysis.
#' 
#' @param x a matrix representation of multi dimensional scaling result, rows
#' correspond to samples
#' @param exp an expression matrix, with columns corresponding to samples 
#' (these have to be in the same order as in \code{x})
#' @param annotation a \code{data.frame} describing the samples, its row names 
#' should match with column names of \code{exp} (Optional)
#' @param components numeric vector of comparisons to annotate
#' @param show_genes logical showing if GO categories or actual genes are shown 
#' in word clouds
#' @param gconvert_target specifies gene ID format for genes showed in word 
#' cloud. The name of the format is passed to \code{\link{gconvert}}, if NULL 
#' original IDs are shown.
#' @param n_genes maximum number of genes shown in a word cloud
#' @param organism the organism that the gene lists correspond to. The format 
#' should be as follows: "hsapiens", "mmusculus", "scerevisiae", etc.
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}}
#' @return A gosummaries object.
#' 
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' library(vegan)
#' 
#' data("metagenomic_example")
#' 
#' # Run Principal Coordinate Analysis on Bray-Curtis dissimilarity matrix 
#' pcoa = cmdscale(vegdist(t(metagenomic_example$otu), "bray"), k = 3)
#' 
#' # By turning off the GO analysis we can show the names of taxa
#' gs = gosummaries(pcoa, metagenomic_example$otu, metagenomic_example$annot, 
#'                  show_genes = TRUE, gconvert_target = NULL, n_genes = 30)
#' 
#' plot(gs, class = "BodySite", fontsize = 8)
#' }
#' 
#' @export
gosummaries.matrix = function(x, exp = NULL, annotation = NULL, components = 1:min(ncol(x), 10), show_genes = FALSE, gconvert_target = "NAME", n_genes = ifelse(show_genes, 30, 500), organism = "hsapiens", ...){
    # Check assumptions
    if(is.null(exp)){
        stop("expression matrix has to be specified")
    }
    
    if(nrow(x) != ncol(exp)){
        stop("expression matrix has to have the same number of columns as the MDS matrix rows")
    }
    
    # Create wordcloud data
    wc_data = list()
    for(i in components){
        component_name = sprintf("Component %d", i)
        wc_data[[component_name]] = spearman_mds(x[, i], exp, n_genes)
    }
    
    # Create gosummaries object
    if(!show_genes){
        gl = plyr::llply(wc_data, plyr::llply, function(x){
            as.character(x$Term)[1:n_genes]
        })
        
        gosummaries = gosummaries.default(gl, organism = organism, ...)
    }
    else{
        unique_ids = unlist(plyr::llply(wc_data, plyr::llply, function(x){
            as.character(x$Term)[1:n_genes]
        }))
        i2g = convert_gene_ids(unique_ids, gconvert_target, organism)
        
        wc_data = plyr::llply(wc_data, plyr::llply, function(x){
            filter_wc_data(x, i2g, n_genes)
        })
        
        gosummaries = gosummaries_base(gl = NULL, wc_data = wc_data, 
                                       wc_algorithm = "top", 
                                       score_type = "p-value",
                                    wordcloud_legend_title = "Spearman p-value")
    }
    
    # Add histogram data 
    gosummaries = add_pca.gosummaries(gosummaries, x, annotation)
    gosummaries = add_to_slot.gosummaries(gosummaries, "Percentage", 
                                          rep("   ", length(components)))
    
    return(gosummaries)
}

 
#' Prepare gosummaries object based on PCA results 
#' 
#' The PCA results are converted into a gosummaries object, by extracting genes with the largest positive and negative weights from each component. 
#' 
#' The usual visualisation of PCA results displays the projections of sample 
#' expression on the principal axes. It shows if and how the samples cluster, 
#' but not why do they behave like that. Actually, it is possible to go 
#' further and annotate the axes by studying genes that have the largest 
#' influence in the linear combinations that define the principal components. 
#' For example, high expression of genes with large negative weights pushes 
#' the samples projection to the negative side of the principal axis and large 
#' positive weigths to the positive side. If a sample has highly expressed 
#' genes in both groups it stays most probably in the middle. If we annotate 
#' functionally the genes with highest positive and negative weights for each 
#' of the principal axes, then it is possible to say which biological 
#' processes drive the separation of samples on them.   
#' 
#' This function creates a gosummaries object for such analysis. It expects 
#' the results of \code{\link{prcomp}} function. It assumes that the PCA was 
#' done on samples and, thus, the row names of the rotation matrix can be 
#' interpreted as gene names. For each component it annotates \code{n_genes} 
#' elements with highest positive and negative weights.
#' 
#' The function can also display genes instead of their GO annotations, while 
#' the sizes of the gene names correspond to the PCA loadings. The 
#' corresponding parameters are described in more detail in 
#' \code{\link{gosummaries.MArrayLM}}.
#'     
#' @param x an object of class \code{prcomp}
#' @param annotation a \code{data.frame} describing the samples, its row names 
#' should match with column names of the projection matrix in x
#' @param components numeric vector of components to include 
#' @param show_genes logical showing if GO categories or actual genes are 
#' shown in word clouds
#' @param gconvert_target specifies gene ID format for genes showed in word 
#' cloud. The name of the format is passed to \code{\link{gconvert}}, if NULL 
#' original IDs are shown.
#' @param n_genes shows the number of genes used for annotating the component, 
#' in case gene names are shown, it is the maximum number of genes shown in a 
#' word cloud
#' @param organism the organism that the gene lists correspond to. The format 
#' should be as follows: "hsapiens", "mmusculus", "scerevisiae", etc
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}}
#' @return  A gosummaries object.
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' data(tissue_example)
#' 
#' pcr = prcomp(t(tissue_example$exp))
#' gs_pca = gosummaries(pcr, annotation = tissue_example$annot)
#' 
#' plot(gs_pca, classes = "Tissue", components = 1:3, fontsize = 8)
#' }
#' 
#' # Read metabolomic data
#' data(metabolomic_example)
#' 
#' pca = prcomp(t(metabolomic_example$data))
#' 
#' # Turn off GO enricment, since it does not work on metabolites
#' gs = gosummaries(pca, annotation = metabolomic_example$annot, 
#'                  show_gene = TRUE, gconvert_target = NULL)
#' plot(gs, class = "Tissue", components = 1:3, fontsize = 8)
#' 
#' @method gosummaries prcomp
#' 
#' @export
gosummaries.prcomp = function(x, annotation = NULL, components = 1:10, show_genes = FALSE, gconvert_target = "NAME", n_genes = ifelse(show_genes, 30, 500), organism = "hsapiens", ...){
    
    gl = list()
    for(i in components){
        comp_name = sprintf("Principal component %d", i)
        comp_weights = x$rotation[, i]
        genes = rownames(x$rotation)
        
        gl[[comp_name]] = list(
            gl1 = genes[order((comp_weights))][1:n_genes],
            gl2 = genes[order(-(comp_weights))][1:n_genes]
        ) 
    }
    
    if(show_genes){
        unique_ids = unique(unlist(llply(gl, llply, function(x){
            x[1:min(length(x), 2 * n_genes)]
        })))
        
        i2g = convert_gene_ids(unique_ids, gconvert_target, organism)
        
        # Create tables for wordclouds
        wc_data = list()
        for(i in components){
            comp_weights = x$rotation[, i]
            genes = rownames(x$rotation)
            
            n_up = sum(comp_weights > 0)
            ind_up = order(-(comp_weights))[1:min(n_up, 2 * n_genes)]
            wc_data_up = data.frame(
                Term = genes[ind_up], 
                Score = comp_weights[ind_up]
            )
            
            n_down = sum(comp_weights < 0)
            ind_down = order((comp_weights))[1:min(n_down, 2 * n_genes)]
            wc_data_down = data.frame(
                Term = genes[ind_down], 
                Score = -comp_weights[ind_down]
            )
            
            wc_data_up = filter_wc_data(wc_data_up, i2g, n_genes)
            wc_data_down = filter_wc_data(wc_data_down, i2g, n_genes)
            
            wc_data_up = filter_pca_wc_data(wc_data_up)
            wc_data_down = filter_pca_wc_data(wc_data_down)
            
            wc_data[[sprintf("Principal component %d", i)]] = list(
                wc_data_down,
                wc_data_up
            ) 
        }
        
        # Create gosummaries object
        gosummaries = gosummaries_base(gl = gl, wc_data = wc_data, 
                                       wc_algorithm = "top", 
                                       score_type = "count", 
                                       wordcloud_legend_title = "PCA weight")
    }
    else{
        gosummaries = gosummaries.default(gl, organism = organism,  ...)
    }
    
    percentages = round((x$sdev ** 2)[components] / sum(x$sdev ** 2) * 100)
    percentages = paste(percentages, "%", sep = "")
    gosummaries = add_to_slot.gosummaries(gosummaries, "Percentage", 
                                          percentages)
    gosummaries = add_pca.gosummaries(gosummaries, x, annotation)
    
    return(gosummaries)
}

 
#' Prepare gosummaries object based on k-means results 
#' 
#' The gosummaries object is created based on the genes in the clusters, it is 
#' possible to add corresponding gene expression data as well.
#' 
#' The k-means clustering of expression matrix naturally defines a set of gene 
#' lists that can be annotated functionally and displayed as a GOsummaries  
#' figure. This functon takes in a \code{kmeans} object and and converts it to 
#' a \code{gosummaries} object that can be plotted. If expression matrix is 
#' attached then the panel shows the expression values for each gene as 
#' boxplots, if not then number of genes is displayed
#' 
#' It is advisable to filter some genes out before doing the clustering since 
#' the very large gene lists (more than 2000 genes) might fail the annotation 
#' step and are usually not too specific either.  
#'
#' @param x an object of class \code{kmeans}
#' @param exp an expression matrix, with row names corresponding to the names 
#' of the genes in clusters (Optional)
#' @param annotation a \code{data.frame} describing the samples, its row names 
#' should match with column names of \code{exp} (Optional)
#' @param components numeric vector of clusters to annotate
#' @param organism the organism that the gene lists correspond to. The format 
#' should be as follows: "hsapiens", "mmusculus", "scerevisiae", etc.
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}} 
#' @return  A gosummaries object.
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' \dontrun{
#' data(tissue_example)
#' 
#' # Filter genes and perform k-means
#' sd = apply(tissue_example$exp, 1, sd)
#' exp2 = tissue_example$exp[sd > 0.75,]
#' exp2 = exp2 - apply(exp2, 1, mean)
#' kmr = kmeans(exp2, centers = 6, iter.max = 100)
#' 
#' # Create gosummaries object  
#' gs_kmeans = gosummaries(kmr, exp = exp2, annotation = tissue_example$annot)
#' plot(gs_kmeans, panel_height = 0, components = 1:3, fontsize = 8)
#' plot(gs_kmeans, classes = "Tissue", components = 1:3, fontsize = 8)
#' }
#' 
#' @method gosummaries kmeans
#' 
#' @export
gosummaries.kmeans = function(x, exp = NULL, annotation = NULL, components = 1:length(x$size), organism = "hsapiens", ...){
    
    gl = list()
    for(i in components){
        gl[[sprintf("Cluster %d", i)]] = names(x$cluster[x$cluster == i])
    }
    
    gosummaries = gosummaries.default(gl, organism = organism, 
                                      ordered_query = FALSE,  ...)
    
    if(!is.null(exp)){
        gosummaries = add_expression.gosummaries(gosummaries, exp, annotation)
    }
    
    return(gosummaries)
}

 
#' Prepare gosummaries object based on limma results
#' 
#' The gosummaries object is created based on the differentially expresed 
#' genes, each contrast defines one component.
#' 
#' The usual differential expression analysis involves making several 
#' comparisons between treatments ehere each one yields an up and down 
#' regulated gene list. In a GOsummaries figure each comparison is displayed as 
#' one component with two wordclouds. If expression matrix is attached then the 
#' panel shows the expression values for each gene as boxplots, if not then 
#' number of genes is displayed 
#' 
#' It is possible to show the gene names instead of GO annotations in the 
#' wordclouds. The word sizes in wordclouds are defined by the limma p-values. 
#' As the gene identifiers in expression matrices are usually rather 
#' unintelligible then they are automatically converted into gene names using  
#' \code{\link{gconvert}} function. It is possible to show also the original 
#' identifiers by setting \code{gconvert_target} to NULL. This can be useful if
#' the values do not correspond to genes, but for example metabolites.  
#'
#' @param x an object of class \code{MArrayLM}
#' @param p.value p-value threshold as defined in topTable
#' @param lfc log fold change threshold as defined in topTable
#' @param adjust.method multiple testing adjustment method as defined in 
#' topTable
#' @param exp an expression matrix, with row names corresponding to the names 
#' of the genes in clusters (Optional)
#' @param annotation a \code{data.frame} describing the samples, its row names 
#' should match with column names of \code{exp} (Optional)
#' @param components numeric vector of comparisons to annotate
#' @param show_genes logical showing if GO categories or actual genes are 
#' shown in word clouds
#' @param gconvert_target specifies gene ID format for genes showed in word 
#' cloud. The name of the format is passed to \code{\link{gconvert}}, if NULL 
#' original IDs are shown.
#' @param n_genes maximum number of genes shown in a word cloud
#' @param organism the organism that the gene lists correspond to. The format 
#' should be as follows: "hsapiens", "mmusculus", "scerevisiae", etc.
#' @param \dots GO annotation filtering parameters as defined in 
#' \code{\link{gosummaries.default}}
#' @return A gosummaries object.
#' @author  Raivo Kolde <raivo.kolde@@eesti.ee>
#' @examples
#' 
#' \dontrun{
#' data(tissue_example)
#' 
#' # Do the t-test comparisons
#' mm = model.matrix(~ factor(tissue_example$annot$Tissue) - 1)
#' colnames(mm) = make.names(levels(factor(tissue_example$annot$Tissue)))
#' 
#' contrast = limma::makeContrasts(brain - cell.line, 
#'                                 hematopoietic.system - muscle, 
#'                                 cell.line - hematopoietic.system, 
#'                                 levels = colnames(mm))
#' 
#' fit = limma::lmFit(tissue_example$exp, mm)
#' fit = limma::contrasts.fit(fit, contrast)
#' fit = limma::eBayes(fit)
#' 
#' gs_limma = gosummaries(fit)
#' gs_limma_exp = gosummaries(fit, exp = tissue_example$exp, 
#'                            annotation = tissue_example$annot)
#' 
#' plot(gs_limma, fontsize = 8)
#' plot(gs_limma, panel_height = 0, fontsize = 8)
#' plot(gs_limma_exp, classes = "Tissue", fontsize = 8)
#' }
#' 
#' @method gosummaries MArrayLM
#' 
#' @export
gosummaries.MArrayLM = function(x, p.value = 0.05, lfc = 1, adjust.method = "fdr", exp = NULL, annotation = NULL, components = 1:ncol(x), show_genes = FALSE, gconvert_target = "NAME", n_genes = 30, organism = "hsapiens",  ...){
    
    # Calculate the gene list
    gl = list()
    perc = list()
    for(i in components){
        flevels = rownames(x$contrasts)
        
        tt = limma::topTable(x, coef = i, p.value = p.value, lfc = lfc,
                             adjust.method = adjust.method, number = Inf)
        
        if(nrow(tt) == 0){
            tt = data.frame(logFC = numeric(0), AveExpr = numeric(0), t = numeric(0), adj.P.Val = numeric(0), P.Value = numeric(0), B = numeric(0))
        }
        
        tt$ID = rownames(tt)
        
        gl_up = as.character(tt$ID[tt$logFC > 0])
        gl_down = as.character(tt$ID[tt$logFC < 0])
        
        g1 = paste(flevels[x$contrasts[, i] < 0], collapse = ", ")
        g2 = paste(flevels[x$contrasts[, i] > 0], collapse = ", ")
        
        title = sprintf("G1: %s; G2: %s", g1, g2)
        perc[[title]] = sprintf("G1 > G2: %d\nG1 < G2: %d", length(gl_down), 
                            length(gl_up))
        
        gl[[title]] = list(
            gl1 = gl_down,
            gl2 = gl_up
        ) 
    }
    
    # Either add gene names or GO categories
    if(show_genes){
        # Convert IDs
        # unique_ids = unique(unlist(gl))
        unique_ids = unique(unlist(llply(gl, llply, function(x){
            x[1:min(length(x), 2 * n_genes)]
        })))
        
        i2g = convert_gene_ids(unique_ids, gconvert_target, organism)
        
        # Create tables for wordclouds
        wc_data = list()
        flevels = rownames(x$contrasts)
        for(i in components){
            
            tt = limma::topTable(x, coef = i, p.value = p.value, lfc = lfc, 
                                 adjust.method = adjust.method, number = Inf)
            if(nrow(tt) == 0){
                tt = data.frame(logFC = numeric(0), AveExpr = numeric(0), t = numeric(0), adj.P.Val = numeric(0), P.Value = numeric(0), B = numeric(0))
            }
            
            tt$ID = rownames(tt)
            
            tt = tt[!is.na(tt$adj.P.Val), ]
            
            wc_data_up = tt[tt$logFC > 0, c("ID", "adj.P.Val")]
            wc_data_down = tt[tt$logFC < 0, c("ID", "adj.P.Val")]
            
            colnames(wc_data_up) = c("Term", "Score")
            colnames(wc_data_down) = c("Term", "Score")
            
            wc_data_up = filter_wc_data(wc_data_up, i2g, n_genes)
            wc_data_down = filter_wc_data(wc_data_down, i2g, n_genes)
            
            g1 = paste(flevels[x$contrasts[, i] < 0], collapse = ", ")
            g2 = paste(flevels[x$contrasts[, i] > 0], collapse = ", ")
        
            title = sprintf("G1: %s; G2: %s", g1, g2)
            
            wc_data[[title]] = list(
                wc_data_down,
                wc_data_up
            ) 
        }
        
        # Create gosummaries object
        gosummaries = gosummaries_base(gl = gl, wc_data = wc_data, 
                                       wc_algorithm = "top")
        gosummaries = add_dummydata.gosummaries(gosummaries)
    }
    else{
        gosummaries = gosummaries.default(gl, organism = organism, ...)
    }
    
    # Add additional data
    if(!is.null(exp)){
        gosummaries = add_expression.gosummaries(gosummaries, exp, annotation)
    }
    
    gosummaries = add_to_slot.gosummaries(gosummaries, "Percentage", perc)
    
    return(gosummaries)
}

##