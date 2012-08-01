findCoordinates = function (width, height){
	.Call("findCoordinates", width, height, PACKAGE = "GOsummaries")
}

findCoordinates_left = function (width, height){
	.Call("findCoordinates_left", width, height, PACKAGE = "GOsummaries")
}


#' Plot a wordcloud given words and frequencies
#' 
#' Uses the algorithm from wordcloud package to calculate the positions of the words. 
#' then uses grid graphics to plot the words on screen. The shape of the wordcloud 
#' depends on the shape of the plotting window
#'
#' @param words vector of words to draw
#' @param freq frequencies for words, has to be the same length as words vector
#' @param rot.per percentage of vertical words
#' @param max_min relative scales to adjust the size difference between largest and 
#' smallest word, by default the largest word is written with 10 times as large font than 
#' the smallest
#' @param scale a fraction of the available space on figure that will be covered with the bounding boxes of words
#' @param min.freq minimal frequency of words to be displayed
#' @param max.words maximal number of words to be displayed
#' @param random.order plot words in random order. If false, they will be plotted in
#'  decreasing frequency
#' @param colors vector of colors fro the words. This vector will be extrapolated into as 
#' many colors as needed, starting with the first color for lower frequencies and ending 
#' with last color for higher frequencies.
#' @param random.colors if true, assigns random color for the words.
#' @param algorithm algorithm to find positions of words possible values: "circle", "leftside" and "rightside".
#' @param tryfit if TRUE the algorithm checks if all words fit to the figure, if not it
#' tries gradually smaller values of scale parameter until everything fits
#' @param add if TRUE adds the picture to existing plot.
#' 
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' @examples
#'  plotWordcloud(c("Audi", "Volkswagen", "Opel", "Porsche", "Mercedez", "BMW"), 8:3)
#' 
#' @export
plotWordcloud = function(words, freq, rot.per = 0.3, max_min = c(1, 0.1), scale = 0.4, min.freq = 3, max.words = Inf, random.order = FALSE, colors = "black", random.colors = FALSE, algorithm = "circle", tryfit = T, add = F){
	# Empty the drawng area
	if(add){
		grid.newpage()
	}
	
	# Check if word and freq are same length
	if(length(words) != length(freq)){
		stop("The word and freq length do not match")
	}
	
	# Create the data frame for words
	d = data.frame(words, freq)
	
	# Order words according to frequencies
	ord = order(-d$freq)
	d = d[ord, ]
	
	# Filter the words and frequencies
	d = d[1:nrow(d) <= max.words, ]
	d = d[d$freq >= min.freq, ]
	
	# Create colors
	if(length(colors) == 1){
		d$colors = colors
	}
	else if(random.colors){
		d$colors = sample(colorRampPalette(rev(colors))(nrow(d)))
	}
	else if(length(colors) == nrow(d)){
		d$colors = colors
	}
	else{
		d$colors = rep(colorRampPalette(rev(colors))(length(unique(d$freq))), table(d$freq))
	}
	
	# Randomize order
	if(random.order){
		d = d[sample(1:nrow(d)), ]
	}
	
	# Normalize the text sizes
	d$words = as.character(d$words)
	normedFreq <- d$freq/max(d$freq)
	d$size <- (max_min[1]-max_min[2])*normedFreq + max_min[2]
	
	# Decide angle
	n = nrow(d)
	d$angle = ifelse(runif(n) < rot.per, 90, 0)
	
	# Calculate the word sizes
	w = unit(rep(1, n), "strwidth", as.list(as.character(d$words)))
	d$width = convertWidth(w, "npc", valueOnly = T)
	h = unit(rep(1, n), "strheight", as.list(as.character(d$words)))
	d$height = convertHeight(h, "npc", valueOnly = T)
	
	tailed = grepl("g|j|p|q|y|_", d$words)
	d$height[tailed] = d$height[tailed] * 1.3
	
	# Rotate words
	if(any(d$angle == 90)){
		a = d$width[d$angle == 90]
		d$width[d$angle == 90] = convertUnit(unit(d$height[d$angle == 90], "npc"), "npc", axisFrom = "y", axisTo = "x", valueOnly = T)
		d$height[d$angle == 90] = convertUnit(unit(a, "npc"), "npc", axisFrom = "x", axisTo = "y", valueOnly = T)
	}
	
	# Calculate area of all words and apply scaling factor to sizes
	area = sum(d$height * d$width * d$size ** 2)
	d$size = d$size / sqrt(area / scale)
	
	# If some words are bigger than window scale down all words to fit there
	maxsize = max(c(d$height * d$size, d$width * d$size))
	if(maxsize > 1){
		d$size = d$size / maxsize * 0.8
	}
	
	# Resize the words
	d$width = d$width * d$size
	d$height = d$height * d$size
	
	
	# Calculate the coordinates of words
	if(tryfit) {
		dontfit = T
		while(dontfit){
			if(algorithm == "leftside"){
				a = findCoordinates_left(d$width, d$height)
			}
			if(algorithm == "rightside"){
				a = findCoordinates_left(d$width, d$height)
				a[, 1] = 1 - a[, 1]
			}
			if(algorithm == "circle"){
				a = findCoordinates(d$width, d$height)
			}
			
			if(any(a[, 2] == 3)){
				d$size = d$size * 0.95
				d$width = d$width * 0.95
				d$height = d$height * 0.95
			}
			else{
				dontfit = F
			}
		}	
	}
	else{
		if(algorithm == "leftside"){
			a = findCoordinates_left(d$width, d$height)
		}
		if(algorithm == "rightside"){
			a = findCoordinates_left(d$width, d$height)
			a[, 1] = 1 - a[, 1]
		}
		if(algorithm == "circle"){
			a = findCoordinates(d$width, d$height)
		}
	}

	d$x = a[, 1]
	d$y = a[, 2]
	
	# Issue warning when something is ńot drawn
	if(any(d$y == 3)){
		warning(sprintf("Words not drawn: %s", paste(as.character(d$words[d$y == 3]), collapse = ", ")))
	}
	
	d = d[!(d$y == 3), ]
	
	if(nrow(d) < 1){
		warning("No words to be drawn. Try adjusting scales if some of them did not fit.")
		return()
	}
	
	# Shift the talied words upwards
	tailed = grepl("g|j|p|q|y|_", d$words)
	if(any(d$angle == 90)){
		if(algorithm == "circle"){
			d[tailed & d$angle == 90, "x"] = d[tailed & d$angle == 90, "x"] - 0.1153846 * d[tailed & d$angle == 90, "width"]
		}
	}
	d[tailed & d$angle == 0, "y"] = d[tailed & d$angle == 0, "y"] + 0.1153846 * d[tailed & d$angle == 0, "height"]
	
	
	# Draw the words
	if(algorithm == "leftside"){
		hjust = ifelse(d$angle == 90, 0.5, 0)
		vjust = ifelse(d$angle == 90, 1, 0.5)
		grid.text(d$words, d$x, d$y, rot = d$angle,  hjust = hjust, vjust = vjust, gp = gpar(cex = d$size, col = d$colors))
	}
	if(algorithm == "rightside"){
		hjust = ifelse(d$angle == 90, 0.5, 1)
		vjust = ifelse(d$angle == 90, 0, 0.5)
		grid.text(d$words, d$x, d$y, rot = d$angle, hjust = hjust, vjust = vjust, gp = gpar(cex = d$size, col = d$colors))
	}
	if(algorithm == "circle"){
		grid.text(d$words, d$x, d$y, rot = d$angle, gp = gpar(cex = d$size, col = d$colors))
	}
}
