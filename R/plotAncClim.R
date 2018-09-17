## This code is part of the phyloclim package
## © C. Heibl 2009 (last update 2018-05-23)

#'@title Chronograms with Climatic Data on the Y-Axis
#'@description Plot the history of climatic tolerance for a clade \emph{sensu}
#'  Evans et al. (2009).
#'@param x A list with elements \code{tree}, \code{data}, and (optional)
#'  \code{central.density} (see details).
#'@param layer \emph{Do not use.}
#'@param clades A list containing vectors of tip labels which define  the clades
#'  to highlight.
#'@param col A vector containing colors for different clades (see Examples).
#'@param density Logical, if \code{TRUE}, the central density intervals for
#'  recent species are plotted.
#'@param tipmode Integer: \code{tipmode = 0} means no tiplabels (and no central
#'  density intervals) are plotted; values \code{1}, \code{2}, and \code{3}
#'  define different ways of plotting the tiplabels (see Examples).
#'@param nchar An integer giving the number of characters to which the tiplabels
#'  are truncated.
#'@param cex Numeric \bold{c}haracter \bold{ex}pansion factor for tiplabels;
#'  multiplied by \code{\link{par}}(\code{"cex"}) yields the final character
#'  size. \code{NULL} is are equivalent to \code{1.0}.
#'@param tipspace The fraction of the x-axis that is reserved for tiplabel
#'  plotting. If no value is given (default), \code{plotAncClim} calculates the
#'  fraction as \emph{1 - (4 / nbtips)}, but this can be overridden by
#'  specifiying a fixed value for \code{tipspace} in \emph{]0, 1[}.
#'@param cladespace A positive reel number; \bold{tentative}: a factor
#'  controlling the space between tiplabels of different clades.
#'@param lwd The line width, a positive number, defaulting to 1.
#'@param ylab A character string, giving a label for the y-axis, i.e., for the
#'  bioclimatic dimension of the plot.
#'@details The main argument \code{x} is a list consisting of at least the first
#'  two of the following elements: (1) \code{tree} is an ultrametric
#'  phylogenetic tree stored as object of class \code{phylo}; (2) \code{data} is
#'  an object of class \code{matrix}; its columns correspond to bioclimatic
#'  variables and its rows corresond to node numbers such as defined by class
#'  \code{phylo} objects; (3) \code{central.density} must only be included if
#'  \code{density = TRUE} -- it is a list, which for every bioclimatic variable,
#'  contains a matrix that contains the some minimum and maximum quantile of the
#'  respective bioclimatic variable for every tip in the tree.
#'@references  Evans, M. E. K., S. A. Smith, R. S. Flynn, and M. J. Donoghue.
#'  2009. Climate, niche evolution, and diversification of the 'bird-cage
#'  evening primroses' (\emph{Oenothera}, sections \emph{Anogra} and
#'  \emph{Kleinia}). \emph{Am. Nat.} \bold{173}: 225-240.
#'@seealso \code{\link{pno}}, \code{\link{pno.weighted.mean}},
#'  \code{\link{anc.clim}}
#' @examples
#'  # load phylogeny and PNOs of Oxalis sect. Palmatifoliae
#'  data(tree)
#'  data(PNO)
#'
#'# choose summer precipitation for analysis
#'clim <- PNO$PrecipitationWarmestQuarter
#'
#'# estimate ancestral tolerances
#'ac <- anc.clim(target = tree, pno = clim, n = 100)
#'
#'# visualize results with default branch coloration
#'plotAncClim(ac)
#'
#'# alternative clade colors are given according to the order
#'# in which tip labels appear from left to right
#'plotAncClim(ac, col = c("red", "purple", "blue"))
#'
#'# the 'tipmode' argument
#'plotAncClim(ac, tipmode = 0)
#'plotAncClim(ac, tipmode = 1)
#'plotAncClim(ac, tipmode = 2, nchar = 5)
#'plotAncClim(ac, tipmode = 3, nchar = 4)
#'@importFrom grDevices rainbow
#'@importFrom graphics lines plot points strheight strwidth text
#'@export

plotAncClim <- function(x, layer, clades = NULL, col, density = TRUE, 
     tipmode = 1, nchar = 3, cex, tipspace, cladespace = 1, 
     lwd, ylab = ""){
	
	# Get data
  # Not sure how the layer argument was used
	# ----------------------------------------
	tr <- x$tree
	tips <- tr$tip.label
	nbtips <- length(tips)
	if (missing(layer)) {
	  clim <- x$means
	} else {
	  clim <- x$means[, layer]
	}
	
	if (density)
	  if (missing(layer)) {
	    cd <- x$central.density
	  } else {
	    cd <- x$central.density[[layer]]
	  }
		
	if ( missing (cex) ) cex <- 1
	if ( missing (lwd) ) lwd <- 1
		
	# if no clades are given use terminal sisters 
	# for coloration instead. Species in a grade
	# will receive their own color	
	if ( is.null(clades) ) {
		ts <- terminal.sisters(tr)
		nts <- tr$tip.label[!tr$tip.label %in% ts]
		for (i in seq(along = nts))
			ts <- rbind(ts, rep(nts[i], 2))
		lts <- vector(mode = "list", length = dim(ts)[1])
		for (i in seq(along = lts))
			lts[[i]] <- ts[i, ]
		clades <- lts
	}
	
	# abbreviate taxon labels:
	# -------------------------
	tips <- substring(tr$tip, 1, nchar) 
		
	# calculate x-coordinate for nodes bases on node ages
	# and 'max_age'
	# -------------
	nodeages <- c(rep(0, nbtips), -branching.times(tr))
	max_age <- -max(branching.times(tr))
	
	# coordinates for plotting: node ages versus clima values
	# --------------------------------------------------------
	xy <- cbind(nodeages, clim)
	
	# plot coordinate system
	# ----------------------
	if ( !density ){
	  yrange <- range(xy[, 2])
	} else {
	  yrange <- c(min(cd[1, ]), max(cd[2, ]))
	}		 
	if (tipmode == 3) yrange[2] <- yrange[2] + 0.05 * diff(yrange)
	if (tipmode == 0) tipspace <- 0
	if (missing(tipspace)){
		tipspace <- 1 - (4 / nbtips) 	
	}
	tipspace <- -tipspace * max_age
	plot(x = c(max_age, tipspace), y = yrange,
		type = "n",
		xlab = "Time (Ma)", ylab = ylab,
		bty = "c", xaxp = c(floor(max_age), 0, 5))
	
	# calculate space needed for tip labels:
	# -------------------------------------
	tipxpos <- max(strwidth(tips, units = "f", cex = cex)) * 1:nbtips
	tipxpos <- tipxpos * cladespace
	
	# calculate coordinates for tiplabels
	# -----------------------------------
	tex <- cbind(tipxpos, xy[1:nbtips, 2])
	rownames(tex) <-  tr$tip.label
	for (i in 2:length(clades)){
		id <- which(rownames(tex) %in% clades[[i]])
		extraspace <- (-tex[1, 1]  + tex[1, 1] * 1) * (i - 1)
	    tex[id, 1] <- tex[id, 1] + extraspace
	}
	rownames(tex) <-  tips
	maxtip <- max(tex[, 1])
	tex[, 1] <-  tex[, 1] / maxtip * tipspace
	totalspace <- max(tex[, 1])
			
	# edge colors:
	# ------------	
	n <- noi(tr, clades)
	n <- lapply(n, descendants, tree = tr, internal = TRUE)
	lincol <- rep("grey", dim(xy)[1])
	if (missing(col))
		col <- rainbow(length(n))
	for (i in seq(along = n)){
		lincol[n[[i]]] <- col[i]
	}
	
	# plot edges:
	# ----------------
	colind <- vector() # used to order tip colors!
	for (i in seq(along = tr$edge[, 1])){
		ind <- tr$edge[i, ]
		lines(xy[ind, 1], xy[ind, 2], lwd = lwd, col = lincol[ind[2]])
		colind <- c(colind, ind[2])
	}
	
	# plot zero line:
	# ----------------
	if (min(xy[, 2]) < 0 & max(xy[, 2] > 0))
	lines(x = c(min(xy[, 1]), 0), y = rep(0, 2), lwd = 0.8 * lwd, 
	    lty = "12", col = "gray50")
	
	# density:
	# ----------------
	if (density & tipmode != 0){
		for (i in seq(along = cd[1, ]))
		lines(rep(tex[i, 1], 2), cd[, i], 				
		    col = lincol[1:nbtips][i], lwd = lwd, lty = "12")
	}
	
	# plot tiplabels:
	# ----------------
	if (tipmode != 0)
	    if (tipmode == 1)
	        text(tex[, 1], tex[, 2], rownames(tex), cex = cex, 		
	            adj = 0.5, col = lincol[1:nbtips], font = 4)
	    
	    if (tipmode == 2){
	    	points(tex, col = lincol[1:nbtips], pch = 18)
	    	id <- mean(yrange)
	    	id <- tex[, 2] < id
	    	offset <- strheight(".", cex = cex, font = 4)
	    	text(tex[id, 1], cd[2,id] + offset, rownames(tex)[id], 
	    	    cex = cex, adj = c(0, 0.5), font = 4, srt = 90, 
	    	    col = lincol[1:nbtips][id])
	        text(tex[!id, 1], cd[1, !id] - offset, rownames(tex)[!id], 
	            cex = cex, adj = c(1, 0.5), srt = 90, 
	            col = lincol[1:nbtips][!id], font = 4)
	    }
	    if (tipmode == 3){
	        points(tex, col = lincol[1:nbtips], pch = 18)
	        text(tex[, 1], max(cd[2,]), rownames(tex), cex = cex,
	            adj = c(0, 0.5), col = lincol[1:nbtips], font = 4, 
	            srt = 45)
	    }
}
