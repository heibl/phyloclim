## This code is part of the phyloclim package
## © C. Heibl 2009 (last update 2020-01-18)

#' @title  Plot Predicted Niche Occupancy Profiles
#' @description Plot predicted niche occupancy profiles (PNOs). PNOs can be
#'   obtained in a geographical information system by summing the cumulative
#'   probabilies of each climatical value for a species distribution model
#'   (SDM).
#' @param x A data frame or matrix with columns corresponding to species and
#'   rows corresponding to values along an environmental gradient. The first
#'   columns contains the environmental variable, the remaining colums
#'   probabilities of suitability.
#' @param subset A vector of mode \code{"character"} which can be used to
#'   restrict the calculation of weighted means to those columsn in \code{x}
#'   whose column names match \code{subset}; defaults to \code{NULL}.
#' @param thinning An integer that can be used to thin fuzzy PNOs prior to
#'   plotting; defaults to \code{NULL}.
#' @param xlab A character string given the label for the x-axis.
#' @param tail_threshold A numeric that can be used cut long tails of PNOs;
#'   defaults to \code{0}.
#' @param wm Logical indicating if weighted mean should added for each species.
#' @param legend.pos Controls the position of the legend. Might eihter be a list
#'   object containing x and y coordinates (such as e.g. returned by
#'   \code{\link{locator}}) of the \bold{topleft corner} of the legend box or
#'   one of the following: \code{"topleft"} (default), \code{"bottomleft"},
#'   \code{"topright"}, or \code{"bottomright"}. If \code{legend.pos == NULL}
#'   the plotting of the legend is suppressed.
#' @param legend.cex Numeric controlling the size of the legend.
#' @references \insertRef{evanssmith2009}{phyloclim}
#' @seealso \code{\link{pno}}
#' @examples
#' # load PNOs for Oxalis sect. Palmatifoliae
#' data(PNO)
#'
#' # plot predicted niche occupany for annual mean temperature
#' plotPNO(x = PNO$AnnualMeanTemperature,
#'         xlab = "Annual Mean Temperature (degree C)")
#'
#' # same plot, but with weighted means added
#' plotPNO(x = PNO$AnnualMeanTemperature,
#'         xlab = "Annual Mean Temperature (degree C)", wm = TRUE)
#' @importFrom grDevices rainbow
#' @importFrom graphics legend lines plot
#' @export

plotPNO <- function(x, subset = NULL, thinning = NULL, xlab = NULL, 
                    tail_threshold = 0, wm = FALSE, legend.pos = "topleft", legend.cex = 1){
	
	# calculate weighted means:
	# -------------------------
	wmean <- pno.weighted.mean(x, subset = subset)
		
	# subset matrix
	# ---------------------
	if (!is.null(subset))
		x <- x[, c(1, which(names(x) %in% subset))]
		
	for (i in 2:(dim(x)[2])){
		nf <- sum(x[, i])
		x[, i] <- x[, i] / nf
	}
	
	# delete 'zero tails'
	# ---------------------
	if (dim(x)[2] > 2)
	    zeros <- which(apply(x[, -1], 1, sum) <= tail_threshold)
	else 
	    zeros <- which(x[, -1] <= tail_threshold)
	if (length(zeros) > 0)
			x <- x[-zeros, ]
			
	# thin matrix
	# --------------------
	if (!is.null(thinning)) 
		x <- x[seq(1, dim(x)[1], length.out = thinning), ]
	
	
	col <- rainbow(dim(x)[2] - 1)
	max_val  <-  max(x[, -1])
	plot(x[, 1], x[, 2], type = "l", col = col[1], 
	    ylim = c(0, 	max_val),
		main = "Predicted niche occupancy",
		xlab = xlab,
		ylab = ""
		)
		
	if (dim(x)[2] > 2)
	    for (i in 3:(dim(x)[2]))
		    lines(x[, 1], x[, i], col = col[i - 1])
		
	# plot legend:
	# --------------------
	if (!is.null(legend.pos)){
		if (is.list(legend.pos))
		    legend(x = legend.pos$x, y = legend.pos$y, legend = colnames(x)[2:dim(x)[2]], fill = col,
		           cex = legend.cex)
		else {
			if (legend.pos == "topleft")
			    lxy <- legend(x = min(x[, 1]), y = max(x[, -1]), 
			                  legend = colnames(x)[2:dim(x)[2]], fill = col, cex = legend.cex)
			
			if (legend.pos == "bottomleft"){
				lxy <- legend(x = min(x[, 1]), y = min(x[, -1]), legend = colnames(x)[2:dim(x)[2]], 
				    fill = col, plot = FALSE)$rect
				xx <- min(x[, 1])
				yy <- min(x[, -1]) + lxy$h
			    legend(x = xx, y = yy, fill = col, legend = colnames(x)[2:dim(x)[2]], cex = legend.cex)
			}
			if (legend.pos == "topright"){
				lxy <- legend(x = min(x[, 1]), y = min(x[, -1]), legend = colnames(x)[2:dim(x)[2]], 
				    fill = col, plot = FALSE)$rect
				xx <- max(x[, 1]) - lxy$w
				yy <- max(x[, -1]) 
			    legend(x = xx, y = yy, fill = col, legend = colnames(x)[2:dim(x)[2]], cex = legend.cex)
			}
			if (legend.pos == "bottomright"){
				lxy <- legend(x = min(x[, 1]), y = min(x[, -1]), legend = colnames(x)[2:dim(x)[2]], 
				    fill = col, plot = FALSE)$rect
				xx <- max(x[, 1]) - lxy$w
				yy <- min(x[, -1]) + lxy$h
			    legend(x = xx, y = yy, fill = col, legend = colnames(x)[2:dim(x)[2]], cex = legend.cex)
			}
		}
	}
	
		
	# plot weighted means:
	# --------------------
	if (wm){
		for (i in seq(along = wmean))
			lines(rep(wmean[i], 2), range(x[, -1]), col = col[i], lty = 3, lwd = 3)
	}
}