## This code is part of the phyloclim package
## © C. Heibl 2009 (last update 2018-05-24)

#' @importFrom ape getMRCA

nbConnectingNodes <- function(phy, npair){
	ntips <- length(phy$tip.label)
	nds <- getMRCA(phy, npair)
	nds <- descendants(phy, nds, internal = TRUE)
	if (identical(sort(nds), sort(npair)))
		nb <- 1										else {
		nds <- nds[nds > ntips]
		check <- function(x, npair)
			any(npair %in% descendants(phy, x))
		id <- sapply(nds, check, npair = npair)
		nds <- nds[id] 
		nb <- length(nds) + 1 
	}
	nb
}