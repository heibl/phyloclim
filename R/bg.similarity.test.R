## This code is part of the phyloclim package
## © C. Heibl 2009 (last update 2020-01-15)

#' @rdname niche.tests
#' @importFrom methods slot
#' @importFrom raster extract unstack mask sampleRandom writeRaster
#' @importFrom rgeos gConvexHull
#' @importFrom sp coordinates read.asciigrid SpatialPointsDataFrame
#' @importFrom utils write.table
#' @export

bg.similarity.test <- function(p, env, n = 99, study.area.y = "mcp", conf.level = .95, app, dir){
  
  ## Checks and definitions
  ## ----------------------
  study.area.y <- match.arg(study.area.y, c("env", "mcp"))
  
  ## Prepare file system 
  ## -------------------
  DIR <- ifelse(missing(dir), tempdir(), dir)
  if (dir.exists(DIR))
    unlink(DIR, recursive = TRUE)
  dir.create(DIR)
  dir.create(ODIR <- file.path(DIR, "out/"))
  dir.create(PDIR <- file.path(DIR, "proj"))
  
  p <- SpatialPointsDataFrame(coords = p[, 2:3], data = p)
  names(p) <- c("species", "long", "lat")
  p <- p[order(p$species), ]
  n_occ <- table(as.character(p$species)) # corresponds to 'o' (p.2872)
  if (length(n_occ) != 2) stop("'p' must contain exactly two species")
  
  ## Append covariate data to presence points
  ## ----------------------------------------
  attributes <- extract(x = env, y = p)
  p <- cbind(p, attributes)
  NAs <- which(is.na(slot(p, "data")), arr.ind = TRUE)
  NAs <- unique(NAs[, 1])
  if (length(NAs)){
    p <- p[-NAs, ]
    warning(length(NAs), " presence points with missing environmental data removed")
  }
  
  ## Sample n(spec1) and n(spec2) points from the study area, do it n times
  ## ----------------------------------------------------------------------
  spec_vect <- as.character(p$species)
  randomPresence <- function(i, n_occ, env1, env2, spec_vect){
    name <- paste(spec_vect, i, sep = "_")
    s <- rbind(sampleRandom(env1, size = n_occ[1], sp = TRUE),
               sampleRandom(env2, size = n_occ[2], sp = TRUE))
    SpatialPointsDataFrame(coords = s, data = data.frame(name, coordinates(s), slot(s, "data")))
  }
  
  if (study.area.y == "mcp"){
    
    ## Define the "study areas" as convex hull of the presence points
    S1 <- gConvexHull(p[p$species == names(n_occ)[1], ], byid = FALSE, id = NULL)
    env1 <- mask(env, S1)
    S2 <- gConvexHull(p[p$species == names(n_occ)[2], ], byid = FALSE, id = NULL)
    env2 <- mask(env, S2)
    rp <- lapply(1:n, FUN = randomPresence, n_occ = n_occ, 
                 env1 = env1, env2 = env2, spec_vect = spec_vect)
  } else {
    
    ## study.area.y == "env"
    rp <- lapply(1:n, FUN = randomPresence, n_occ = n_occ, 
                 env1 = env, env2 = env, spec_vect = spec_vect)
  }
  rp <- do.call(rbind, rp)
  names(rp)[1:3] <- c("species", "long", "lat")
  
  ## Produce plot to check observed and sampled presence points
  ## ----------------------------------------------------------
  # pdf(file.path(DIR, "presence-points.pdf"))
  # plot(env[[1]])
  # plot(rp[grep(names(n_occ)[1], rp$species), ], pch = 21, add = TRUE, col = "blue")
  # plot(rp[grep(names(n_occ)[2], rp$species), ], pch = 21, add = TRUE, col = "red")
  # plot(p[p$species == names(n_occ)[1], ], pch = 21, add = TRUE, bg = "skyblue")
  # plot(p[p$species == names(n_occ)[2], ], pch = 21, add = TRUE, bg = "orange")
  # plot(S1, add = TRUE, border = "blue")
  # plot(S2, add = TRUE, border = "red")
  # dev.off()
  
  ## Bind observed and random presence points together
  ## -------------------------------------------------
  p <- rbind(p, rp)
 
  ## Sample 99999 background points from env
  ## ---------------------------------------
  bg <- sampleRandom(env, size = 9999, na.rm = TRUE, sp = TRUE)
  bg <- data.frame("background", coordinates(bg), slot(bg, "data"))
  
  
  ## Save input files
  ## ----------------
  write.table(bg, file.path(DIR, "background.csv"), 
              row.names = FALSE, col.names = TRUE, sep = ",")
  write.table(slot(p, "data"), file.path(DIR, "samples.csv"),
              row.names = FALSE, col.names = TRUE, sep = ",")
  fn <- file.path(PDIR, paste(names(env), "asc", sep = "."))
  env <- unstack(env)
  for (i in seq_along(fn)){
    writeRaster(x = env[[i]], filename = fn[i], format = "ascii", 
                overwrite = TRUE, NAflag = -9999)
  }
	
  ## Call MAXENT:
  ## ------------
  togglelayertype <- ifelse(length(grep("cat_", names(env))) > 0, "-t cat_", "")
  CALL <- paste("java -jar", app,   	
                "-e ", file.path(DIR, "background.csv"),
                "-s ", file.path(DIR, "samples.csv"),
                "-j ", PDIR, 
                "-o ", ODIR, 	
                togglelayertype,
                "-r removeduplicates nopictures autorun")
  system(CALL, wait = TRUE)
	
  ## Calculate D and I for observed data
  ## -----------------------------------
  fns <- paste0(ODIR, names(n_occ), "_proj.asc")
  x <- read.asciigrid(fns[1])
  y <- read.asciigrid(fns[2])
  di <- di.enm(x = x, y = y)
  
  # Calculate D and I for null distributions 
  # ----------------------------------------
  di.x.randomY <- sapply(X = 1:n, FUN = di.enm, x = x, y = fns[2])
  di.x.randomY <- t(di.x.randomY)
  di.y.randomX <- sapply(X = 1:n, FUN = di.enm, x = fns[1], y = y)
  di.y.randomX <- t(di.y.randomX)
  
  ## CIs for null distributions
  ## --------------------------
  conf.limits <- c((1 - conf.level) / 2, 1- (1 - conf.level) / 2)
  ci.x.randomY <- apply(di.x.randomY, 2, quantile, probs = conf.limits)
  ci.y.randomX <- apply(di.y.randomX, 2, quantile, probs = conf.limits)
  
#   h0.x.randomY <- ci.x.randomY[1, ] < di & di < ci.x.randomY[2, ]
#   h0.y.randomX <- ci.y.randomX[1, ] < di & di < ci.y.randomX[2, ]
  
  # The null hypothesis that measured niche overlap between species is explained
  # by regional similarities or differences in available habitat (and not by
  # niche conservatism), is rejected if the actual similarity between two
  # species falls outside of the 95% confidence limits of the null distribution.

  # Remove MAXENT output:
  # ---------------------
  if (missing(dir)) unlink(DIR, recursive = TRUE)
	
	# Create output object:
	# ---------------------
  out <- list(
    method = "background similarity test",
    species = names(n_occ),
    null = paste0("niche models are either more similar \n", 
                 paste(rep(" ", 24), collapse = ""), 
                 "or more different than expected by chance"),
    statistic = di,
    ci.x.randomY = ci.x.randomY,
    ci.y.randomX = ci.y.randomX,
    nd.x.randomY = di.x.randomY,
    nd.y.randomX = di.y.randomX
    )
  class(out) <- "ntest"
  out
}
