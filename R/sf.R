poly1 <- function(x) {
  nr <- length(x)/2; 
  structure(list(matrix(x, ncol = 2)[c(seq_len(nr), 1), ]), 
            class = c("XY", "POLYGON", "sfg"))
}

fc_interval_sf <- function(m) { 
  xx <- lapply(split(m[, 1:2], rep(m[, 5], 2)), poly1)
  uu <- unlist(lapply(xx, sf::st_is_valid))
  
  x <- lapply(split(xx[uu], m[!duplicated(m[,5]), 3][uu]), 
              sf::st_geometrycollection)
  
  sf::st_union(sf::st_sfc(x))
}


sfill_contour <- function(x, levels = NULL, ..., maxpixels = 2e16) {
  UseMethod("sfill_contour")
}
sfill_contour.BasicRaster <- function(x, levels = NULL, ..., maxpixels = 2^15) {
  if (ncell(x) > maxpixels) {
    warning("setting maxpixels", immediate. = TRUE)
    x <- raster::sampleRegular(x, size = maxpixels, asRaster = TRUE)
  }
  if (is.null(levels)) levels <- pretty(c(raster::cellStats(x, min), 
                                          raster::cellStats(x, max)))
  
  z <- t(raster::as.matrix(x)[nrow(x):1, ])
  y <- rev(yFromRow(x))
  x <- xFromCol(x)
  sfgs <- vector("list", length(levels) - 1)
  for (ilevel in seq_along(sfgs)) {
    p <- fcontour(x, y, z, levels[c(ilevel, ilevel + 1L)])

    m <- cbind(x = unlist(p[[1]]), 
           y = unlist(p[[2]]), 
           lower = rep(unlist(p[[3]]), lengths(p[[1]])), 
           upper = rep(unlist(p[[4]]), lengths(p[[1]])), 
           g = rep(seq_along(p[[1]]), lengths(p[[1]]))) 

    sfgs[[ilevel]] <- fc_interval_sf(m)
    print(sprintf("level %i of %i", ilevel, length(sfgs)))
  }
  sfc_vec <- sf::st_set_crs(sf::st_cast(do.call(c, sfgs)), raster::projection(x))
  sf::st_sf(geometry = sfc_vec, min = head(levels, -1L), max = tail(levels, -1L))
}
