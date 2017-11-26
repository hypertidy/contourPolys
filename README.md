# contourPolys


polygons from contourLines


```R
library(raster)
r <- aggregate(raster(volcano), fact = 2, fun = median) %/% 20
r <- tabularaster:::set_indextent(r)
plot(r, asp = "")
#val <- cellStats(r, min)
#r <- extend(r, 1, value = val-1)

x <- list(x = xFromCol(r), y = rev(yFromRow(r)), z = t(as.matrix(flip(r, "y"))))
cl <- contourLines(x, levels = sort(na.omit(unique(values(r)))))
bound <- as(as(extent(r) - res(r)/2, "SpatialPolygons"), "SpatialLines")@lines[[1]]@Lines

boundary_coords <- rbind(cbind(xmin(r), yFromRow(r)), 
                         cbind(xmax(r), yFromRow(r)), 
                         cbind(xFromCol(r), ymin(r)), 
                         cbind(xFromCol(r), ymax(r)))

nearest_point <- function(coords, pt) {
  distances <- sqrt((coords[,1] - pt[1])^2 + 
                      (coords[,2] - pt[2])^2)
  
  coords[which.min(distances), , drop = FALSE]
  
}
for (i in seq_along(cl)) {
  xxs <- cl[[i]]$x
  yys <- cl[[i]]$y
  
  ## create polygon if bounded
  if (bounded) {
    
  } else {
    ## otherwise insert extra coordinates to nearest side
    pt1 <- nearest_point(boundary_coords, cbind(xxs[1], yys[1]))
    pt2 <- nearest_point(boundary_coords, cbind(xxs[1], yys[1]))
  }
}


## get all start coordinates 
points(
  do.call(rbind, lapply(cl, function(x) cbind(rev(x$x)[1], rev(x$y)[1])))
)
#cl <- contourLines(x, levels = sort(setdiff(na.omit(unique(values(r))), val )))


plot(r)
purrr::walk(purrr::map(cl, function(x) x[c("x", "y")]), lines)
l <- sp::SpatialLines(list(sp::Lines(c(bound, lapply(cl, function(x) sp::Line(cbind(x$x, x$y)))), "1")))

plot(r)

plot(l, add = T)
```
