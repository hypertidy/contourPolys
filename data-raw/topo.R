topo <- crop(SOmap:::Bathy, new("Extent", xmin = 1854636.90418169, xmax = 6300684.27722, 
                                ymin = -736436.317086573, ymax = 1840242.32060924))

usethis::use_data(topo)

library(contourPolys)
library(raster)
library(sf)
levels <- c(-6000, -4000, -2000, 0, 2000, 4000)

fc <- contourPolys::fcontour_sf(xFromCol(topo), rev(yFromRow(topo)), t(as.matrix(flip(topo, "y"))), c = levels)
library(sf)
fragments_sf <- st_sf(st_sfc(fc[[1]], crs = projection(topo)), region = unlist(fc[[2]]))
usethis::use_data(fragments_sf)
# saveRDS(fragments_sf, file = "fragments_sf.rds")
# system("zip fragments_sf.zip fragments_sf.rds")
