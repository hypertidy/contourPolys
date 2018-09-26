topo <- crop(SOmap:::Bathy, new("Extent", xmin = 1854636.90418169, xmax = 6300684.27722, 
                                ymin = -736436.317086573, ymax = 1840242.32060924))

usethis::use_data(topo)