require(methods)
require(xts)
require(raster)

setOldClass("spatial")
setOldClass("xts")

# creating the DTM project for a subcatchment
setClass("Catchment",
         list("name"="character",
           "dem"="Raster",
           "drn"="spatial",
           "id"="numeric",
           "name"="character",
           "type"="character",
           "parent"="character",
           "downstream"="character",
           "gauging_stations"="list",
           "aws"="list",
           "props"="list"  # anything else such as channel width
           
           )
)

setMethod("plot", c(x="Catchment", y="missing"),
          function(x, y, ...)
          {
            plot(x@dem, ...)
            plot(x@drn, col="blue")
          }
)

setClass("Gauge",
         list("location"="spatial",
           "id"="numeric",
           "name"="character",
           "code"="character",
           "props"="list"   # anything else e.g rating_curve
         )
)         


# e.g. a rainfall record that can be associated with a gauge
setClass("Observations",
         list("dt"="numeric",
              "values"="xts",
              "id"="numeric",
              "units"="character")
)

# parse a directory data into a subcatchment
create_catchment <- function(dn, cb, dem, drn, type, parent="", downstream="")
{
  message("Creating ", cb$wb_name[1], "...")
  fn <- file.path(dn, "dem.tif")
  if(file.exists(fn))
  {
    dem <- raster(fn)
  }
  else
  {
    message("Masking DEM to catchment boundaries...")
    dem <- mask(dem, cb)
  }
  
  if(shapefile_exists(dn, "drn"))
  {
    drn <- readOGR(dn, "drn")
  }
  else
  {
    message("Masking DRN to catchment boundaries...")
    drn <- gIntersection(drn, cb, byid=T)
  }

  catchment <- new("Catchment",
                   type=cb$SubCatchTy[1],
                   name=name,
                   id=cb$OBJECTID[1],
                   parent=parent,
                   drn=drn,
                   dem=dem+0,     # ensure raster is held in memory or the object may be invalid on other machines
                   downstream=downstream
  )
  
  return(catchment)
}