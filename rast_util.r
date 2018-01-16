# set of useful matrix <-> raster operatyions
# convert the data in the raster object into a list of simple vars and matrix
require(raster)
require(rgdal)
#load.source("util.r", chdir=T)
#load.source("sp_util.r")

# # scale up the raster by the given factor, keep origin or centre, dims and
# values
expand.raster <- function(rast, fact, keep.origin = T, crs=rast@crs)
{
   vals <- rast[]

   ext <- extent(rast)
   # scales the extent around original centre
   ext2 <- ext*fact

   rast2 <- raster(ncol=ncol(rast), nrow=nrow(rast),ext=ext2,
                   #xmn=ext2@xmin, xmx=ext2@xmax, ymn=ext2@ymin, ymx=ext2@ymax,
                   crs=crs)

   # reestablish origin
   if(keep.origin)
   {
     rast2<-shift(rast2, x=ext@xmin-ext2@xmin, y=ext@ymin-ext2@ymin)
   }
   rast2 <- raster::setValues(rast2, getValues(rast))
   return(rast2)
}

# write the given raster to a disk file, overwriting silently
overwrite.raster <- function(rast, fn)
{
  res <- rast
  overwrite <- file.exists(fn)
  if(is.raster.stack(rast))
  {
    #		res <- writeRaster(rast, fn)
    if(!file.exists(fn))
    {
      res <- writeRaster(rast, fn, overwrite=T)
    }
    # 		else
    # 		{
    #
    # 		}
  }
  return(res)
}


# encoding a multiband raster into a single layer whose cells
# are the sum of layer values multiplied by successive powers of two?
encode_multiband <- function(rast)
{
	vals <- rast
	for(i in 1:nlayers(rast))
	{
		vals[[i]] <- vals[[i]]*10^i
	}

	return(sum(vals))
}

decode_multiband <- function(rast)
{
	
}

# similar to addLayer but allows the new layer names to be specified on creattion
add.named.layers <- function(x, ...)
{
	layers <- as.list(...)
	nms <- c(names(x), names(layer))

	for(i in 1:length(layers))
	{
		x <- addLayer(x, layers[[i]])
	}

	names(x)<- nms

}

# aggregate the given raster by the specified factor if > 1 and rast not null
aggregate.null <- function(rast, fact, fun=mean)
{
	fact <- round(fact)
	if(is.raster.stack(rast) & fact > 1)
	{
		rast <- aggregate(rast, fact, fun=fun)
	}
	return(rast)
}

# existance test
is.raster.stack <- function(obj)
{
	if(is.null(obj)){return(FALSE)}
	if(inherits(obj, "RasterStack") | inherits(obj, "RasterLayer") | is(obj, "RasterBrick")){return(TRUE)}
	return(FALSE)
}



load.rast.if.exists <- function(fn, dir=NULL)
{
	if(!is.null(dir))
	{fn <- file.path(dir, fn)}
	if(file.exists(fn)){try(return(raster(fn)))}
	return(NULL)
}

merge.rasts.dir <- function(dn)
{

	rast.fn <- dir(dn, "*.tif", recursive=T, full.names=T)
	rasts <- lapply(rast.fn,
									function(fn)
									{
										rast <- NULL
										try(rast <- raster(fn))
										return(rast)
									})
	rasts <- delete.NULLs(rasts)
	merged <- merge.rasts(rasts)
	return(merged)
}




# iterate through every raster located in given directory and merge with previous
# results.
merge_rasts <- function(rasts=NULL, maxcell=1e8)
{
	if(is.null(rasts)){return(NULL)}
	if(!is.list(rasts)){rasts <- list(rasts)}

	ints <- 0:length(rasts)*maxcell
	celltots <- sapply(rasts, ncell)
	group.nums <- findInterval(cumsum(celltots), ints)
	group.nums <- split(1:length(rasts), group.nums)
	
  res <- lapply(group.nums,

        function(irasts)
        {
            cat("Merging ", length(irasts), "...\n")
            merged <<- NULL

            lapply(rasts[irasts],
               function(rast)
               {
                   if(is.null(merged))
                   {
                       merged <<- rast
                   }
                   else
                   {
                       try(merged <<- raster::merge(merged, rast))
                   }
               }
            )
            return(merged)
        }
    )
  rm(merged)
  gc()

	return(res)
}

# raster may be refer to disk file which can have a significant performance impact
# return a raster (stack) with values of each layer  provided copied to memory
rast.mem.copy <- function(rast)
{
  if(raster::inMemory(rast)){return(rast)}
#  rast.copy <- stack(rast)
  layers <-
    lapply(names(rast),
         function(ln)
         {
            layer <- rast[[ln]]
            raster::setValues(layer, getValues(layer))
         }
  )
  if(length(layers)>1)
  {
    return(stack(layers))
  }
  else
  {
    return(layers[[1]])
  }
}

DEMFromRaster <- function(rast,res=NULL,title=NULL)
{
	# if input already a matrix then return (no info available about res in this case so assume this is known elsewhere)
	if(is.matrix(rast))
	{
		res <- res
		title <- title
		dem <- rast
		# BL is assummed to be 0,0 but we can calculate the extent from the matrix dims and grid size
		extent <- extent(0, (ncol(rast)-1)*res,
						 0, (nrow(rast)-1)*res)
	}
	else
	{
		# otherwise need to convert the raster data to matrix and flip and transpose to match
		# the R convention for delaing with matrix "rasters"

		title <- rast@title
		# assuming here that the x and y resolutions are equal
		res <- xres(rast)
		# needs some manipulation to get the data in the right order
		# i.e. raster treats columns as y, rows as x and proceeds from top to bottom i.e first value
		# is top left of map (NW), the complete opposite to the standard R routines for displaying
		# elevation grids
		# use either rast values or
		dem<- matrix(nrow=rast@ncols, ncol=rast@nrows, getValues(rast))
		extent <- extent(rast)
		# invert vertically
		dem<-dem[,ncol(dem):1]
	}

	return(list("title"=title,"res"=res,"extent"=extent,"dem"=dem))
}

# area not NA within elevaton raster
dem.area <- function(dem)
{

	return(xres(dem)*yres(dem)*length(which(!is.na(dem[]))))
}


# create a raster from gridded data plus georeferencing and other data
# proj4str is the proj4 string specifying the projection to use - use OSGB36by default
# converting tabulated elevation values to a raster DEM
matrix.to.raster.DEM <- function(dem, title="", dx=1, dy=dx,
                                 xmin=0, ymin=0, 
                                 proj4str=NA) #"+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs")
{
  # determine max x and y from the mins and the resolution (assumed equal in both directions)
  ymax = ymin+ncol(dem)*dy
  xmax = xmin+nrow(dem)*dx
	# invert vertically and transpose: raster displays x / y and row / col in opposite sense
	dem<-t(dem[,ncol(dem):1])
  # create and return the raster

  ras <- raster(dem, xmx=xmax, ymx=ymax, xmn=xmin, ymn=ymin,crs=proj4str)
  ras@title=title
	return(ras)
}

#' Upslope contributing area and wetness index calculation
#'
#' @description Determine upslope contributing area based on an elevation raster and, optionally, compute the topographic wetness index.
#' @export upslope.area
#' @import raster
#' @import topmodel
#' @param dem   Elevation raster (in m), using a  projected coordinate system with identical x and y resolutions.
#' @param log   Return the natural log of the values.
#' @param atb   If TRUE, include both the upslope contributing area and the topographic wetness index \eqn{ln(a/tan(beta))}. Otherwise calculate just the upslope area.
#' @param deg   Minimum intercell slope to identify with a sink (degrees).
#' @param fill.sinks Fill sinks before calculation using the threshold angle given by deg.
#' @note This is a wrapper to the function implemented in the TOPMODEL package by Wouter Buytaert.
#' @author Peter Metcalfe and Wouter Buytaert
#' @references Quinn, P. F., Beven, K. J., & Lamb, R. (1995). The In (a/tan/beta) index: How to calculate it and how to use it within the Topmodel framework. Hydrological processes, 9(2), 161-182.
#' @examples
#' require(dynatopmodel)
#' data(brompton)
#'
#' a.atb <- upslope.area(brompton$dem, atb=T)
#' sp::plot(a.atb, main=c("Upslope area (log(m^2/m))", "TWI log(m^2/m)"))
upslope.area <- function(dem, log=T, atb=F, deg=0.1, fill.sinks=T)
{
  # check, but not too close
  if(round(xres(dem))!=round(yres(dem)))
  {
    stop("Raster has differing x and y cell resolutions. Check that it is in a projected coordinate system (e.g. UTM) and use raster::projectRaster to reproject to one if not. Otherwise consider using raster::resample")
  }
  # any sinks still present may give strange results
#  sink(file="e:/junk/sink.txt")
#  on.exit(sink(NULL))
	if(fill.sinks)
	{
	  # use capture.output to supress the function console output
		capture.output(dem <- invisible(raster::setValues(dem, topmodel::sinkfill(raster::as.matrix(dem), res=xres(dem), degree=deg))))
	}
  topidx <- topmodel::topidx(raster::as.matrix(dem), res=xres(dem))

	a <- raster::setValues(dem, topidx$area)
  if(log)
  {
    a <- log(a)
  }
  if(atb)
  {
    atb <- raster::setValues(dem, topidx$atb)
    # add the topographic index ln(a/tanB)
    a <- addLayer(a, atb)
    names(a)<-c("a", "atb")
  }
  return(a)
}

# use function from TOPMODEL
flow.lens <- function(dem,
                      src=NULL,  #  starting cells, defaults to all cells i dem
											agg=1, # initial aggregation factor
											max.agg=4,

											outlet=NULL)  # A vector containing the row and column indices of the pixel representing the catchment outlet. if a single value, treated as the index of a DEM cell
{
		lens <- raster::setValues(dem, NA)
	if(length(outlet)>0)
	{
		outlet.sp <- xyFromCell(dem, outlet, spatial=T)
	}

	dem.agg <- dem
	while(agg <= max.agg & max(c(0,lens[]), na.rm=T)==0)
	{
		if(agg>1)
		{
			message("Trying a aggregated dem to determine flow lengths")
			# try a coarser
			dem.agg <- raster::aggregate(dem, agg)
			#	reaches <- aggregate(reaches, )
		}

		if(length(outlet)>0)
		{
			outlet <- extract.cells(dem.agg, outlet.sp)
			iout<- rowColFromCell(dem.agg, outlet)
		}
    else{iout <- NA}

		dem.agg <- fill.sinks(dem.agg, deg=0.1)
		lens <- raster::setValues(dem.agg, flowlength(as.matrix(dem.agg), outlet=iout))
    if(!is.null(src))
    {
      lens[setdiff(1:ncell(dem), src)]<- NA
    }
		agg <- agg+1
	}
	agg <- agg-1
	# disaggregate
	if(agg>1)
	{
		message("Disaggregating back to original resolution...")
		lens <- raster::disaggregate(lens, agg, method="bilinear")
	}
	return(raster::xres(dem)*lens)
}


# demsrc - original dem
# demaff - dem aggregated by a factor
# return the cells in targ that correspond to the given cell(s) in src, assumming
# that targ is an aggregate (or disaggregate) of src by a integer factor
# MapCells <- function(src, targ, cell)
# {
# 	# check
# 	factc <- round(ncol(src)/ncol(targ),1)
# 	factr <- round(nrow(src)/nrow(targ),1)
# 	fact <- round(sqrt(ncell(src)/ncell(targ)),1)
# 	if(factr != factc){ stop("Source and aggregated rasters have inconsistent dimensions")}
# 	# could check extents at this point also, but aggregation could lead to differences
# 	cells <- sapply(as.list(cell),
# 		   function(c)
# 			{
# 				ext <- CellExtent(src, c)
# 				cellsFromExtent(targ, ext)
# 			}
# 	)
# 	res <- unique(as.vector(cells))
# 	return(res)
# }

MapCell <- function(src, targ, cell)
{
	fact <- round(sqrt(ncell(targ)/ncell(src))) # this assummes src is larger than targ, i.e. has been aggregated from it
	rc<- raster::rowColFromCell(src, cell)
	# expand to cell row /col number indexes
	rownr <- ((rc[1]-1)*fact+1):(rc[1]*fact)
	colnr <- ((rc[2]-1)*fact+1):(rc[2]*fact)
	return(cellFromRowColCombine(targ, rownr, colnr))
}

# quicker version
# return indexes of cells in target raster contained by specified cells in source raster
MapCells <- function(src, targ,
					 cells=1:ncell(src)) # cell indeex
{
	fact <- round(sqrt(ncell(targ)/ncell(src))) # this assummes src is larger than targ, i.e. has been aggregated from it

	#sapply(cells, MARGIN=1, FUN=function(cell){MapCell(src,targ,cell)})
	res<-NULL
	for(cell in cells)
	{
	#		cat(cell, "\n")
			ccells <- MapCell(src,targ,cell)
		 	res <- c(res,ccells)
	}

	return(res)
}

add_layers_with_names <- function(rast, add.layers)
{
  nms <- names(rast)
  nms.add <- names(add.layers)
  
  for(ilayer in 1:nlayers(add.layers))
  {
    rast <- addLayer(rast, add.layers[[i]])
  }
  
  names(rast)<- c(nms, nms.add) 
  return(rast)
}

#
# methods for 3d DEM visualisations
require(fields)
drape.plot.raster <- function(z,    # e.g. dem with z coords
                              z2=z, # matrix supplying draping colours - dims and resolution must match z
                              cols=rev(terrain.colors(25)),
                              box=F, shade=0.5,...)
{
  x <- seq(extent(z)@xmin, extent(z)@xmax, length.out=ncol(z))
  y <- seq(extent(z)@ymin, extent(z)@ymax, length.out=nrow(z))
  z <- as.matrix(z)
  # flip and transpose matrix vals to get in format expected by persp plot
  z <- t(z[nrow(z):1,])
  # matrix supplying colours
  z2 <- t(as.matrix(z2)[nrow(z2):1,])

   drape.cols <- fields::drape.color(z2, col=cols, midpoint=F)$color.index
 # persp(x=x, y=y, z=z, col=drape.cols, box=box, shade=shade,...)

  fields::drape.plot(x=x, y=y, z=z, z2=z2, col=cols, #drape.cols,
             box=box, shade=shade,...)

}

# area of raster that is not NA - e.g catchmet area. Assummes projected coordinate system
area.notNA <- function(rast)
{
  if(is.null(rast@crs) | !is.projected(rast))
  {
    warning("Projected coordinate system required when calculating area. Raster uses geographical system")
  }
  return(xres(rast)*yres(rast)*length(which(!is.na(rast[]))))

}

# return extent object of bounds of raster cell with given index
CellExtent <- function(dem, cell)
{
	ptmid <-xyFromCell(dem,cell, spatial=T)
	# extend to edge of cell
	xymin<- ptmid@coords - c(xres(dem),yres(dem))/2
	xymax<- ptmid@coords + c(xres(dem),yres(dem))/2
	ext <- extent(matrix(cbind(xymin,xymax),ncol=2))
	return(ext)
}

cells.to.poly <- function(rast, cells)
{
	res <- NULL
	for(cell in cells)
	{
		poly <- poly.extent(CellExtent(rast, cell), crs=rast@dem)

		res <- gUnion(res, poly)
	}
	return(res)
}

# draw a box around the specified cell. If fillCol specified colour cell background
HighlightCell <- function(dem, cell, col=NA, label=F, border=par("fg"),...)
{
	ext <- CellExtent(dem,cell)
#	plot(ext,add=add,...)
#	if(!is.na(fillcol))
#	{
		rect(xleft=ext@xmin, xright=ext@xmax, ybottom=ext@ymin, ytop=ext@ymax,
         col=col, border=border, ...)
#	}
  # label with cell number, centred
  if(label)
  {
    xbar <- (ext@xmin+ext@xmax)/2
    ybar <- (ext@ymin+ext@ymax)/2
    text(labels=cell, x=xbar, y=ybar,...)
  }

#	rect(xleft=ext@xmin, xright=ext@xmax, ybottom=ext@ymin, ytop=ext@ymax,...)
}

HighlightCells <- function(dem,cells, col=NA, label=F,
						   border=par("fg"), ...)
{
	if(length(cells)>0)
	{
		# produce a colur vector of same length as the vector of cells to highlight
		fillcols <- rep(col, length.out=length(cells))
		bcols <- rep(border, length.out=length(cells))
		for(cellno in 1:length(cells))
		{
			cell<-cells[cellno]
			col <- fillcols[cellno]
			bcol <- bcols[cellno]
			HighlightCell(dem,cell,
						  col=col, label=label,
						  border=bcol,...)

		}
	}

}

PlotDemADrn <- function(dem, drn, a=NULL, sel=extent(dem),
						grid=F, grid.freq=3,
						grid.col="gray",
						compass = F,
						contour=T,
						nlevels=10,
						compass.pos = "topleft",
						scalebar=F, scalebar.dist=500, ...)
{
	dem <- crop(dem,sel)
	if(!is.null(a))
 	{
		a <- crop(a,sel)
		plot(a,...)
	}
	else
	{
		plot(dem, ...)
	}
	if(grid)
	{
	#	dem <- aggregate(dem,grid.freq)
	#	cells <- which(!is.na(getValues(dem)))

	#	HighlightCells(dem,cells,col=grid.col)
		grid(col=grid.col)
	}
	if(scalebar)
	{
	#	labs<-pretty(0:scalebar.dist,n=5)
		scalebar(d=scalebar.dist,below="m",
			#	 label=labs,
				 divs=4,
				 type="bar")
	}
#	drnexp2<-gIntersection(drnexp, poly.extent(sel, asSpatial=T, dem@crs))
	plot(drn,col="blue",add=T,lwd=2)
	if(contour)
	{
	contour(dem,add=T,nlevels=nlevels)
	}
	if(compass)
	{
		compassRose(x=compass.pos, cex=0.75)
	}
}

# adjust extent object so that its bounds match the grid squares on the raster
SnapExtent <- function(dem, ext)
{

	yres <- yres(dem)
	# shift by remainder of difference of extents
	yshift <- (ext@ymin - extent(dem)@ymin) %% yres
	ymin <-  ext@ymin - yshift
	yshift <- (extent(dem)@ymax-ext@ymin) %% yres
	ymax <- ext@ymax - yshift
	xres <- xres(dem)
	xshift <- (ext@xmin - extent(dem)@xmin) %% xres
	xmin <-  ext@xmin - xshift
	xshift <- (extent(dem)@xmax-ext@xmax) %% xres
	xmax <- ext@xmax + xshift

	return(extent(xmin, xmax, ymin, ymax))
}

# return indexes of cells covered by shape(s)
extract.cells <- function(dem, drn,...)
{
	target <- extract(dem, drn, cellnumbers=T,...)
	if(is.list(target))
	{
		return(do.call(rbind, target)[,1])
	}
	else
	{
		return(target[,1])
	}
}

require(rgeos)

# return a polygon that emcompasses all of the given raster cells
CellBoundingPoly <- function(dem, cells, width=xres(dem)-1)
{
	pts <- SpatialPoints(xyFromCell(dem,cells), dem@crs)
	src <- gBuffer(pts,
				   width=width, byid=T)

	# return one large polygon
	return(gUnionCascaded(src))
}

# create a polygon comprised of the given cells
cells.to.poly <- function(dem,cells)
{
	cells <- unique(cells)
	polys <- lapply(cells,
		function(cell)
		{
			ext <- CellExtent(dem, cell)
			poly.extent(ext, crs=dem@crs, id=cell)
		}
	)
	poly <- do.call(rbind,polys)
	# list of polys
	res<-gUnionCascaded(poly)
	return(res)
}

IntercellDists <- function(dem, from, to)
{
	srcRC <- rowColFromCell(dem, from)
	destRC <- rowColFromCell(dem, to)
	dRC <- srcRC-destRC
	dxy <- dRC * c(xres(dem), yres(dem))
	return(sqrt(rowSums(dxy^2)))
}

require(topmodel)
#
# fill.sinks <- function(dem, deg=0.01, silent=T, fail.if.not.complete=F)
# {
# 	output <- capture.output(res <- raster::setValues(dem,
#                                             topmodel::sinkfill(as.matrix(dem), res=xres(dem), deg=deg)))
# 	if(!silent){message(output)}
# #	if(fail.if.not.complete){stop(output)}
# 	return(res)
#
# }

fill.sinks <- function(dem, deg=0.01,
                       silent=T,
                       ipass=3,   # perform sinkfill a maximum of this times or until all sinks filled
                       fail.if.not.complete=F)
{
    DEM <- as.matrix(dem)
    res <- xres(dem)
    stopifnot(is(DEM, "matrix"))
    if (min(as.vector(DEM[!is.na(DEM)])) < -9000)
        stop("DEM contains unrealistic values (< -9000)")
 #   DEM[is.na(DEM)] <- -9999
    nrow <- dim(DEM)[1]
    ncol <- dim(DEM)[2]

    i <- 1
    sinks.remain <- T
    while(i <= ipass & sinks.remain)
    {
        DEM[is.na(DEM)] <- -9999
        result <- .C("sinkfill", PACKAGE = "topmodel", as.double(DEM),
                 result = double(nrow * ncol + 2), as.integer(nrow), as.integer(ncol),
                 as.double(res), as.double(deg))$result
        result[result > 999998] <- NA
        DEM <- matrix(result[3:(nrow * ncol + 2)], nrow = nrow)
        # 100 is max number of iterations, so if reached thsi then have to run the sinkfill again
        if(result[1]< 100 & result[1]>0)
        {
            sinks.remain <- F

        }
        else if(!silent)
        {
            cat("Sinkfill pass #", i, " No. of sinks remaining = ", result[2], "\n")

        }

        i <- i + 1
    }
    if (result[1] == -1)
    {
        warning("incomplete sink removal")
    }
    if (result[1] == 100)
    {
        msg <- paste("Maximum no. iterations reached (100). No. sinks remaining=", result[2])
        message(msg)
    }

  #  mat <- matrix(result[3:(nrow * ncol + 2)], nrow = nrow)

    return(setValues(dem, DEM))
}





# return a matrix of indexes of<- lls immediately downslope of the given cell in the first column
# and correspoinding flow proportion allocated according to weighted midpoint slope in teh second
# downslope.flow.alls <- function(rast, cur, thresh=0)  #-0.01)  # thresh is maximum slope to allow downslope flow
#   # +ve value allows flow to go "uphill"
# {
#   fd <- terrain(rast, opt="slope")
#   r2 <-
#   w <- matrix(rast, data=raster::xres(rast)*1.4, type="circle")
#   focal(rast, w= w,
#         fun = function(x)
#         {
#           # encode dests
#         }
#         )
#   # encoded by powers of two
#
#   fd.cells <- fd[cur]
#
#   return(adj[,c(1:2,4)])
#
#
#   #	curElev <- rast[cur]
#   # adjacent cells
#   adj <- adjacent(rast,cur,directions=8, pairs=T)
#   #	rastAdj <- rast[adj]
#
#   dz <- rast[adj[,2]] - rast[adj[,1]]
#   bad <- which(is.na(dz) | dz>0)
#   adj <- cbind(adj, dz, NA)
#   if(length(bad)>0)
#   {
#     #	browser()
#     adj <- adj[-bad,]
#     #		dz <- rast[adj[,2]] - rast[adj[,1]]
#   }
#   icell<-1
#   pb <- winProgressBar(max=length(cur), title="flow allocations")
#   on.exit(close(pb))
#   # identify cells numbers downslope from this, if any
#   #	downslope <- which(dz<0)
#   #	ndest <- length(downslope)
#   for(cell in cur)
#   {
#     # downslope destinations  current cell
#     i.adj.cell <- which(adj[,1]==cell)
#     dz.cell <- adj[i.adj.cell,3]
#     #  	if(ndest>0)
# {
# {		# centre point of this cell
#   #locn <- matrix(rep(xyFromCell(rast, cur),ndest),ncol=2, byrow=T)
#   #		locn <- matrix(rep(xyFromCell(rast, cur),length(adj)),ncol=2, byrow=T)
#   # geographical location of downslope cells
#   #		downslopelocn <- xyFromCell(rast, downslope)
#   #		adjlocn  <- xyFromCell(rast, adj)
#   # vector lengths to adjacent upslope cells (midpoints)
#   # pythagoras
#   #dists <- sqrt(rowSums((downslopelocn-locn)*(downslopelocn-locn)))
#   # 		dists <- IntercellDists(rast, adj[,1],adj[,2])     #sqrt(rowSums((adjlocn-locn)^2))
#   # 		# calculate slopes (check not NaN
#   # 		slopes <- dz/dists
#   # 		# contour lengths perp to lines joining midpints - see Quinn et al 1991
#   # 		# these are  easily determined by their relationship with the intercells dists
#   # 		clens <- dists/3
#   # 		# allow flow in downslope directions only, or in directions with slopes
#   # 		# greater than threshold in downslope direction
#   # 		clens[which(slopes>thresh)] <- 0
#   # 		# flow allocations are weighted mean of specific flow
#   # 		sumclens <- sum(abs(clens*slopes))
#   #  sumdz <- sum(dz)
#   #		if(sumclens==0){sumclens<-0}
#   #	p <- abs(clens*slopes/sumclens)
# }

# p <- abs(dz.cell/sum(dz.cell))
# adj[i.adj.cell,4]<-p
# setWinProgressBar(pb, value=icell, label=round(100*icell/length(cur)))
# icell <- icell+1
# # 		return(cbind(adj[downslope,2], p[downslope]))
# }
#   }
# # table: first col is source, second destination, third proportion of flow in that direction
# return(adj[,c(1:2,4)])
# # no downslope cells
# return(NULL)
# }

# clone a raster layer n times, create a stack and apply names
CloneLayer <- function(r, n, names=NULL)
{
	b<-brick(lrep(r, n))
	names(b)<- names
	return(b)
}

height.diff<- function(x)
{
  diff <- x[1] - x[2]
  if(diff > 0){return(diff)}
  else{return(0)}
}

# move to a cell that's downslope from this one
# if random slect teh cell with a probability proportional to the
# contour length * slope
# otherwise, always take the sttepest direction
GetNextDownslopeCell <- function(rast, cur, random)
{
	# numbered 1 to 8,
	adj <- adjacent(rast,cur,directions=8,pairs=F)
	downslope <- adj[which(rast[adj]<rast[cur])]
	ndest <- length(downslope)
	if(ndest>0)
	{
		locn <- matrix(rep(xyFromCell(rast, cur),ndest),ncol=2, byrow=T)
		# geographical location of downslope cells
		downslopelocn <- xyFromCell(rast, downslope)
		# vector lengths to adjacent upslope cells (midpoints)
		#  	browser()
		# pythagoras
		dists <- sqrt(rowSums((downslopelocn-locn)*(downslopelocn-locn)))
		# calculate slopes (check not NaN)
		slopes <- (rast[downslope]-rast[cur])/dists

		# contour lengths perp to lines joining midpoints - see Quinn et al 1991
		# these are  easily determined by their relationship with the intercells dists
		clens <- dists/3

		# (b) randomly select downstream path weighted according to slope
		# (a) choose direction with largest drop (or shallowest uphill...)
		if(random)
		{
			p <- abs(clens*slopes/sum(abs(clens*slopes)))
			flowdest <- WeightedRandomSelection(p)

		}
		else
		{
			# always choose steepest slope (equivalent to single-direction flow algorithm)
			flowdest <- which(slopes==min(slopes))[1]
		}
		nextcell <-  downslope[flowdest]
		return(nextcell)
	}
	#
	return(NULL)
}
################################################################################
# routines for rendering 3d DEM plots and line overlays etc
# ------------------------------------------------------------------------------

# iterate through each line in the lines lying on the DEM and transform to a 2d
# projection of 3d location according to affine transform pmat obtained e.g. from persp functiom

Get3dTransformedDRN <- function(dem, drn, pmat)
{
  coord.list <- GetCoords(drn)
  #coord.list<- unlist(coord.list, recursive=F)

  if(is.list(coord.list))
  {coord.list<-unlist(coord.list, recursive=F)}

  trans<-lapply(coord.list,
                function(crds)
                {

                  # if(is.list(coords)){coords<-coords[[1]]}
                  # browser()
                  x <- crds[,1]
                  y <- crds[,2]
                  z <- extract(dem, crds)
                  crds.dash <- trans3d(x,y,z,pmat)
                })

}


# draw rectangle around a cell in a perpsective plot
HighlightCell3d <- function(dem,    # the dem that
                            cell=NULL,   # cell number or x y coords
                            xy = NULL,
                            pmat,    # affine transformation matrix obtained from call to perps or fields::drape.plot
                            label=NULL,
                            ...)  # col gives fill colour, border the edge colour
{
  if(is.vector(xy))
  {
    # should be a vector of two elements
    cell <- cellFromXY(dem, xy)

  }
  ext <- CellExtent(dem,cell)
  z <- dem[cell]
  xy.dash<-trans3d(c(ext@xmin, ext@xmax), c(ext@ymin, ext@ymax), rep(z, 2), pmat=pmat)
  x.dash <- xy.dash$x
  y.dash <- xy.dash$y

  rect(xleft=x.dash[1], xright=x.dash[2],
       ybottom=y.dash[1], ytop=y.dash[2], ...)
  #	}
  # label with cell number, centred
  if(!is.null(label))
  {
    xbar <- sum(x.dash)/2
    ybar <- sum(y.dash)/2
    text(labels=label, x=xbar, y=ybar,...)
  }
}

Catchment3dOverview <- function(dem,
                                drn=NULL,
                                dem.cols=NULL,    # dem supplying draping colours

                                phi=60, theta=0, border=NA, box=F,
                                add.ctrs=T,
                                nlevels=10, ctrs.lty=2, ctrs.lwd=1, ctrs.col="black",
                                drn.lty=1, drn.lwd=2, drn.col="navy",
                                drn.dash=NULL,  # can supply precalculated transformed drn to save time
                                ...)
{
  if(is.null(dem)){return(NULL)}
  if(!is.null(dem.cols))
  {
    # use supplied raster as basis for colouring faets in perspective plot
    pmat<- drape.plot.raster(dem, dem.cols, box=box, phi=phi, theta=theta, border=border, d=5,...)
  }
  else
  {
    pmat<-persp(dem, box=F, phi=phi, theta=theta, border=border, ...)
  }
  if(add.ctrs)
  {
    # convert dem cvontour lines to spatial lines
    ctrs <- rasterToContour(dem, nlevels=nlevels)
    # debugonce(Get3dTransformedDRN)
    ctrs.dash <- Get3dTransformedDRN(dem, ctrs, pmat)

    lapply(ctrs.dash,
           function(coords)
           {
             lines(coords$x, coords$y, col=ctrs.col,lwd=ctrs.lwd, lty=ctrs.lty)
           })
  }

  if(!is.null(drn))
  {
    cat("Transforming river network...")
    # river network overlain on dem to 2d affine projection using projection matrix determined above

    drn.dash <- Get3dTransformedDRN(dem, drn, pmat)
  }
  if(!is.null(drn.dash))
  {

    lapply(drn.dash,
           function(coords)
           {
             lines(coords$x, coords$y, col=drn.col,lwd=drn.lwd, lty=drn.lty)
           })

  }
  return(list("pmat"=pmat, "drn.dash"=drn.dash)) #ctrs.dash"=ctrs.dash))
}

# set raster cells to the given values, excluding NAs
setValues_notNA <- function(rast, vals)
{

	res <- setValues(rast, vals)
	# this will NA all NA cells in the original
	res <- res + rast - rast
}

# raterizing the vector data with bounds given by the given raster
# use the attribute value given by val
rasterize_layer <- function(layer, rast, val=NULL, zero.val=NA)
{
  # 
  message("Extracting cells occupied by layer...")
  rl <- extract(rast, layer, cellnumbers=T)

  # val gives the attribute value applied to the result raster
  # if not supplied then the rast is 1 where overlapped by shape and zero.val elsewhere
  # cells no occupied by anything are set to the default value supplied
  rast[] <- zero.val
  # firt column the cell number occupied by the layer, second the index of the shape
  # third the atttribute valuel
  rlw <- lapply(1:length(rl),
                function(i)
                {
                  if(is.null(val))
                  {
                    # the value simply indicates whether the shape overlaps the cell
                    val.layer <- 1
                  }
                  else
                  {
                    val.layer <- layer@data[i,val]
                  }
                  cbind(rl[[i]][,1], i, val.layer) 
                  
                })
  # rast.out <- rast
  rlw <- do.call(rbind, rlw)
  
  rast[rlw[,1]] <- rlw[,2]
  
  return(rast)
  
}
