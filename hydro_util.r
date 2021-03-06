require(xts)

#' Time of maximum observation
#'
#' @description Determine the time of the maximum item in the supplied time series.
#' @author Peter Metcalfe
#' @export time_at_peak
#' @param ts Time series
#' @param icol Column index if a multi-column time series
#' @author Peter Metcalfe
#' @examples
#'
#' require(dynatopmodel)
#'
#' data(brompton)
#'
#' with(brompton$storm.run, time_at_peak(qsim))
#'
time_at_peak <- function(ts, icol=1)
{
  tms <- index(ts)
  imax <- which.max(ts)
  return(tms[imax])

  # ttp <- apply(ts, MARGIN=2,
  #   function(x)
  #   {
  #     imax <- which.max(x)
  #     return(tms[imax])
  #   })
  #
  # return(unlist(ttp))
}

#' Time between the peak rainfall and the peak discharge
#'
#' @description Return the lag, in hours, between the peak in the rainfall record and that of the discharge
#' @export time_to_peak
#' @param rain Time series of rainfall.
#' @param qsim Time series of discharges.
#' @param units Units in which to return the value
#' @author Peter Metcalfe
#' @seealso time_at_peak
#' @examples
#'
#' require(dynatopmodel)
#'
#' data(brompton)
#'
#' with(brompton$storm.run, time_to_peak(rain, qsim))
time_to_peak <- function(rain, qsim, units="hour")
{
  qt <- time_at_peak(qsim)
  rt <- time_at_peak(rain)

  ttp <- difftime(qt, rt, units)

  return(as.numeric(ttp))
}


#' Nash Sutcliffe Efficiency of a model's output against observations
#'
#' @description Returns the the NSE (NSE, Nash and Sutcliffe, 1970) of the simulated values against the given observations.
#' @param qsim Time series or vector of simulated values
#' @param qobs Time series or vector of observations
#' @param digits No. decimal places in returned value
#' @return A number <= 1 indicating the goodness of fit of the simulated series against observations (1= perfect fit). Values of >0.8 are generally regarded as "behavioural"
#' @author Peter Metcalfe
#' @import zoo
#' @export NSE
#' @references Nash, J., & Sutcliffe, J. V. (1970). River flow forecasting through conceptual models part I-A discussion of principles. Journal of hydrology, 10(3), 282-290.
#' @examples
#'\dontrun{
#' require(dynatopmodel)
#'
#' data(brompton)
#'
#' # Goodness of fit for the storm simulation
#'
#' NSE(brompton$storm.run$qsim, brompton$storm.run$qobs)
#' }
NSE <- function(qsim, qobs, digits=2)
{
  #qobs <- qobs[index(qsim)]

  qobs <- as.vector(qobs)
	qsim <- as.vector(qsim)

	# shrink to the smallest array
	len <- min(length(qobs), length(qsim))
	qobs <- qobs[1:len]
	qsim <- qsim[1:len]

	igood <- which(!is.na(qobs))
	# test only against non-null observations
	qsim <- qsim[igood]
	qobs <- qobs[!is.na(qobs)]

	if (length(qobs) == 0 || length(qsim) == 0)
	{
	  stop("No non-null observations found")
	}

	res <- 1 - (sum((qobs - qsim)^2)/sum((qobs - mean(qobs))^2))

	res <- round(res, digits)
	return(res)
}

# return a list of efficiency measures for the given run against
run.gof <- function(qsim, qobs, digits=2,...)
{
  qobs <- as.vector(qobs)
  qsim <- as.vector(qsim)

  # fitted linear model
  #  run.vals <- data.frame(qobs=qobs, qsim=qsim)
  #  run.lm <- lm(qsim ~ qobs, data=run.vals)
  # shrink to the smallest array
  len <- min(length(qobs), length(qsim))

  NSE <-  round(NSeff(qobs[1:len], qsim[1:len]), digits)
  R2 <- R2(qobs, qsim)
  logNSE <- NA
  if(all(qobs[1:len]>0 & qsim[1:len]>0))
  {
    logNSE<- round(NSeff(log(qobs[1:len]), log(qsim[1:len])), digits)
  }
  # ratio of the root mean square error to the standard
  # deviation of measured data (RSR),
  res <- list("NSE"=NSE, "logNSE"=logNSE, "R2"=R2)
  return(res)
}

# coefficient of determination
R2 <- function(qobs, qsim)
{
  res<- sum((qsim-mean(qobs))^2)/sum((qobs-mean(qobs))^2)
  return(res)
}


# identify and return as a raster the miniumum downslope flowlengths from all points in teh DEM to teh target cells
#
MinFlowLengthsToTarget <- function (dem,targ)
{
  # test
  if(isLonLat(dem))
  {
    stop("Projected CRS required but using ", dem@crs)
  }

  #plot(dem)

  mat <- as.matrix(dem)

  # minimum lengths to target cells
  lmin <- mat
  lmin[] <- NA

  rowcols<- rowColFromCell(dem,targ)

  n<- length(targ)
  # cells in target are assummed to have nil flow length
  cat("Identifying minimum flow lengths to ", n, " cells\n")
  # identify flowlengths for a cell in target set, then compare with existing
  # min lengths: update wwith new min length
  pb<-txtProgressBar(min=1,max=n,title="",label="", style=3)
  for(rw in 1:n)
  {
    setTxtProgressBar(pb, rw)
    # flow lengths to this cell (NA if none)
    lens<-flowlength(mat, rowcols[rw,])
    lmin <- pmin(as.vector(lens), lmin, na.rm=TRUE)
  }

  dim(lmin)<-dim(mat)

  # create and return min flow length
  # mul by rasterresoultion to get actual length
  ldem <- raster::setValues(dem, lmin)

  return(ldem*xres(dem))
}



# example call
#targ<-extract(dem,wn,cellnumbers=TRUE)[[1]][,1]
# lens<-MinFlowLengthsToTarget(dem,targ)
# cb<-lens>=0
#
# cbound<-rasterToPolygons(cb,dissolve=T)
# plotKML(cbound)
#
# acatch <- round(length(which(!is.na(cb[])))*31*31/1e6,2)
#
#   *31*31/1e6,2))
# aforest<-round(length(targ)*31*31/1e6,2)

# IdentifyAllCatch <- function(dem, nmax=ncell(dem))
# {
#   targ<- 1:ncell(dem)
#   xy <- xyFromCell(dem,targ)
#   gr <- graph.empty() + vertex(targ, "x"=xy[,1], "y"=xy[,2])
#
#
#   mat <- as.matrix(dem)
#
#   rowcols<- rowColFromCell(dem,targ)
#
#   grps<-list
#   for(rw in 1: length(targ))
#   {
#     lens<-flowlength(mat, rowcols[rw,])
#
#     # cells linked downslope to this one
#     linkgr<- which(!is.na(lens))
#
#     if(length(linkgr)>0)
#     {
#       cat("Found connections to cell ", rowcols[rw,], "\n")
#       # add edges and vertices
#   #    browser()
#       conns <- rbind(rep(targ[rw], length(linkgr)),linkgr)
#       gr <- gr + edges(as.vector(conns))
#       # show connected parts of graph
#     #  plot(gr-which(degree(gr)==0), vertex.size=2, vertex.label=NA, edge.width=0.1,
#      #      edge.arrow.size=0.1)
#      # readline()
#     #  grps<-c(grps, c(targ[rw],linkgr))
#     }
#
#   }
#   gr<-gr-which(degree(gr)==1)
#   return(gr)
# }

# saturation vapour pressure at temperature
SatVP <- function(temp)
{
	return(exp((16.78*temp - 116.9)/(temp + 237.3)))
}

# calculates the relative humidity given the wet and dry bulb temperatures.
# source: http://home.fuse.net/clymer/water/wet.html
RelativeHumidity <- function(tdb, Twb, P=101.3)
{
	# 1. Default pressure is taken as 101.3 kPa (kiloPascals).

	# 2. A conversion factor is calculated: A = 0.00066 (1.0 + 0.00115 Twb)
	A <- 0.00066 * (1.0 + 0.00115* Twb)

	# 3. The saturation vapor pressure is calculated at Twb: eswb = e[(16.78 Twb - 116.9)/(Twb + 237.3)]
	eswb <- SatVP(Twb)

	# 4. The water vapor pressure is calculated: ed = eswb - A P (tdb - Twb)
	ed <- eswb - A*P*(tdb - Twb)

	# 5. The saturation vapor pressure is calculated at tdb:
	esdb <- SatVP(tdb)

	#
	# 6. The relative humidity is then calculated: RH = 100 ed/esdb
	return(100*ed/esdb)

	# The algorithm was from a U.S. Water Conservation Laboratory page which no longer exists.
}


# read AWS data from given file and use wet and dry bulb data it contains to
# calculate the rel. humidty. Use row date-times to construct an extended time-series (xts) object
ReadRHseries <- function(fn)
{
	dat <- read.table(fn, header=TRUE, sep="\t")

	RH <- RelativeHumidity(dat$DRY_BULB_TEMP, dat$WET_BULB_TEMP)

	tms <- as.POSIXct(dat$HOUR_ENDED, "GMT", format="%d-%B-%Y %H:%M:%S")
	ser <- xts(RH, order.by=tms)
	return(ser)
}

PlotRHSeries<-function(ser, ...)
{
	plot.xts(ser, ylab="Rel. humidity (%)", major.format="%b %Y", cex.axis=0.75, cex.main=0.75, cex.lab=0.75, ...)

}

