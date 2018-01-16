#' Default list of parameters to control the graphical output during model simulation and via disp.output
#'
#' @description list of parameters to control the graphical output during model
#'   simulation. Parameters with names corresponding to the graphical
#'   parameters returned by par() will be applied to the plot.
#' @export get.disp.par
#' @param ... Further named arguments supplied will overwrite default values
#' @return List of parameters with default values. These include:
#' graphics.show Whether to show graphic output. Default TRUE
#' graphics.window.length Width of display window in hours. Default is 120 days.
#' graphics.interval Interval between refreshing the graphical output, in hours.
#' max.q Max discharge (mm/hr) for display
#' max.rain numeric Max rainfall (mm/hr) for display
#' lwd.rain numeric Line size for rainfall plot
#' int.q Interval between ticks / line on the y axis, in mm/hr
#' int.time Period between ticks on the time axis, a numerical value in hours or one of "day", "week", "month"
#' prop  Proportion of screen occupied by the rainfall hyteograph
#' cex   Overall scaling factor of the plot
#' las.time  Alignment of time axis labels
#' xmar Size of margin on right and left of plot, inches
#' ymar Size of margin above and below plot
#' col  Colours for plot lines: in order, simulated values, observed values
#' @note In this version many of the options are now obsolete. The most relevant are graphics.show and graphics.interval.
#' @examples
#' # Enable graphics output and set display interval to 6 hours
#' par <- get.disp.par(graphics.show=TRUE,
#'        graphics.interval=6)
#'
#' @author Peter Metcalfe
get.disp.par <- function(...)
{
  res <- list("main"= "",
    "sub"= "",
    col.qsim = "blue",
    col.qobs = c("green", "red", "purple", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666"),
    col.evap ="brown",
    lwd.qsim = 2,
    lty.qsim = 1,
    lwd.qobs=1,
    lty.qobs=3,
    int.q=0.2,
    int.evap=0.2,
    int.time=24,
    int.rain=2,
    lwd.rain=3,
    legend.show=TRUE,
    # whether to display run stats at end
    stats.show=TRUE,
    # displayed max discharge in mm/hr
    max.q=4,
    # displayed max rainfall in m/hr
    max.rain=5,
    # minimum displayed ea in m/hr
    min.evap=1,
    # where to send text output - could be a file or console (blank)
    "text.out"="",
    # format for times on axis labels and at the current time
    "time.fmt"="%d-%b",
    "time.text.fmt"="%Y-%m-%d %H:%M",
    "time.pos.fmt"="%H:%M",
    "graphics.window.length"=120*24,      # width of display window (hr)
    "graphics.delay"=10,       # wait before display, expressed in time steps
    "family"="serif",
    "lab.q" = "Specific discharge (mm/hr)",
    "graphics.show" = TRUE,   #  whether or not to show graphic output
#    "graphics.spatial.show"=FALSE,   #  show graphic output
#   "graphics.spatial.output"="qbf", # what to show
    "graphics.interval"=12,     # graphics display interval, in time steps (hrs)
    "graphics.out"= ".",        # save location
    "graphics.save"= FALSE,         # whether to save graphics
    "graphics.save.width"=1024, # dims of output window, pixels
    "graphics.save.height"=768,
    # by how much to offset the y-axis
    "legend.yoff"=-0.05,
    # proportion of screen occupied by the rainfall
    "prop" = 0.25,
    "cex"=1,                        # overall scaling factor of the plot
    # alignment of time axis labels
    "las.time"=3,
    # where to place the labels
    "time.axis.side" = "top",
    # szie of margin on right and left of plot
    "xmar" = 4,
    # size of margin above and below plot
    "ymar" = 4,
    units.discharge = "mm/hr",
#     "graphics.spatial.window.id"=0,
#     "graphics.spatial.3d" = TRUE,   # plot catchment deficit data in 3d
#     "graphics.spatial.theta"=45,  # view angles
#     "graphics.spatial.phi"=30,
#     "graphics.spatial.max.percent" =25,  # max val in % to be coloured in spatial view
#     "graphics.spatial.expand"=0.5,  # vertical expansion factor for 3d view
    "graphics.save.interval"= 48,  # number of time steps between saving
    "cex.axis"=1,          # expansion factor applied to plot axes titles and labels
    "graphics.fn.fmt"="%Y%m%d-%H.jpg",  # format for filenames for saved graphics
    "spatial.interval"=100)

  # merge in any values syupplied as parameters
  res <- merge.lists(res, list(...))

  return(res)
}

# colours for observed values
get.qobs.cols <- function()
{
  cols <- c("green", "red", "purple", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666")
  return(cols)  #[1:nobs])
}

#' get.run.par
#' Initalise model run parameters. Note this function is maintained for backward compatibility only
#' @param dt Numeric Time step in seconds
#' @param tms xts time series or an object that has a POSIXct index
#' @param units string  Units for dt.
#' @param ... Any other parameters returned by get.run.par
#'
#' @return Structure to maintain run information.
#' @details The returned value includes a simulation times calculated from the supplied time range and interval
get.run.par <- function(tms=NULL,
												dt=NULL,
												units="secs",...)
{

  res <- list(
#   "start"=,  # start time / dat
#  "end"=index(rain[length(rain)]),
  "sim.delay"=10,       # hours to run in order to "bed in" model
  "debug"=TRUE,         # debug mode (applies breakpoints etc)
  "unsat.evap" = FALSE,     # whether evapotranspiration takes place from unsat zone at maximum allowable rate  (see Beven 2012)
  "id" = "DTM",
  "start" = "2000-01-01",
  "end" = "2001-12-13",
  #"dt"=900,
  "ntt"=4,                # number of inner time steps
  "log.out"="drm.log")		# location for logging output

  # mwrge in other values
  res <- merge.lists(res, list(...))

  # if a time series supplied (or via an xts object), use this to instantiate the start and end
  if(!is.null(tms))
  {
  	if(is.zoo(tms))
  	{
  		# times supplied via a time series. in this case the default is to use the series' interval
  		tms <- index(tms)
  		if(is.null(dt))
  		{
  			res$dt <- as.numeric(difftime(tms[2], tms[1], units="secs"))
  		}
  	}
  	else
  	{
  		if(is.null(dt))
	  	{
  			# 15 minute interval
	  		res$dt <- 15*60
	  	}
# 	    if(is.vector(tms))
# 		  {
# 		  	start <- tms[1]
# 		  	end <- tms[2]
# 		  }
  	}
  	dt <- res$dt
  	# use the the range supplied or infered to build time series
  	res$start <- tms[1]
  	res$end <- tms[length(tms)]
	}

 	if(units=="secs")
 	{
	  if(dt < 5)
	  {
	  	message("Note: time step converted to seconds")
	  	dt <- dt*3600
	  }
 		dt.units <- dt
 	}
  else if(units=="hours")
  {
	 	if(dt > 5)
		  {
		  	message("Note: time step converted to hours")
		   dt <- dt/3600
	 	}
  	res$dt.units <- dt
  	dt <- dt*3600
  }

	# always in seconds
  res$dt <- dt
  res$start <- as.POSIXct(res$start)
  res$end <- as.POSIXct(res$end)
  # quick way to subset other datasets
  res$sel <- with(res, paste0(start, "::", end))
  res$tms <- seq(res$start, res$end, by=res$dt)

  # could apply a different time interval

  res$nobs <- length(res$tms)

  return(res)
}

# default set of hydrological parameters - applied to HSUs where not supplied
# use underscores rather than periods as ESRI shape format doesn'r support this charcter
# in field names of attached data
def.hsu.par <- function()
{
  list("gauge.id"=1,
  	"srz_max"=0.1,
       "ln_t0"=7,     # saturation transmissivity
       "m"=0.01,
       "srz0"=0, "td"=1,
       "vchan"=1000,
       "vof"=100,
       "k0"=1e8,     # surface conductivity - not used in DynaTM. Used in infiltration excess calcs in Buytaert's implementation of TOPMODEL. Set to large value for no infiltration excess
       "CD"=0.1,     # capillary drive (see Morel-Seytoux and Khanji, 1974, Beven, 1984), as above.
       "atb.bar"=0,  # areal average of topographic index
       "sd_max"=0.5, # max storage deficit
  	"pe_fact"=1,
  	"vof_fact"=1,
  	"rain_fact"=1,   # scaling factor for rainfall input
  	"mann.n"=0.01,    # roughness factor for rough grass
  	"S0"=0.1,       # nominal surface slope
  	"ex_max"=1)      # maximum surface excess storage
}

# default set of parameters applied to channel HSU
def.chan.par <- function(vals=list())
{
  res <- list("srz_max"=0, # rainfall goes straight to the unsat zone infiltration
       "ln_t0"=8,
       "m"=0.01,
       "srz0"=0,
       "td"=1e-8 ,
       "vchan"=1500,
       "vof"=10,
       "sd_max"=2)    # depth of rectangular channel?
  merge.lists(res, vals)
}
