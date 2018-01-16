##########################################################################################
# Routines for reading time series data (e.g. obs, p.e and discharges) from disk
##########################################################################################
check.obs <- function(input, msg=FALSE)
{

#   # return number of missing values
  return(length(which(is.na(input))))
}

# fill in missing values by interpolation and update flags
fix.obs<- function(obs, maxgaps=24, missing.data=NA)
{
  bad <- which(is.na(obs))
  if(length(bad)>0)
  {
  	# if a non-na value  missing.data will replace the gaps
  	obs[bad]<- missing.data
	  nbad <- length(bad)
	  cat(nbad, " missing obsfall data... attempting interpolation...")

	  # linear interpolation
	  obs <- na.approx(obs, maxgap=maxgaps)

	  # mark interpolated data
#	  obs[bad]$status <- -1  #"I"

	  nbad2 <- length(which(is.na(obs)))
	  nfixed <- nbad - nbad2

	  cat(nbad2, " missing data remaining\n")

	  message(nfixed, " missing data interpolated, ", nbad2, " remaining")
  }
  return(obs)
}

# loading tabulated data from the specified file
load_data_table <- function(fn, sep="\t", header=FALSE, 
         skip=1, nfield=NULL, 
         fields=NULL, fixed=TRUE,
         fact=1)
{
  # read all lines
  dat <- readLines(fn)
  istart <- skip
  
  # check for repeated separators that indicate empty data
  # dat <- gsub(paste0(sep, sep), paste0(sep, " ", sep), dat)
  if(is.null(nfield))
  {
    if(header)
    {
      # this assume that the first line contain the field information
      fields <- strsplit(dat[skip], sep)[[1]]
      nfield <- length(fields)    
      istart <- skip+1
      nms <- fields
    }
  }
  else
  {
    nfield <- length(dat[skip])
    
    # first line with data sets the expected size per 
    istart <- skip+1
    nms <- paste0("X", 1:nfield)
    # could keep reading lines until finding a line with the requisite number of fields
    #    nfield <- length(fields)    
  }
  
  if(length(dat)<istart)
  {
    # data.frame(matrix(ncol=length(fields), nrow=0, byrow=TRUE))
    return(NULL)
  }

  pb <- txtProgressBar(min=istart-1, max = length(dat), style = 3)
  
  dat.good <- lapply(istart:length(dat),
                     
                     function(i)
                     {
                       setTxtProgressBar(pb, value = i)
                       
                       fields <- strsplit(dat[i], sep)[[1]] 
                       
                       if(length(fields)==nfield)
                       {
                         return(matrix(fields, nrow=1)) 
                       }
                       # if a different no. items ignore
                       return(NULL)
                     }
  )

  # calculating the number of rows that failed
  nbad <- (length(dat)-istart)-length(dat.good)
  
  if(nbad>0)
  {
    message("Couldn't read ", nbad, " records; removed")
  }
  message("Loaded ", length(dat.good), " record(s), merging....")
  dat.good <- data.frame(do.call(rbind, dat.good))
  
  # always now set
  names(dat.good)<- nms
  
  # identify numeric columns
  for(i in 1:ncol(dat.good))
  {
    if(all(is.number(dat.good[,i])))
    {
      dat.good[,i] <- as.numeric(dat.good[,i])
    }
  }
  
  return(dat.good)
}

# load the calibration results file and ensure all lines have the same length (set by the first line)
# remove any that fail. skip gives the numeber of initila lines to ignore
load_results_table <- function(fn, sep="\t", header=FALSE, 
                            itm=1, ival=2, icol=1:2,  # time column; value column(s); cols to return (including the time)
                            as.xts=TRUE,
                            skip=1, nfield=NULL, 
                            fields=NULL, 
                            fixed=FALSE,
                            fact=1,
                            start=NULL, end=NULL,
                            fmt="%Y-%m-%d %H:%M:%S",...)
{
  # read the first line
  dat <- readLines(fn)
  istart <- skip
  
  # check for repeated separators that indicate empty data
  # dat <- gsub(paste0(sep, sep), paste0(sep, " ", sep), dat)
  
  if(is.null(nfield))
  {
    if(header)
    {
      # this assume that the first line contain the field information
      fields <- strsplit(dat[skip], sep)[[1]]
      nfield <- length(fields)    
      istart <- skip+1
    }
    
  }
  else
  {
    # first line with data
    istart <- skip+1
    # could keep reading lines until finding a line with the requisite number of fields
    #    nfield <- length(fields)    
  }
  
  # always returning the column containing tm 
  # cols <- unique(c(col, itm))
  
  pb <- txtProgressBar(min=istart, max = length(dat), style = 3)

  iend <- length(dat)
  if(!is.null(end))
  {
    end <- as.POSIXct(end)
    # search for first time that matches 
    tm.str <- format(end, format=fmt)
    iend.test <- grep(tm.str, dat)
    if(length(iend.test)>0)
    {
      iend <- iend.test
    }
  }
  
  if(!is.null(start))
  {
    start <- as.POSIXct(start)
    # search for first time that matches 
    tm.str <- format(start, format=fmt)
    istart.test <- grep(tm.str, dat)
    if(length(istart.test)>0)
    {
      istart <- istart.test
    }
  }

  
  dat.good <- lapply(istart:iend,
                     
                     function(i)
                     {
                       setTxtProgressBar(pb, value = i)
                       
                       fields <- strsplit(dat[i], sep)[[1]] 
                       
                       add <- length(fields)==nfield | fixed
                       
                        if(add)
                        {
                         # return subset
                         return(fields[icol])
                       }
                       
                     }
  )
  
  nbad <- (length(dat)-istart)-length(dat.good)
  
  if(nbad>0)
  {
    message("Couldn't read ", nbad, " records; removed")
  }
  message("Loaded ", length(dat.good), " records, merging....")
  dat.good <- data.frame(do.call(rbind, dat.good))
  
  if(header & is.null(fields))
  {
    nms <- fields
    names(dat.good)<- nms
  }
  
  if(as.xts)
  {
    # first col by default holds the time uindex
    tms <- strptime(dat.good[,itm], format=fmt)
    # assuming firther cols are numeric 
    message("Observations in period ", dat.good[1,itm], " to ", dat.good[iend,itm])
    vals <- as.numeric(as.character(dat.good[,ival]))*fact
    res <- xts(vals, order.by=tms)
    return(res)
  }
  
  return(dat.good)
}

# create time series from given input , start time and interval
# if already a time series, truncate series by teh start and check at same resolution
# otherwise build the
get.obs <- function(obs, sim.start, sim.end=NULL, dt=1)  # dt is in hours
{
	if(length(obs)==0)
	{
		# must have both limits if nothing supplied here
		tms <- seq(sim.start, sim.end, by=dt*3600)
		message("Time series input created ") #, obj.name(obs))
		obs <- xts(rep(0, length(tms)), order.by=tms)   #  could attempt to
	}
	else if(is.zoo(obs))
	{
		#obs <- obs
	}
	# if arrays without time info supplied attamept to create time series
	else if(is.matrix(obs) & !is.null(sim.start))
	{
		tms <- seq(sim.start, length.out=nrow(obs), by=dt)
		obs <- xts(obs, order.by=tms)

	}
	else if(is.vector(obs) & !is.null(sim.start))
	{
		tms <- seq(sim.start, length.out=length(obs), by=dt)
		obs <- xts(obs, order.by=tms)
	}
	# fix any NAs - warn if this
	which.na <- which(is.na(obs[]))
	if(length(which.na)>0)
	{
		warning(paste(length(which.na), " NA observations replaced by zeroes"))
		obs[which.na] <-0	 # deals with multi dimensional obeject
	}

	return(obs)
}

# read time string and parse to date time.
# ParseDateTime <- function(time_str, tz="GMT", fmt="%Y-%m-%d %H:%M:%S",
# 						  revorder=TRUE)
# {
# 	time <-strptime(time_str, fmt, tz=tz)
#
# 	if(any(is.na(time)))
# 	{
# 		nerrs <- length(which(is.na(time)))
# 		warning(paste(nerrs, " strings not parsed succesfully into date-times"))
# 	}
# 	   return(time)
#
# }

# write a subset of obsfall (or other time series) record back to the specified location
write.obs <- function(fn, dat,
                          start = index(dat)[1],
                          end = index(dat)[nrow(dat)], sep="\t",
                          which.col=1:ncol(dat), ...)
{
#  if(is.null(start)){start<-}
  tz <- slot(index(dat), "tzone")
  end <- as.POSIXct(end, tz=tz)
  start <- as.POSIXct(start, tz=tz)
  tms<- index(dat)
  sel <- which(tms<=end & tms>=start)
  subdat <- dat[sel, which.col]
 # tms <- index(dat)[sel]
  write.table(data.frame(subdat), fn, quote=FALSE, sep=sep,...)

}

# test val to see if it can be co-erced
is.number <- function(s)
{
  w <- as.numeric(options("warn"))
  options("warn"=-1)
  val <- FALSE
  try(val <- !is.na(as.numeric(s)), silent=TRUE)
  options("warn"=w)
  return(val)

}

# test the first line of the file to see if its is a hader of meta data
# return a value > 0 indicating the number of header line found
header.lines <- function(fn, sep, max.n=2, skip=0)
{
  res <- -1  # default value inicates no input
  lines <- as.list(readLines(fn, max.n+skip))
  lines <- lines[(skip+1):(max.n+skip)]
  lines.vals <- sapply(lines,
    function(l)
    {
      # keep splitting the line
      s <- strsplit(l, split=sep)
      s2 <- sapply(s, strsplit, split=sep)
      return(any(sapply(s2, is.number)))

    }
  )
  seqs <- rle(lines.vals)
  # lengthg of the first sequence of FALSE in file, indicating a header line comprised of
  if(!seqs$values[1]){
    return(seqs$lengths[1])
  }
  return(res)
}


# read a single tabulated observation datum at date-times given in the specified column
read.obs <- function(fn,
                          #  start=NULL, # return subset of
                          #  end=NULL,
                            maxgaps=24, # largest gap in data that can be interpolated
                            itm=1,   # col name, col number or vector (start, end, number is inferred from number of observations. If missing, generate a time series from start and end)
                            val=NULL,  # col number(s) or name(s) containing observations. if > number of colums then truncate
                            fmt="%Y-%m-%d %H:%M:%S",   # date time format string
                            tz="GMT",           #  tz defaults to UTC, blank uses current tz - watch out for daylighgt saving time changes
                            header=FALSE, #is.character(tm) | is.character(val),           # whether column names are present, can also point to a line in the input
                            metalines=0,        # any metadata lines - can also point to a file location. synonym with skip
                            skip=NA,
                            units="m",          # expected units
                            dt=1,               # default time step in hours
                            sep="\t",
										        missing.data = NA,  # missing data values: set to e.g. 0 to treat as no reading
                            include.status=FALSE,...)   # include an integer status code for each observation in another colums
{
  skip <- max(metalines, skip, na.rm=TRUE)
#  try(header <- header.lines(fn, sep, skip=metalines)>0)
  # determine and override type of separator
  sep.rep <- switch(tools::file_ext(fn),
         "csv"=",",
         "tsv"="\t",
         "dat"="\t")
  if(!is.null(sep.rep)){sep<-sep.rep}

  colnames <- NULL
  if(is.numeric(header))
  {
      # numeric value  indicate ro on which to find the colum header names
      colnames <- readLines(fn, n=header)[header]
    header <- FALSE
    colnames <- strsplit(colnames, sep) [[1]]
  }

  # is the first line a numeric
#  first <- split(first.line, split=sep)
#  try(if(is.character(first.line[[2]])){header=TRUE})

	obs <- read.table(fn,
                     fill=FALSE,
	                 header=header, # if using named columns then this must be TRUE
	                 sep=sep,   # or delimited if sep=","
	                 skip=skip,...)    # how mant lines to skip at head

  if(length(itm)>0)
  {

    # time column specified
    # e.g. vector of start, end times
#     times <- seq(start, end, length.out = nrow(obs))
	  # check for any column headers, if not found then assumme that first col is time and second the reading
	  time_str<-  obs[,itm]  #paste(obs$year,"-",obs$month,"-",obs$day," ",obs$hr,":00:00",sep="")
    if(length(itm)>1)
    {
      # collapse into a single column assuming e.g first is date and second the time
      time_str<-do.call(paste, time_str)
    }
#     start <- obs[1, tm[]
#     end <- obs[nrow(obs), tm]
	  # parse to local times (POSIXlt), then calendar times ()
	  times <- strptime(time_str, fmt, tz=tz)   #, tz="GMT")
  }
  else
  {
    # assume given time step (1 hr by default) and create a dummy time series
    start <- as.POSIXct("2000-01-01")
    end <- as.POSIXct(start+dt*3600*(nrow(obs)-1))
    times <- seq(start, end, by=dt*3600)
    message(paste("Time index created for series from ", start, " to ", end, " by ", dt, "hours"))
  }
  # cam specify all colukms byu supplying a large range of column indexes containing observations
  # will be truncated to the column count
  if(is.numeric(val)& max(val)>ncol(obs))
 	{
  	val <- min(val, ncol(obs)):ncol(obs)
  }
	else
	{
	  # if val is null or invalid then use all columns
	  val <- setdiff(1:ncol(obs), itm)
	}
    if(!is.null(colnames))
    {
        nms <- colnames[val]
    }
    else
    {
        nms <- colnames(obs)[val]
    }
  # ensure numeric, char arguments will become NAs
  obs <- apply(as.matrix(obs[,val]), MARGIN=2, function(x){as.numeric(x)})
  # strip out invalid cols
  OK.cols <- apply(obs, MARGIN=2, function(x){!all(is.na(x))})
  obs <- obs[, OK.cols]
	# second column is the hourly measured obsfall. Create time series
	obs <- xts(obs, order.by=as.POSIXct(times))
	# determine ends of period
#	if(is.null(start)){start<-times[1]}
#	if(is.null(end)){end<-times[length(times)]}

#	sel<-paste(start, "::", end, sep="")
  # remove data whose time entries couldn't be parsed
	badtimes <- which(is.na(times))
	if(length(badtimes)>0)
  {
	  obs <- obs[-badtimes,]
  }
	if(length(obs)==0)
	{
    stop("No valid times located in data source - check time string format ", fmt, "\n")
	}
 # obs <- apply(obs, MARGIN=2, FUN=fix.obs, maxgaps=maxgaps, missing.data=missing.data)
  n <- nrow(obs)
  cat(paste("Read ", length(obs), " records for period ",
      index(obs[1]), " to ",  index(obs[n]), " no. missing=", check.obs(obs), "\n"))

  if(units=="mm" & max(obs>1)){
    message("Data may be in mm")
  }
  colnames(obs)<-nms[OK.cols]

  # just the observations
  return(obs)
}


# read in observed discahrge data in IH format then aggregate to the time step used in simulation
# defaults to 1hour. Assmumes standard CEH format with 14 lines of metadata
ReadDischargeData <- function(fn, dt=1, start=NULL, end=NULL,
							  aggregate =FALSE,
							  headerlines=14)
{
	flow <- read.csv(fn, header=FALSE, skip= headerlines)
	# remove index column if it exists
	#   if(ncol(flow)>2)
	#   {
	#     flow <- flow[,-1]
	#   }
	# "unfactor" the first column
	flow[,1] <- as.character(flow[,1])

	# parse the first col into dates - should be in UTC format e.g. 1969-11-21T22:00:00
	# see http://www.cl.cam.ac.uk/~mgk25/iso-time.html
	#datetimes<- strptime(flow[,1],"%Y-%m-%d",tz="GMT")
	datetimes<- strptime(flow[,1],"%Y-%m-%dT%H:%M:%S")  #,tz="GMT")
	#   inperiod <- which(datetimes >= start & datetimes <= end)
	#   # select subset
	#   datetimes <- datetimes[inperiod]
	vals <- flow[,2]
	names(flow) <- c("date", "flow")

	cat(nrow(flow), " records loaded... creating time series...\n")

	#vals <- which(date)
	# create a time series (daily mean flow in m^3/hr)
	sel<-paste(start, "::", end, sep="")
	Qobs <- xts(vals,order.by=datetimes)

	if(!is.null(sel)){Qobs <- Qobs[sel]}

	if(aggregate)
	{
		# convert to aggregated hourly data (means)
		l <- 3600*dt
		cat("aggregating to every", dt, " hour(s)...\n")
		Qobs <- aggregate(Qobs, time(Qobs) - as.numeric(time(Qobs)) %% l, mean)}
	return(xts(Qobs))
}


# load the calibration results file and ensure all lines have the same length (set by the first line)
# remove any that fail. skip gives the numeber of initila lines to ignore
# version that can deal with missing lines
# fixed TRUE if should attempt to truncate / expand any data to the given number of fields
load_results_table <- function(fn, sep="\t", header=FALSE, 
                               itm=1, ival=2, icol=1:2,  # time column; value column(s); cols to return (including the time)
                               as.xts=TRUE,
                               skip=1, nfield=NULL, 
                               fields=NULL, 
                               fact=1,
                               fixed=FALSE,
                               fmt="%Y-%m-%d %H:%M:%S",...)
{
  # read the first line
  dat <- readLines(fn)
  istart <- skip
  
  # check for repeated separators that indicate empty data
 # dat <- gsub(paste0(sep, sep), paste0(sep, " ", sep), dat)
  
  if(is.null(nfield))
  {
    if(header)
    {
      # this assume that the first line contain the field information
      fields <- strsplit(dat[skip], sep)[[1]]
      nfield <- length(fields)    
      istart <- skip+1
      
    }

  }
  else
  {
    # first line with data
    istart <- skip+1
    # could keep reading lines until finding a line with the requisite number of fields
#    nfield <- length(fields)    
  }

  if(length(dat)<=istart)
  {
   # data.frame(matrix(ncol=length(fields), nrow=0, byrow=TRUE))
    return(NULL)
  }
  
  # always returning the column containing tm 
 # cols <- unique(c(col, itm))

  pb <- txtProgressBar(min=istart, max = length(dat), style = 3)
  
  dat.good <- lapply(istart:length(dat),
     
         function(i)
         {
           setTxtProgressBar(pb, value = i)
           
           fields <- strsplit(dat[i], sep)[[1]] 
           
           
           if(length(fields)==nfield | fixed)
           {
             fields <- sapply(fields[icol],
                    function(f)
                    {
                      
                      if(is.number(f))
                      {
                        return((as.numeric(f)))
                      }
                      return(f)
                    }
             )
             return(fields)
           }
           
         }
  )

  nbad <- (length(dat)-istart)-length(dat.good)
  
  if(nbad>0)
  {
    message("Couldn't read ", nbad, " records; removed")
  }
  message("Loaded ", length(dat.good), " records, merging....")
  dat.good <- data.frame(do.call(rbind, dat.good))
  
  if(header & is.null(fields))
  {
    nms <- fields
    names(dat.good)<- nms
  }
  
  if(as.xts)
  {
    # first col by default holds the time uindex
    tms <- strptime(dat.good[,itm], format=fmt)
    # assuming firther cols are numeric (ensure vals are converted from levels)
    vals <- as.numeric(as.character(dat.good[,ival]))*fact
    res <- xts(vals, order.by=tms)
    return(res)
  }
  
  return(dat.good)
}
