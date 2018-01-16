#load.source("build_hru.r")
#load.source("read_obs.r")
# p.e. calculations
#load.source("evap.r")
#load.source("defs.r")
#load.source("results.r")
#load.source("calib_DTM.r")

#load.source("dtm.main.r", chdir=T)

require(raster)
require(xts)

get.sim.range <- function(proj)
{
  #  require(intervals)
  obs <- proj$obs$obs
  pe <- proj$obs$pe
  qobs <- proj$obs$qobs
  try(proj$sim.start <- as.POSIXct(proj$sim.start), silent=T)
  try(proj$sim.end <- as.POSIXct(proj$sim.end), silent=T)

  if(length(proj$sim.start)==0)
  {
    s1 <- NA
    s2 <- NA
    s3 <- NA
    try(s1 <- start(obs), silent=T)
    try(s2 <- start(pe), silent=T)
    try(s3 <- start(qobs), silent=T)

    proj$sim.start <- max(c(s1, s2, s3), na.rm=T)
    if(!is.null(proj$sim.start))
    {
      message(paste("Start of simulation inferred from input as ", proj$sim.start))
    }
  }
  if(length(proj$sim.end)==0)
  {
    e1 <- NA
    e2 <- NA
    e3 <- NA
    try(e1 <- end(obs), silent=T)
    try(e2 <- end(pe), silent=T)
    try(e3 <- end(qobs), silent=T)
    proj$sim.end <- min(c(e1, e2, e3), na.rm=T)
    if(!is.null(proj$sim.end))
    {
      message(paste("End of simulation inferred from input as ", proj$sim.end))
    }
    if(proj$sim.end < proj$sim.start)
    {
      stop("Error: sim.start after end. Check supplied data and values")
    }
  }
  return(proj)
}



# use catchment area calculated from dem to determine specific discharges from the
#input given in cu.m/sec
convert.to.specific.discharges <- function(proj, q)
{
	if(max(q))
	# catchment area
  a <- with(proj, sum(length(which(!is.na(dem[])))*xres(dem)*yres(dem)))
  # assumme in cu.m/s
  res <- 3600*q/a
  if(max(res, na.rm=T)>1){warning("Very large specific discharges calculated: check input not in mm")}

  return(res)
}

load.obs <- function(proj,  obs.dir=proj$obs$dir)
{
  if(!is.null(obs.dir) & file.exists(obs.dir))
  {
    proj$obs$dir <- obs.dir
    rain <- NULL
    pe <- NULL
    qobs <- NULL
    qobs.fn <- dir(obs.dir, "qobs*.*$")
    rain.fn <- dir(obs.dir, "rain*.*$")
    pe.fn <- dir(obs.dir, "pe*.*$")
    # specfying large number  - will be lmited to number present. load data for all locations
    try(rain <- read.obs(fp(obs.dir, rain.fn[1]), val=2:10, sep="\t"), silent=F)
    try(qobs <- read.obs(fp(obs.dir, qobs.fn[1]), val=2:10, sep="\t"), silent=T)
    try(pe <- read.obs(fp(obs.dir, pe.fn[1]), val=2:10), silent=T)
    proj$obs$rain <- rain
    proj$obs$qobs <- qobs
    proj$obs$pe <- pe
  #  try(proj$qmax <- max(qobs, na.rm=T), silent=T)
#    try(proj <- fix.run.dates(proj))
  }
  return(proj)
}

# second go at load project using more logical structure
create.proj <-
    function(data.dir,  # location for DEM, drn etc
             id=data.dir,
             cuts=NULL,      # cuts applied for first discretisation
             area.thresh=1,
             chan.width=1,
             disc.dir=file.path(data.dir, "disc"),   #  root location for discretisations

    				 # specify one or more discretisations by name. other directories wil be ignored
    				 discs=NULL,
    				 # obseravtions e.g rain, pe and discharge. can be from multiple sites identified by id
    				 obs = NULL,
    				 # default location
             obs.dir=fp(data.dir, "obs"),
             #rebuild=F,
             ...)
    {
        # ensure relative paths converted to absolute
        data.dir <- path.expand(data.dir)
        #	sink(type = "message")
        #  on.exit(sink(NULL))
        if(!file.exists(data.dir))
        {
            if(readline(paste("Specified location ", data.dir, " not found, create?"))=="y")
            {
                dir.create(data.dir, recursive=T)
            }
            else
            {
                stop("Can't find directory")
            }
        }
        proj <- new.project()
        proj$name <- id
        proj$dir<-data.dir
        # root folder where discretisation live
        proj$disc.dir<-disc.dir
        # add in extra parameters
        proj <- merge.lists(proj, list(...))

        #  try(proj$dem <- rast.mem.copy(raster(file.path(data.dir, "dem.tif"))))
        #try(proj$drn <- rgdal::readOGR(data.dir, "drn"))
        # iterate through any shape files located and add verbatim
        for(shp in dir(data.dir, "*.shp"))
        {
            # name without extention
            shp.nm <- sub("*.shp", "", shp)
            cat("Loading shape file ", shp.nm, "...")
            try(proj[[shp.nm]] <- rgdal::readOGR(data.dir, shp.nm, verbose=F))
            if(!is.null(proj[[shp.nm]])){cat("done\n")}
            else{cat("failed\n")}
        }

        # might be a routing table applicable to every dicsretisation
        proj$routing <- load.disc.data(file.path(data.dir, "routing.dat"), header=F, rebuild=F)
        # catchment layers that can be used to provide data for discretisations
        catch <- NULL
        # add raster verbatim and copy to memeory to prevent hard coding disk referencve
        for(rast.fn in dir(data.dir, "*.tif", full.names=T))
        {
            rast <- NULL
            try(rast <- raster(rast.fn))
            #   try(rast <- rast.mem.copy(rast))
            # name for list element is raster file name with out path and extension
            rast.nm <- sub("*.tif", "", basename(rast.fn))
            proj[[rast.nm]] <- rast
            message(paste("Adding raster ", rast.fn))
            if(is.null(catch))
            {
                catch <- stack(rast)
            }
            else
            {
                try(catch <- addLayer(catch, rast))
            }
            #  try(proj$reaches <- raster(file.path(data.dir, "reaches.tif")))
        }
        proj$catch <- catch
        # build
        #  try(proj$dem <- raster::setValues(proj$dem, sinkfill(as.matrix(proj$dem), res=xres(proj$dem), degree=0.1)))
        # any files ending with ".par" are assummed to contained R source that can be loaed and parsed as is
        par.fn <- dir(data.dir, "*.par$", full.names=T)
        for(fn in par.fn)
        {
            #  envir <- as.environment(proj)
            cat("Trying settings from ", fn, "...\n")
            try(source(fn, local=proj, echo=T), silent=F)
        }
        # rainfall, pe and observed flows
        proj <- load.obs(proj, obs.dir)
        # attempt to set the range
        proj <- get.sim.range(proj)
        # load cuts and other discretisation info from subdirectories of specified directory or use default location "disc"
        # list of discretisations excluding any observations
        if(length(discs)>0)
        {
        	# supplied one or discretisation (named by dir locn)
        	dns <- lapply(as.list(discs), function(x)file.path(disc.dir, x))
        }
        else
        {
        	# try to load a discretisation from every subdirectory located in disc dir
        	dns <- setdiff(list.dirs(disc.dir, recursive=F), c(disc.dir, obs.dir))
        }
        proj$disc <- lapply(dns,
                            function(dn)
                            {
                                disc <- NULL
                                cat("trying directory ", dn, "...")
                                # look in every subfolder of the discretisation root dir and attempt to load a
                                # discretisation. if it fails the error will be returned in disc
                                try(disc <- disc.from.dir(dem=proj$dem, drn=proj$drn,
                                                          reaches=proj$reaches, dn=dn, cuts=NULL,
                                                          #rebuild=rebuild,...
                                                          chan.width=chan.width,
                                                          catch=catch,
                                                          routing=proj$routing,
                                                          area.thresh=area.thresh, dn.out=NULL)
                                )
                                return(disc)
                                #if(!is.null(disc)){cat("success")}
                            }
        )
        # add in any new discretisations specified via teh cuts function arguments
        if(!is.null(cuts))
        {
            proj <- add.disc(proj, cuts, ...)
        }
        # remove NuLL and invalid discs then order by no. HRUs
        try(proj$disc <- delete.NULLs(proj$disc))
        # vector of HRU counts
        try(proj$hru.counts <- get.disc.counts(proj$disc))
        # locate errors when loading discs
        #	try(errs <- lapply(proj$disc, function(disc){if(!is.list(disc)){return(disc)}}))
        # remove failed disc loaded and repiort
        #	cat("\n", length(proj$disc), " discretisation(s) loaded\n")

        # any optional parameters supplied will override the default settings
        proj <- merge.lists(as.list(proj), list(...))
        proj$disp.par$title.main <- proj$name
        return(proj)  # or convert to lits
    }

# create a new discretisation and add to project's list. if it exists already then just add
# to list of discretisations, unless rebuild is true in which case rebuid the components
add.disc <- function(proj, cuts=NULL, i.disc=0, rebuild=F, chan.width=2,
										 ...)
{
 # proj <- merge.lists(proj, list(...))

	# build names from extra parameters
	cuts <- merge.lists(cuts, list(...))

  # name wil be infered from cuts
  dn <- disc.dir.name(cuts, dn=proj$disc.dir)

 # rebuild <- file.exists(dn)
  if(rebuild)
  {
  	message(paste("Rebuilding files in ", dn))
  	proj$reaches <- build.reach.raster(proj$dem, proj$drn, chan.width=chan.width)
  }
  disc <- NULL
  disc <-disc.from.dir(dem=proj$dem,
                           drn=proj$drn,
                           reaches=proj$reaches,
  												 routing=proj$routing,
  												 catch=proj$catch,
                           dn=dn,
                           cuts=cuts,
                           rebuild=rebuild,
                           chan.width=chan.width,
                           area.thresh=proj$area.thresh)  #, silent=F)

  if(is.null(disc))
  {
#     if(i.disc<=0 | i.disc > length(proj$disc))
#     {
#       # append to list, otherwsie specified disc is replaced
#       i.disc <- length(proj$disc) + 1
#     }
#     # same cuts as any other?
#     proj$disc[[i.disc]]<- disc

    warning("Adding discretisation ", paste(cuts, collapse=","), "failed")
  }
 return(disc)
}

rebuild.disc <- function(proj, i.disc=1, disc=NULL, dn=disc$dir,
                         what=c("weights.dat", "groups.dat", "disc.par", "reaches.tif"),...)   # specify what is to be rebuilt
{
  if(is.null(disc)){disc <-proj$disc[[i.disc]]}
  # add in any new or updated data
  args <- list(...)
  disc <- merge.lists(disc, args)
  write.disc(disc)
  dem <- proj$dem
  drn <- proj$drn
  # remove specified files
  lapply(file.path(dn, what),
         unlink)

  disc <- disc.from.dir(disc$dir, dem=dem, drn=drn,
                        area.thresh=disc$area.thresh,
                        chan.width=disc$chan.width,
                        cuts=disc$cuts)

  return(disc)
}

add.layers <- function(proj, reload=F)
{
	nms <-  dir(proj$dir, "\\.shp")
	if(!reload)
	{
		nms <- setdiff(nms, paste0(names(proj), ".shp"))
	}
	#  try(proj$dem <- rast.mem.copy(raster(file.path(data.dir, "dem.tif"))))
	#try(proj$drn <- readOGR(data.dir, "drn"))
	# iterate through any shape files located and add verbatim
	for(shp in nms)
	{
		# name without extention
		shp.nm <- sub("*.shp", "", shp)
		cat("Loading shape file ", shp.nm, "...")
		try(proj[[shp.nm]] <- readOGR(proj$dir, shp.nm))
		if(!is.null(proj[[shp.nm]])){cat("...done\n")}
		else{cat("...failed\n")}
	}
	return(proj)
}

# save contents of a project to disk in tif, shp so it can be read by other software
write.proj <- function(proj,
                       dir = proj$dir,  # location for DEM, drn etc
                       disc.dir = file.path(dir, "disc"),
                       obs.dir = file.path(dir, "obs"),
											 save.disc=T,
                     #  copy.source = T,  # if true the original files creating the discretisation will be copied to dir
                       save.obs=F)    # obseravtions
{
  if(!file.exists(dir)){dir.create(dir, recursive=T)}
  proj$disc.dir <- disc.dir
  proj$obs.dir <- obs.dir
  proj$data.dir <- dir

  if(!file.exists(proj$disc.dir)){dir.create(proj$disc.dir, recursive=T)}

  proj$sim.start <- as.character(proj$sim.start)
  proj$sim.end <- as.character(proj$sim.end)
  # save options can be later reloaded with source
  try(dump(list=c("disp.par", "run.par", "dt", "ntt", "qmax", "qt0", "sim.start", "sim.end"), envir=as.environment(proj),
       file=file.path(proj$dir, "proj.par")))

  #  write.table(proj$disc$w, file.path(proj$disc$dir, "weights.dat"), quote=F, col.names=T, row.names=F, sep="\t")

  #  write.table(proj$disc$groups, file.path(proj$disc$dir, "groups.dat"), quote=F, row.names=F, sep="\t")
  if(save.obs)
  {
    write.zoo(proj$obs$qobs, file=file.path(proj$obs$dir, "qobs.dat"), index.name="time", sep="\t")
    write.zoo(proj$obs$rain, file=file.path(proj$obs$dir, "rain.dat"), index.name="time", sep="\t")
    write.zoo(proj$obs$pe, file=file.path(proj$obs$dir, "pe.dat"), index.name="time", sep="\t")
  }
  if(save.disc)
  {
	  # all the discretisations
	  proj$disc<-lapply(proj$disc,
	         function(disc)
	         {
	            disc$dir <- disc.dir.name(disc$cuts, dn=proj$disc.dir)

	            write.disc(disc)
	            return(disc)
	         }
	  )
  }
#  if(copy.source)
#  {
    try(writeRaster(proj$dem, file.path(dir, "dem.tif")))
    # channel proportions should be saved with the discretisation as they depend on the channel width specified
    # locations are dicretisation independent
    try(writeRaster(proj$drn.rast, file.path(dir, "reaches.tif")))
    try(writeOGR(proj$drn, dir, "drn", driver="ESRI Shapefile"))
    try(writeOGR(proj$sites, dir, "sites", driver="ESRI Shapefile"), silent=T)
 # }
  #
  return(proj)
}

# specify the location of dem, drn etc and discretisation separately
proj.from.dirs <- function(disc.dir,
                           root.dir=NULL, donor=NULL, obs.dir=NULL)
{

  proj <- create.proj(fn=NULL, data.dir = disc.dir, root.dir=root.dir)

  if(!is.null(donor))
  {
    #proj$groups <- donor

    proj$groups[,par.names] <-  donor$groups[nrow(donor$groups),par.names] #), nrow(proj$groups))
    proj$run.par <- donor$run.par
    proj$disp.par <- donor$disp.par
    proj$dt <- donor$dt
    proj$ntt <- donor$ntt
    proj$rain <- donor$rain
    proj$qt0 <- donor$qt0
    proj$qmax <- donor$qmax
  }
  # look for data
  return(proj)
}

# new, blanK project
new.project <- function()
{
  res <- environment()
  res$disp.par <- disp.par()
  res$run.par <- def.run.par()
  res$hsu.par<- def.hsu.par()
  res$description = "Dynamic TOPMODEL project file"
  res$run.title <-""
  res$disp.title <- ""
  res$notes <- ""
  res$dt<-1
  res$ntt<-1
  return(res)
}


# load a project file containing information required for dynamic TOPMODEl run
# and check that it's consistent (DEM same size and resolution, weighting and routing matrices consistent etc)
create.proj.old <-function(fn,
                           data.dir=NULL,  # set this to set the data dir, or read from file
                           root.dir=NULL,
                           read.obs=T)
{
  #proj <- new.project()
  # copy values form the project to current environment
  # ensure options$stringsAsFactors is set to F
  # default values - can be overwritten when loading project file
  # default time format - system default - CEH format is 10-sep-2001 00:00:00
  time.fmt = "%Y-%m-%d %H:%M:%S"
  dt <- 1
  cm <-NULL
  ntt <- 1
  vof <- 10
  vchan <- 150
  #  graphics.interval <- 10
  #  graphics.spatial.interval <- 10
  qobs <- NULL
  rain <- NULL
  drn <- NULL
  routing <- NULL
  weights <- NULL

  routing.fn <- NULL
  # parameters for parsing rain data
  rain.sep<-"\t"
  rain.metadata <- 0
  # do data contain col headings?
  rain.header <- T
  rain.fn <- NULL
  rain.header <- T
  rain.fn <- NULL

  # column indexes for time and values columns, if specified
  rain.tm <- 1  #"HOUR_ENDED"
  rain.val <- 2  #"RAINFALL"
  pe.tm <- 1
  pe.val <- 2
  pe.sep <- "\t"
  qobs.tm <- 1
  qobs.val <- 2
  qobs.sep <- "\t"

  # these can be overwritten
  rain.fmt <- time.fmt
  qobs.fmt <- time.fmt
  pe.fmt <- time.fmt

  dem.fn <- "dem.tif"
  cm.fn <- "cm.tif"
  dists.fn <- "dists.tif"
  groups.fn <- "groups.dat"
  weights.fn <- "weights.dat"
  ichan<-NULL
  drn.ln <- "drn"
  #   sim.delay=10
  start <- "2001-01-01"
  end <- "2001-01-31"
  disp.par <- disp.par()
  run.par <- def.run.par()
  hsu.par <- def.hsu.par()

  run.title <-""
  disp.title <- ""
  description = "Dynamic TOPMODEL project file"
  notes <- ""

  # probably total discharge, converted to specific within programme
  qobs.metadata<- 14
  qobs.fn <- ""
  qobs.sep <-"\t"
  qobs.header <- F
  qobs.tm <- 1  #"HOUR_ENDED"
  qobs.val <- 2  #"RAINFALL"

  # initial specific discharge m/hr. small non zero value
  qt0 <- 1e-4

  # evapotranspiration: limits or calculated
  pe.min <- 0
  pe.fn <- ""
  pe.max <- 8.4
  time.zone <- "GMT"
  pe <- NULL

  # qobs.fmt <- time.fmt
  routing.fn <- "routing.dat"

  #  id <- ""
  # log output defaults to console if empty, otherwise supply a valid file name
  # log.out <- "dtm.log"
  text.out <- ""
  if(!is.null(fn))
  {
    cat("Opening project file ", fn, "....\n")
    setwd(dirname(fn))
    # quickest way is to simply use the file as a  R source file - that way comments etc will be ignored
    # can use var=val syntax too
    source(fn, keep.source=F, local=TRUE)
    src.fn <- fn
    # note that data.dir specified here wil overwrite parameter
  }

  #else{return(environment())}
  #   if(!is.null(d.dir))
  #   {
  #     data.dir <- d.dir
  #   }

  # use default date time format. UTC unless time zone specified
  run.par$"start" <- as.POSIXct(start,tz=time.zone)  # start time / date (default: start of rain series)
  run.par$"end" <- as.POSIXct(end, tz=time.zone)

  # no more file information available so return
  if(is.null(fn) & is.null(data.dir)){return(environment())}
  save.dir <- getwd()

  if(exists.not.null("data.dir"))
  {
    setwd(data.dir)
    #  data.dir <- fp(getwd(), data.dir)
  }
  else if(exists("data.dir"))
  {
    #
    if(file.exists(file.path(getwd(), data.dir)))
    {
      data.dir <- file.path(getwd(), data.dir)
    }
    else if(!file.exists(data.dir))
    {
      # assume relative paths are wrt current dir
      cat("Specified data dir ", data.dir, " not found, trying to create...\n")

      # dir.create(data.dir, recursive=T)

    }
    setwd(data.dir)
  }
  else
  {
    # if unspecified, assumme everything relative to file directory
    data.dir <- getwd()
  }
  # default location for DEM, DRN etc
  if(!exists.not.null("root.dir")){
    root.dir <- file.path(getwd(), "..")
  }
  # restore current working directory howver we exit
  on.exit(setwd(save.dir))

  # copy current environment vals to matching values in the display and run
  # options list.Flag missing values?
  run.par <- CopyVals(from=environment(), to=run.par)
  # build run parameters
  disp.par <- CopyVals(from=environment(), to=disp.par)
  setwd(root.dir)

  if(exists.not.null("dem.fn", warn="DEM not found or not specified"))
  {
    dem <- try(raster(dem.fn))
    run.par$dem <- dem
    catch.area <- sum(!is.na(dem[]))*xres(dem)*yres(dem)
  }

  # discretisation specific data
  setwd(data.dir)

  if(exists.not.null("cm.fn", warn = "HSU classification raster not specfied or not found"))
  {
    cm <- try(raster(cm.fn))
    # add topogrphic data
    run.par$cm <- cm
  }

  # can supply routing explictly
  if(is.null(routing) & exists.not.null("routing.fn"))
  {
    routing <- try(read.table(routing.fn, header=T))
    # determine channels as those groups not identified with areas on raster
    if(is.null(ichan))
    {
      warning("Channel numbers inferred from routing table")
      ichan <- 1:nrow(matrix(routing, ncol=2))
      cat("Channels identified ", ichan, "\n")
      #
    }
  }
  else
  {
    # supply routing explicitly - if proportion not specified then assume all (1)
    if(length(routing)==1){routing<-c(routing, 1)}
    routing <- matrix(routing, ncol=2, byrow=T)
    # can supply this explicitly vector of channel lens
  }

  if(exists.not.null("groups.fn", warn="Areal grouping table not found"))
  {
    # copy defaults specified for this project
    proj.hsu.par <- CopyVals(environment(), def.hsu.par)

    groups <- try(read.table(groups.fn,
                             header=T))
    # add missing params using defaults

    new.groups<-data.frame()
    for(igroup in 1:nrow(groups))
    {
      pars <- as.list(groups[igroup,])
      new.groups<- rbind(new.groups,
                         as.vector(merge.lists(pars, proj.hsu.par)))
    }
    row.names(new.groups)<-new.groups$id
    groups <- new.groups

    #   	if(is.null(groups$dchan))
    #   	{
    #   		# average distance to channel
    #   		if(!file.exists(dists.fn)){stop("channel distances raster required")}
    #   		dists <- raster(dists.fn)
    #   		distsz <-zonal(dists, cm)
    #   		groups <- cbind(groups, dchan=rep(0, nrow(groups)))
    #   		#   groups[1:nrow(groups),]$atb.bar <- 0  # default value for channels
    #   		groups[match(distsz[,1], groups$id), ]$dchan <- distsz[,2]
    #
    #   	}

    # sigma and mean atb
    groups <- AddUpslopeAreas(groups, dem, cm)
  }

  # weighting matrix can be calulated from DEM and classification
  if(exists.not.null("weights.fn"))
  {
    cat("Reading flow transition matrix...\n")
    weights <- as.matrix(read.table(weights.fn))
    # quick check
    if(ncol(weights)!=nrow(groups))
    {
      stop("Group and weighting matrix dimension mismatch")
    }
  }
  else
  {
    # create as new from classification raster
    #  message("Recreating flow transition matrix....")
    #  weights <-make.flow.dist.matrix(dem,cm)
  }
  setwd(root.dir)
  # expect location of rain etc to be relative to dems etc
  if(exists.not.null("rain.fn", warn="Rainfall data nor specified or file not located"))
  {
    cat("Reading rainfall data...\n")
    # expect rain, pe and  in format date - amount variable name indicated
    rain <- read.obs(fn=rain.fn,
                     tm=rain.tm,
                     val=rain.val,
                     sep=rain.sep,
                     fmt=rain.fmt,
                     metalines=rain.metadata,
                     units="m",
                     header=rain.header)

    # convert to m / hr, if required
    if(max(rain, na.rm=T)>1)
    {
      message("Note: rainfall converted to m/hr")
      rain<-rain/1000
    }
  }

  if(read.obs)
  {
    # discharges should have same. If not supplied then create a dummy distribution
    # debugonce("read.obsCEH")
    if(exists.not.null("qobs.fn"))
    {
      cat("Reading discharge data from ", qobs.fn, "...\n")
      # rain, pe, temp etc livein "metobs" subfolder of data dir, flows in  "Qobs"
      qobs <- read.obs(qobs.fn, tm=qobs.tm, val=qobs.val,
                       fmt = qobs.fmt,
                       # DON'T NEED TO specify the columns as in expected order, but cmake sure
                       metalines=qobs.metadata,
                       sep=qobs.sep)
      names(qobs)<-NULL
      # apply conversion function to e.g convert
      if(exists("conv.qobs"))
      {
        qobs <- conv.qobs(qobs, catch.area)
      }
      # set initial discharge if nothing specified
      if(qt0<=0){qt0<-as.numeric(qobs[1,])}
    }

    if(exists.not.null("pe.fn"))
    {
      cat("Loading p.e data from ", pe.fn, "...\n")
      # pe can be supplied or calculated. here a file is supplied
      pe <- read.obs(pe.fn,
                     tm=pe.tm, val=pe.val,
                     fmt = time.fmt,
                     metalines=0,
                     sep=pe.sep)
      #     else
      #     {
      #       cat("Calculating potential evaporation from supplied daily min and max totals...\n")
      #       # create series from min and max
      #       # estimate sinusoidal pe given annual daily range.
      #       pe <- GetEvapSeries(start,
      #                           end,
      #                           dt, # time step, hours in hours
      #                           pe.min,
      #                           pe.max) #/1000 # pe in m/hr
      #     }

      # convert to m / hr, if required
      if(max(pe, na.rm=T)>0.5)
      {
        message("Note: P.E. converted to m/hr")
        pe<-pe/1000
      }
    }
  }
  try(
    drn <- readOGR(root.dir, drn.ln))

  # run.par$drn <- drn
  cat("\n", description, "\n")
  cat(paste0(notes, "\n"))
  return(environment())

}

# update contents of project folder, for examplel groups info, and
SaveProjectContents <- function(proj, dest=NULL)
{
  if(!is.null(dest))
  {
    if(!file.exists(dest))
    {
      dir.create(dest, recursive=T)
    }
    file.copy(file.path(proj$data.dir, proj$dem.fn),
              file.path(dest, proj$dem.fn))
    file.copy(file.path(proj$data.dir, proj$dem.cn),
              file.path(dest, proj$cm.fn))

    proj$data.dir<-dest

  }
  dsn <- proj$data.dir
  write.table(proj$groups, file.path(proj$data.dir , "groups.dat"))

  # dump(ls(proj), envir=proj, control=NULL, file=fn)
}

# expand the project start and end run date-times to include the rain observations

check.time.intervals <- function(proj)
{
  int <- time.interval.intersection(proj$obs$rain, proj$sim.start, proj$sim.end)
  return(length(int)>0)
}

time.interval.intersection <- function(obs, sim.start, sim.end)
{
  int <- which(as.numeric(index(obs))<= as.numeric(sim.end) &
                 as.numeric(index(obs))>= as.numeric(sim.start))
  if(length(int)>0)
  {
    return(index(obs)[int])
  }
  return(NULL)
}

# given a set of observations and a specified run interval
# expand / contrcat the run times to accommodate the observations
fix.run.dates <- function(proj)
{
  obs <- proj$obs$rain

  if(!is.null(obs) & !check.time.intervals(proj))
  {

    warning("No rainfall data within specified run start/ end times: adjusting...")
    len.sim <- proj$sim.start-proj$sim.end

    start <- start(index(obs$rain))
    cat("Setting sim start to ", "\n")
    proj$sim.start <- start
    end <- min(start + len.sim, end(index(obs$rain)))
    cat("Setting sim end to", end, "\n")
    proj$sim.end <- as.POSIXct(end)
  }
  return(proj)
}

aggregate_observations <- function(proj)
{
  try(proj<-fix.run.dates(proj))
  obs <- proj$obs
  # check that the specified start end and end dates contain at least some rainfall data. Other data are
  # less important and take null defaults if not specified
  try(obs$pe <- disaggregate_xts(proj$obs$pe,
                                 ser.start=proj$sim.start,
                                 ser.end=proj$sim.end,
                                 dt=proj$dt, is.rate=T))
  try(obs$rain <- disaggregate_xts(proj$obs$rain,
                                   ser.start=proj$sim.start,
                                   ser.end=proj$sim.end,
                                   dt=proj$dt, is.rate=T))
  # note observed flows required in specific discharge m/hr
  try(obs$qobs <- disaggregate_xts(proj$obs$qobs,
                                   ser.start=proj$sim.start,
                                   ser.end=proj$sim.end,
                                   dt=proj$dt, is.rate=T))

  return(obs)
}

# simple output of results, allows the simulated values to be supplied separately



# plot of simulated discharges etc after
disp.run.results <- function(run,
                             qmax=NULL, legend=F,
                             title = "",
                             disp.par = disp.par(),
                             ...)
{
  qobs <- run$qobs
  layout(matrix(1))
  qmax <- max(run$qsim, run$qobs, na.rm=T)*1000*1.25
  qsim <- run
  pe <- run$evap
  par(family="serif")
  disp.par$legend.show <- legend
  disp.par$title.main<- title
  DisplayDischargeSelection(qsim=run$qsim, evap=run$evap[,"ae"], rain=run$rain, qobs= run$qobs,
                            qmax=qmax,disp.par=disp.par,...)
  #, run.par=run.par)


}

# plot of observed discharges
disp.obs <- function(proj,
                     qmax=NULL, legend=F,
                     title = "",
                     disp.par = disp.par(),
                     ...)
{
  layout(matrix(1))

  qmax <- max(proj$obs$qobs, na.rm=T)*1000*1.25

  par(family="serif")
  disp.par$legend.show <- legend
  disp.par$title.main<- title
  DisplayDischargeSelection(qsim=proj$obs$qobs, evap=proj$obs$pe, rain=proj$obs$rain, qobs= NULL,
                            qmax=qmax, disp.par=disp.par,...)
  #, run.par=run.par)


}


########################################################################
# read and / or create time series input for a dynamic TOPMODEL run
########################################################################
build.obs <- function(dn,               # directory to read from
											est.pe=F,         # if T then calculate
											pe.max=8,         # if pe not found the max daily evapotransporation in mm/day
											pe.min=0          # if not found the min daily PE
											)
{
	proj <- list()
	proj$obs <- list()

	proj <- load.obs(proj, obs.dir = dn)

	if(is.null(proj$obs$rain))
	{
		stop("No rainfall data supplied")
	}
	obs <- proj$obs
	est.pe <- est.pe | is.null(obs$pe)
	if(est.pe)
	{
		# create observations using rainfall data as their basis
		tms <- index(obs$rain)
		dt <- modal(diff(tms))/60
		# create evap, note in m/hr
		obs$pe <- approx.pe.ts(start=tms[1], end=tms[length(tms)],
								 dt=dt,
								 emin=pe.min/1000,
								 emax=pe.max/1000)

	}
	return(obs)
}



# apply the given parameters to groups in all discretisations
apply.params <- function(proj, params, which=1:length(proj$disc))
{
  if(length(proj$disc)==0){return(proj)}
  params <- params[which(names(params) %in% colnames(proj$disc[[1]]$groups))]
  if(length(params)==0){return(proj)}

  proj$disc[which] <- lapply(proj$disc[which],
                             function(disc)
                             {

                               vals <- matrix(rep(unlist(params),nrow(disc$groups)), nrow=nrow(disc$groups), byrow=T)
                               disc$groups[,names(params)]<- vals

                               return(disc)
                             }
  )
  return(proj)
}


# graphics and text output
gr.on <- function(proj, spatial=F)
{
  return(graphics.on.off(proj,T, spatial))
}

gr.off <- function(proj, spatial=F)
{
  return(graphics.on.off(proj, F, spatial))
}

graphics.on.off <- function(proj, val=T, spatial=F)
{
  proj$disp.par$graphics.show <- val
  if(spatial)
  {
    proj$disp.par$graphics.spatial.show <- val

  }
  return(proj)
}

output.off <- function(proj)
{
  proj$disp.par$text.out <- NULL
  return(proj)
}

output.on <- function(proj)
{
  proj$disp.par$text <- stdout()

}




# plot a set of results associated with a topmodel run
plot.run.response <- function(run,   # dtm run output
                               qresp, # response series
                               evap=NULL,
                               ymax=NULL,
                                show.qobs=F,
                               fn=NULL,
                               lwd=1,
                               lty=1,
                               cols=c("blue", rainbow(ncol(qresp)-1)),
                               title="",
                              start=index(qresp)[1],
                              end=index(qresp)[length(index(qresp))],

                              ...)
{
  sel <- paste0(start, "::", end)
  run$qsim <- run$qsim[sel]
  qresp <- qresp[sel]
  if(length(evap)>0)
  {
    run$evap <- evap
    names(run$evap)<- "ae"
  }
  else
  {
    run$evap <- NULL
  }
  nresp<- ncol(qresp)
  if(!show.qobs)
  {
    run$qobs<-NULL
  }

  plot.run(run, qsim=qresp,
           cols=cols, fn=fn,
           title=title,
           ymax=ymax,
           lwd=lwd,
           lty=lty,
           start=start,
           end=end,
           #legend=F,
         #  legend.col=cols,
           ...)


}


# plot the results of a run using another matrix of q and a source project
plot.q <- function(proj,
                   qresp,
                   evap=NULL,
                   ymax=NULL,
                   fn=NULL,
                   lwd=2,
                   lty=1,
                   cols=rainbow(ncol(qresp)),
                   title="", ...)
{
  run <- proj
  run$disc <- NULL
  run$proj <- proj
  run$qobs <- NULL
  run$rain <- proj$obs$rain

  if(evap==F)
  {
    evap <- NULL
  }
  else if(is.null(evap) & !is.null(proj$obs$pe))
  {
    run$evap <- proj$obs$pe
    names(run$evap)<- "ae"
  }
  nresp<- ncol(qresp)

  plot.run(run, qsim=qresp*1000,
           cols=cols, fn=fn,
           title=title,
           ymax=ymax,
           lwd=lwd,
           lty=lty,
           ...)

}

get.calib.dir <- function(proj.ex)
{

  s <- format(proj.ex$sim.start, "s=%Y-%m-%d")
  e <- format(proj.ex$sim.end, "e=%Y-%m-%d")
  return(paste0(s, e, collapse=","))
}

# return the goodness of fit for a simulation run
gof.run <- function(run, pos=1) # specify column for multiple results
                    #s=),
                    #e=)
{
  sel <- paste0(first(index(run$qsim)), "::", last(index(run$qsim)))
  pos <- min(pos, ncol(run$obs$qobs))
  qobs <- run$qobs[sel, pos]
  if(!is.null(run$qobs) & nrow(qobs)==nrow(run$qsim))
  {

    res <- run.gof(run$qsim, qobs)
    #     res <- run.gof(as.vector(subset_zoo(run$qsim, s, e)),
    #                    as.vector(subset_zoo(run$qobs[,pos], s, e)))
    #
    return(res)
  }
}

