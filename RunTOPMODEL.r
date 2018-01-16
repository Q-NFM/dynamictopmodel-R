require(topmodel)

# utility routine that runs the topmodel implementation by Wouter Buytaert
run.topmodel <- function(dtm.parms, qt0, dt, routing,
                         rain, evap,
                         sim.start,
                         sim.end,
                         atb, chans, nclasses=10)
{
  params <- list()
  params$qs0 <- qt0
  params$lnTe <-  dtm.parms$ln_t0
  params$m <- dtm.parms$m
  params$SrMax <- dtm.parms$srz_max
  params$Sr0 <- dtm.parms$srz0*params$SrMax
  params$td <- dtm.parms$td
  params$vch <- dtm.parms$vchan
  params$vr <- dtm.parms$vchan
  params$k0 <- 100   # for infiltration excess, none so set to large value
  params$CD <- 7.235573e-01
  params$dt <- dt

  delays <- routing
  delays[,2]<- cumsum(delays[,2])

  atb[which(chans[[1]][]>0)]<- NA
  atb.cuts<- make.classes(atb[],nclasses)
  sel <- paste0(sim.start, "::", sim.end)
  rain.tm <- rain[sel][]
  evap.tm <- evap[sel][]

  qsim <- topmodel(unlist(params), atb.cuts, delay = delays,
           rain=rain.tm, ET0=evap.tm)

  qsim <- xts(qsim, order.by=index(rain[sel]))

  return(qsim)
}

# Run a multi-parameter set
calib.topmodel <- function(dtm.parms, qt0, dt, routing,
                         rain, evap,
                         sim.start,
                         sim.end,
                         atb, chans, nclasses=10, qobs=NA)
{
  params <- expand.grid(qt0, dtm.parms$ln_t0, dtm.parms$m, dtm.parms$srz_max,
              dtm.parms$srz0, dtm.parms$td, dtm.parms$vchan, dtm.parms$vchan, 100,  0.7, dt)

  params[,5]<- params[,4]*params[,5]

  delays <- routing
  delays[,2]<- cumsum(delays[,2])

  atb[which(chans[[1]][]>0)]<- NA
  atb.cuts<- make.classes(atb[],nclasses)
  sel <- paste0(sim.start, "::", sim.end)
  rain.tm <- as.vector(rain[sel][])
  evap.tm <- as.vector(evap[sel][])
  rain.tm[which(is.na(rain.tm))]<-0
  if(!is.na(qobs))
  {
    qobs.tm <- as.vector(qobs[sel][])
    params[,1] <- as.numeric(obs$qobs[sim.start][1])
  }

  cat("Running ", nrow(params), " sets\n")
  run <- topmodel(as.matrix(params), atb.cuts, delay = delays,
                   rain=rain.tm, ET0=evap.tm, verbose = F, Qobs=qobs.tm)

  qsim <- xts(run$Q, order.by=index(rain[sel]))

  return(list(qsim))
}

