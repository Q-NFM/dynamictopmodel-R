# pre-calc elevation matrix with averaged values and x-y axis values for display
calc.persp.matrix <- function(run.par, disp.par)
{
  dem <- run.par$dem
  if(is.null(dem)){
    stop("DEM raster required for spatial output")
  }
  # these could also be supplied from input
  disp.par$x <- xres(dem)*0:(nrow(dem))
  disp.par$y <- yres(dem)*0:(ncol(dem))

  dem <- as.matrix(dem)
  # elevs are midpoints of cells - produce a matrix with averaged elevs
  # at edges of cells
  z <- dem
  z1 <-rbind(z[1,],z)
  z1 <-cbind(z1[,1],z1)
  z2 <-dem
  z2 <-rbind(z2, z2[nrow(z2),])
  z2 <-cbind(z2, z2[,ncol(z2)])

  disp.par$z <-(z1+z2)/2    # average

  cm <- round(run.par$cm)
  # map class names to order
  # substitute the saturation percents for the hsu ids
  # subst.tab <- cbind(groups, "psat"=pc)
  # doesn't really matter here as long as groups stay in thre right order
  class.ids <- unique(cm, na.rm=T)

  disp.par$cm.map <- as.matrix(subs(cm, data.frame(class.ids, 1:length(class.ids))))

  # plot the dem once to get the transformation matrix
  with(disp.par,
       pmat<-persp(x=x, y=y,shade=0.2,
                   z=z, ,
                   theta=graphics.spatial.theta,
                   expand= graphics.spatial.expand, phi=graphics.spatial.phi,
                   box=F,
                   border="gray")
  )

  if(!is.null(run.par$drn))
  {
    # transform the river channel network to 2d for overlay on plot
    #  message("projecting river network for display...")
    # disp.par$drn.dash < Get3dTransformedDRN(run.par$dem, run.par$drn, pmat)
  }

  return(disp.par)
}

disp.spatial <- function(groups,
                         vals,
                         cm,   # class matrix: dem not required - should be precomputed. Or could add to cm so only one disk access required
                         # disp.par contains the necessary info to show elevation perspective plot
                         # and parameters for plot display e.g. horz and vertical view angles and vert expansion factor
                         disp.par,  #
                         range.vals=range(vals),  # max sd (mm) for areas that are shaded (> is white)
                         time=NULL,  # simulation time
                         ichan=1,  # channel ids
                         nlevels=10,
                         ramp = colorRampPalette(c("white", "blue"), bias=1),
                         main="",  # "Storage deficits (mm)"
                         expand=disp.par$graphics.spatial.expand,
                         theta=disp.par$graphics.spatial.theta,
                         phi=disp.par$graphics.spatial.phi)
{
  groups <- groups[-ichan,]
  # remove channels
  if(length(vals==length(ichan)+length(groups))){
    vals <- vals[-ichan]
  }
  if(is.null(cm)){
    stop("Classification raster required for spatial output")
  }
  # get colours
  satRamp<- ramp(nlevels)
  min.val <- range.vals[1]
  max.val <- range.vals[2]
  range <- max.val - min.val
  # most interested in areas close to saturation (i.e top m) - display these in shades of blue
  # areas with larger deficits are white. any vals > max are set to max
  vals <- pmin(vals, max.val)
  vals <- pmax(vals, min.val)  # -ve values not allowed
  props <- nlevels*(vals-min.val)/range
  pc <- round(props+1)

  cols <- satRamp[pc]

  # determine facet colours for persp plot. fields::fields::drape.colors handles midpoints
  drape.cols <- fields::drape.color(zlim=range(cm[], na.rm=T), cm,
                                    midpoint=F, col=cols)$color.index

  border<-"gray"

  x <- disp.par$x
  y <- disp.par$y
  z <- disp.par$z
  if(length(z)>1e4)
  {
    # No cell  boundaries if large grid
    border <- NA
  }

  # pmat is affine transformation matrix to project 2d to 3d coordinates
  pmat<-persp(x=x, y=y,shade=0.1,
              z=z, col=drape.cols, theta=theta, expand=expand, phi=phi,
              box=F,
              border=border,
              #legend=F,
              d=0.75,
              main=main, cex.main=1)

  fields::image.plot(legend.only=T, col=satRamp, legend.shrink=0.5, legend.mar=3,
                     #   label.breaks = 0:(nlevs-1),
                     zlim=c(0, max.val), horizontal=T)

}


# splitscreen layout
get.layout <- function(disp.par)
{
  #    layout(matrix(c(1,1,1,1,1,1,2,2,3), nrow=3, byrow=TRUE))
  res <- matrix(c(1,1,1,1), nrow=4, byrow=TRUE)  #2

  show.spatial <-  F #disp.par$"graphics.spatial.show"  # & tm >= start

  if(show.spatial)
  {
    if(disp.par$graphics.spatial.window.id<1)
    {
      # layout so that discharges occupy top half, storages etc the bottom
      res <- matrix(c(1,1,1,2,3,3), nrow=2, byrow=TRUE)
    }
    else
    {
      # layout so that discharges occupy top half, storages etc the bottom
      res <- matrix(c(1,1,1,1,2,3), nrow=3, byrow=TRUE)
    }
  }
  return(res)
}

GetSpatialOutput <- function(groups, flows, stores, disp.par, max.q)
{
  if(disp.par$graphics.spatial.output=="qbf")
  {
    # expand to mm/hr
    vals <- 1000*flows$qbf
    range.vals <-  c(0, max.q)
    ramp <- colorRampPalette(c("white", "blue"))
    main <- "Specific base flows (mm/hr)"
    nlevels<-10
  }
  else if(disp.par$graphics.spatial.output=="surface.storage")
  {
    #  if(any(stores$ex>0)){browser()}
    # overland flow
    vals <- 1000*stores$ex#/groups$sd_max
    ramp <- colorRampPalette(c("white", "blue"))
    range.vals <-  c(0, 2)
    nlevels<-10
    main <- "Overland storage (mm)"
  }
  else
  {
    vals.sd <- 1000*stores$sd#/groups$sd_max
    # reversed colour map: close to saturation indicated by blue
    ramp.sd <- colorRampPalette(c("blue", "wheat", "white"))
    range.vals.sd <-  c(0, 150)

    #   vals <- vals.qbf
    vals <- vals.sd
    range.vals <- range.vals.sd
    ramp <- ramp.sd
    nlevels<-20
    main <- "Storage deficits (mm)"
  }
  return(list("vals"=vals, "main"=main, "nlevels"=nlevels, "ramp"=ramp, "range.vals"=range.vals))
}

SpatialSequence <- function(groups, flows, stores, ichan=1,
                            ae, rain, qobs, qsim,
                            s=start(qsim),
                            e=end(qsim),
                            disp.par, title=disp.par$title.main,
                            nframe=5)
{
  # split screen between hydrograph and spatial output (top)
  layout(matrix(c(1,1,1,1,2,2), nrow=3, byrow=TRUE))

  # time points corresponding to points at which frame is shown
  times <- seq(s, e, length.out=nframe)
  dx <- 1/nframe
  for(i in 1:nframe)
  {
    tm <- times[i]
    spatial.output <- GetSpatialOutput(groups, flows[tm], stores[tm], disp.par)
    # poistion frame
    SetFigWindow((i-1)*dx, i*dx, 2/3, 1)
    par(new=T)

    disp.spatial(groups, vals=spatial.output$vals,
                 nlevels=spatial.output$nlevels,
                 cm=disp.par$cm.map,
                 time=tm,
                 main= spatial.output$main,
                 range.vals = spatial.output$range.vals,
                 ichan=ichan,
                 ramp=spatial.output$ramp,
                 disp.par = disp.par)



  }

  disp.discharge.selection(disp.par,
                           qsim, ae, rain,
                           groups, qobs,
                           s, e,
                           ymax, ichan, title)



}
