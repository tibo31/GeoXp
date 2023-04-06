histnbmap<- function(sf.obj, nb.obj, longlat = NULL, nbcol = 10, type = c("count", "percent", "density"),
                     sup = FALSE, criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16,
                     col = "lightblue3", xlab = "", ylab = "count", axes = FALSE, lablong = "", lablat = "") {
 
  ###################################################
  ########## COMMON to ALL FUNCTIONS in GeoXp
  
  envir <- globalenv()
  # Verification of the Spatial Object sf.obj
  class.obj <- class(sf.obj)[1]
  
  if(class.obj != "sf") 
    stop("sf.obj may be a sf object")
  
  # we propose to refind the same arguments used in first version of GeoXp
  if (st_geometry_type(sf.obj, by_geometry = F) %in% c("POINT"))
    my_coords <- st_coordinates(st_geometry(sf.obj))
  else
    my_coords <- st_coordinates(st_point_on_surface(st_geometry(sf.obj)))
  long <- my_coords[, 1]
  lat <- my_coords[, 2]
  
  listvar <- as.data.frame(st_drop_geometry(sf.obj))
  listnomvar <- colnames(listvar)
  
  # for identifying the selected sites
  if (!is.null(identify) && identify %in% colnames(sf.obj))
    label <- sf.obj[[identify]]
  else
    label <- ""
  
  nointer <- FALSE
  nocart <- FALSE
  buble <- FALSE
  z <- NULL
  legmap <- NULL
  legends <- list(FALSE, FALSE, "", "")
  labvar <- c(xlab, ylab)
  
  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix <- ""
  method <- ""
  listgraph <- c("Histogram", "Barplot", "Scatterplot")
  
  # Is there a Tk window already open ?
  if (interactive()) {
    if (!exists("GeoXp.open", envir = envir) ||
        length(ls(envir = .TkRoot$env, all.names = TRUE)) == 2) {
      assign("GeoXp.open", TRUE, envir = envir)
    } else {
      if (get("GeoXp.open", envir = envir)) {
        stop(
          "A GeoXp function is already open. 
          Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")
      } else {
        assign("GeoXp.open", TRUE, envir = envir)
      }
    }
  }
  
  # Windows device
  if(length(dev.list()) == 0 & options("device") == "RStudioGD")
    dev.new()
  # for graphic
  dev.new(noRStudioGD = FALSE)
  num_graph <- dev.list()[length(dev.list())]
  # for map
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]
  # number of devices
  num_supp <- NA
  
  #####################################################
  ##### Arguments proper to each function 
  #  objets spdep
  nb <- nb.obj
  W <- nb2mat(nb, zero.policy = TRUE)
  if (!inherits(nb, "nb"))
    stop("Not a neighbours list")
  c.nb <- card(nb)
  n.nb <- length(nb)
  regids <- attr(nb, "region.id")
  if (is.null(regids))
    regids <- as.character(1:n.nb)
  
  obs <- matrix(FALSE, nrow = n.nb, ncol = n.nb)
  graf <- "Neighbourplot1"
  
  # this code comes from a spdep function 
 
  # this code also ....
  if (sup == TRUE) {
    for (i in 1:n.nb) 
      nb[[i]] <- as.integer(nb[[i]][which(dlist[[i]] == max(dlist[[i]]))])
    W <- nb2mat(nb)   
  }
  
  dlist <- spdep::nbdists(nb=nb, coords = my_coords)
   
  ####################################################
  # selection d'un point sur la carte
  ####################################################
  
  pointfunc <- function() {
    if (graf == "Neighbourplot2") 
      SGfunc()
    
    quit <- FALSE
    graf <<- "Neighbourplot1"
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_carte)
      
      if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
      }
      
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
              lablong = lablong, lablat = lablat, label = label, symbol = pch,
              method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
              legmap = legmap, legends = legends, buble = buble, criteria = criteria,
              nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
        next
      }
      
      if (nrow(sf.obj) > 100 | st_geometry_type(sf.obj, by_geometry = F) == "POINT")
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                          Xpoly = loc[1], Ypoly = loc[2], method = "point")
      else {
        my_points <- st_as_sf(data.frame(x = loc$x, y = loc$y), coords = c("x", "y"),
                              crs = st_crs(sf.obj))
        for (i in 1:length(long)) {
          if (st_intersects(my_points, sf.obj[i, ], sparse = FALSE)) {
            obs[i, ] <<- !obs[i, ] 
            break
          }
        }
      }
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            lablong = lablong, lablat = lablat, label = label, symbol = pch,
            method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
            legmap = legmap, legends = legends, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
                bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
                labvar = labvar, symbol = pch, couleurs = col)
      
    }
  }
  
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunc <- function() {
    if (graf=="Neighbourplot2") 
      SGfunc()
    
    graf <<- "Neighbourplot1"
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    points(long, lat, pch = 16, col = "royalblue")
    
    while(!quit) {
      dev.set(num_carte)
      loc <- locator(1)
      if (is.null(loc)) {
        quit<-TRUE
        next
      }
      
      polyX <- c(polyX, loc[1])
      polyY <- c(polyY, loc[2])
      lines(polyX,polyY)
    }
    
    polyX <- c(polyX, polyX[1])
    polyY <- c(polyY, polyY[1])
    if (length(polyX) > 0) {
      lines(polyX, polyY)
      
      obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                        Xpoly = polyX, Ypoly = polyY, method = "poly")
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            lablong = lablong, lablat = lablat, label = label, symbol = pch,
            method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
            legmap = legmap, legends = legends, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
      
      graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
                bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
                labvar = labvar, symbol = pch, couleurs = col)
    }
  }
  
  ####################################################
  # selection d'une barre de l'histogramme
  ####################################################
  
  barfunc <- function() {
    if (graf=="Neighbourplot1") SGfunc()
    
    graf <<- "Neighbourplot2"
    quit <- FALSE
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_graph)
      loc <- locator(1)
      if (is.null(loc)) {
        quit<-TRUE
        graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
                  bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
                  labvar = labvar, symbol = pch, couleurs = col)
        next
      }
      obs <<- selectstat(var1 = nb, var2 = dlist, obs = obs,
                         Xpoly = loc[1], Ypoly = loc[2], method = "nbhist", nbcol = nbcol)
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            lablong = lablong, lablat = lablat, label = label, symbol = pch,
            method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
            legmap = legmap, legends = legends, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
      
      graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
                bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
                labvar = labvar, symbol = pch, couleurs = col)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
    }
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  
  cartfunc <- function()  { 
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            lablong = lablong, lablat = lablat, label = label, symbol = pch,
            method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
            legmap = legmap, legends = legends, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
    } else {
      tkmessageBox(message = "Spatial contours have not been given",
                   icon = "warning", type = "ok")    
    }
  }
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer <- function() {
    if (length(criteria) != 0) {
      nointer <<- !nointer
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            lablong = lablong, lablat = lablat, label = label, symbol = pch,
            method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
            legmap = legmap, legends = legends, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
    } else {
      tkmessageBox(message = "Criteria has not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function() {
    
    res2 <- choix.bubble(buble, listvar, listnomvar, legends, num_graph, num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          lablong = lablong, lablat = lablat, label = label, symbol = pch,
          method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
          legmap = legmap, legends = legends, buble = buble, criteria = criteria,
          nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- matrix(FALSE, nrow = length(long), ncol = length(long))
    
    # graphiques
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          lablong = lablong, lablat = lablat, label = label, symbol = pch,
          method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
          legmap = legmap, legends = legends, buble = buble, criteria = criteria,
          nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
    
    graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
              bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
              labvar = labvar, symbol = pch, couleurs = col)
  }
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function()  {
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    dev.off(num_graph)
    dev.off(num_carte)
    if(!is.na(num_supp))
      dev.off(num_supp)
  }
  
  quitfunc2 <- function() {
    fig_save <- "fig_GeoXp.pdf"
    map_save <- "map_GeoXp.pdf"
    k <- 1
    while(file.exists(fig_save)) {
      fig_save <- paste("fig_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }
    pdf(fig_save)
    graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
              bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
              labvar = labvar, symbol = pch, couleurs = col)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          lablong = lablong, lablat = lablat, label = label, symbol = pch,
          method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
          legmap = legmap, legends = legends, buble = buble, criteria = criteria,
          nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
    dev.off()
    
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure has been saved in", fig_save, "\n")
    assign("last.select", which(obs), envir = envir)
    
    dev.off(num_carte)
    dev.off(num_graph)
    if(!is.na(num_supp))
      dev.off(num_supp)
  }
  
  ####################################################
  # Graphique de base
  ####################################################

  carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            lablong = lablong, lablat = lablat, label = label, symbol = pch,
            method = "Neighbourplot3", W = W, axis = axes, cex.lab = cex.lab, 
            legmap = legmap, legends = legends, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart)
      
  graphique(var1 = nb, var2 = dlist, obs = obs, num = num_graph, 
                bin = type, graph = "histo.nb", nbcol = nbcol, W = W,
                labvar = labvar, symbol = pch, couleurs = col)
 
  ####################################################
  # Boite de Dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "histnbmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pointfunc)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = polyfunc)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",
           fill = "x")
    
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the histogram", font = "Times 12",
                   foreground = "darkred", background = "white"))
    barre.but <- tkbutton(frame1c, text = "Select a bar  ", command = barfunc)
    tkpack(barre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1c, expand = "TRUE", fill = "x")
    
    
    frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nettoy.but <- tkbutton(frame1b, text = "     Reset selection     " , command = SGfunc)
    tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1b, expand = "TRUE", fill = "x")
    
    
    frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2, text = "Options", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame2, text = "Spatial contours  ", font = "Times 11",
                   foreground = "darkred", background = "white"),
           tklabel(frame2, text = "Preselected sites  ", font = "Times 11", 
                   foreground = "darkred", background = "white"), 
           tklabel(frame2, text = "  Bubbles    ", font = "Times 11",
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = fbubble)
    tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
                   foreground = "blue", background = "white"))
    
    quit.but <- tkbutton(frame3, text = "Save results", command=quitfunc2)
    quit.but2 <- tkbutton(frame3, text = "Exit without saving", command = quitfunc)
    
    tkpack(quit.but, quit.but2, side = "left", expand = "TRUE", fill = "x")
    
    
    tkpack(frame3, expand = "TRUE", fill = "x")
  }
  ####################################################
  
  return(invisible())
}
