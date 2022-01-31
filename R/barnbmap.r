barnbmap <- function(sp.obj, nb.obj, criteria = NULL, carte = NULL, identify = FALSE, 
                     cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", 
                     axes = FALSE, lablong = "", lablat = "") {
  
  # Verification of the Spatial Object sp.obj
  class.obj <- class(sp.obj)[1]
  spdf <- (class.obj == "SpatialPolygonsDataFrame")
  
  if(substr(class.obj, 1, 7) != "Spatial") 
    stop("sp.obj may be a Spatial object")
  
  # Initialisation des objets de spdep
  nb <- nb.obj
  W <- nb2mat(nb,zero.policy=TRUE)
  if (!inherits(nb, "nb")) 
    stop("Not a neighbours list")
  
  # we propose to refind the same arguments used in first version of GeoXp
  coords <- coordinates(sp.obj)
  
  if(substr(class.obj, nchar(class.obj) - 8, nchar(class.obj)) == "DataFrame") {
    listvar <- sp.obj@data
    listnomvar <- names(sp.obj@data)
  } else {
    listvar <- NULL
    listnomvar <- NULL
  }
  
  # Code which was necessary in the previous version
  object <- nb.obj
  
  # for identifyng the selected sites
  if(identify)
    label <- row.names(listvar)
  else
    label <- ""
  
  # initialisation
  nointer <- FALSE
  nocart <- FALSE
  buble <- FALSE
  z <- NULL
  legmap <- NULL
  legends <- list(FALSE, FALSE, "", "")
  labvar <- c(xlab, ylab)
  
  # Transformation data.frame en matrix
  if((length(listvar) > 0) && (dim(as.matrix(listvar))[2] == 1)) 
    listvar <- as.matrix(listvar)
  
  # Initialisation des objets de spdep
  c.nb <- card(nb)
  n.nb <- length(nb)
  regids <- attr(nb, "region.id")
  
  if (is.null(regids)) 
    regids <- as.character(1:n.nb)
  
  long <- coords[, 1]
  lat <- coords[, 2]
  
  obs <- matrix(FALSE, nrow = n.nb, ncol = n.nb)
  
  graf <- "Neighbourplot1"
  
  # Is there a Tk window already open ?
  if (interactive()) {
    if (!exists("GeoXp.open", envir = globalenv())) {
      assign("GeoXp.open", TRUE, envir = globalenv())
    } else {
      if (get("GeoXp.open", envir = globalenv())) {
        stop(
          "A GeoXp function is already open. 
          Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")
      } else {
        assign("GeoXp.open", TRUE, envir = globalenv())
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
  
  
  # for colors in map and new grahics
  if(length(col) == 1)
    col2 <- "blue"
  else 
    col2 <- col

  col3 <- "lightblue3"
  
  
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
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
              lablong = lablong, lablat = lablat, label = label, symbol = pch, 
              method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
              legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
              cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
              classe = card(object), cex.lab = cex.lab)
        next
      }
      
      if (!spdf | length(long) > 75) { 
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                          Xpoly = loc[1], Ypoly = loc[2], method = "point")
        } else {
          if (gContains(sp.obj, SpatialPoints(cbind(loc$x, loc$y),
                                              proj4string = CRS(proj4string(sp.obj))))) {
            for (i in 1:length(long)) {
              if (gContains(sp.obj[i, ], SpatialPoints(cbind(loc$x, loc$y), 
                                                       proj4string = CRS(proj4string(sp.obj))))) {
                obs[i, ] <<- !obs[i, ] 
                break
              }
            }
          } 
        }
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            lablong = lablong, lablat = lablat, label = label, symbol = pch, 
            method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
            legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
            cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
            classe = card(object), cex.lab = cex.lab)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      if(spdf & length(long) > 75 & !buble)
        points(long, lat, pch = 16, col = "royalblue")
      
      graphique(var1 = nb, obs = obs, num = num_graph, graph = "bar.nb", 
                W = W, labvar = labvar, symbol = pch, couleurs = col)
    }
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunc <- function() {
    if (graf == "Neighbourplot2")
      SGfunc()
  
    graf <<- "Neighbourplot1"
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    if(spdf) 
      points(long, lat, pch = 16, col = "royalblue")
    
    while (!quit) {
      dev.set(num_carte)
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
              lablong = lablong, lablat = lablat, label = label, symbol = pch, 
              method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
              legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
              cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
              classe = card(object), cex.lab = cex.lab)
        next
      }
      
      polyX <- c(polyX, loc[1])
      polyY <- c(polyY, loc[2])
      lines(polyX, polyY)
    }
    
    polyX <- c(polyX, polyX[1])
    polyY <- c(polyY, polyY[1])
    if (length(polyX) > 0) {
      lines(polyX, polyY)
      
      obs <<- selectmap(var1 = long, var2 = lat, obs = obs, Xpoly = polyX, 
                        Ypoly = polyY, method = "poly")
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            lablong = lablong, lablat = lablat, label = label, symbol = pch, 
            method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
            legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
            cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
            classe = card(object), cex.lab = cex.lab)
      
      graphique(var1 = nb, obs = obs, num = num_graph, graph = "bar.nb", 
                W = W, labvar = labvar, symbol = pch, couleurs = col)
    }
  }
  
  ####################################################
  # selection d'une barre de l'histogramme
  ####################################################
  
  barfunc <- function() {
    if (graf == "Neighbourplot1") 
      SGfunc()
    graf <<- "Neighbourplot2"
    
    quit <- FALSE
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while(!quit) {
      dev.set(num_graph)
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        graphique(var1 = nb, obs = obs, num = num_graph, graph = "bar.nb", 
                  W = W, labvar = labvar, symbol = pch, couleurs = col)
        next
      }
      obs <<- selectstat(var1 = nb, obs = obs, 
                         Xpoly = loc[1], Ypoly = loc[2],
                         method = "barnb")
      
      # graphiques
      graphique(var1 = nb, obs = obs, num = num_graph, graph = "bar.nb", 
                W = W, labvar = labvar, symbol = pch, couleurs = col)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            lablong = lablong, lablat = lablat, label = label, symbol = pch, 
            method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
            legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
            cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
            classe = card(object), cex.lab = cex.lab)
    }
  }
  
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() { 
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            lablong = lablong, lablat = lablat, label = label, symbol = pch, 
            method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
            legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
            cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
            classe = card(object), cex.lab = cex.lab)
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
      
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            lablong = lablong, lablat = lablat, label = label, symbol = pch, 
            method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
            legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
            cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
            classe = card(object), cex.lab = cex.lab)
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
    
    carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
          lablong = lablong, lablat = lablat, label = label, symbol = pch, 
          method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
          legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
          cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
          classe = card(object), cex.lab = cex.lab)
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- matrix(FALSE, nrow = length(long), 
                   ncol = length(long))
    
    # graphiques
    carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
          lablong = lablong, lablat = lablat, label = label, symbol = pch, 
          method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
          legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
          cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
          classe = card(object), cex.lab = cex.lab)
    
    graphique(var1 = nb, obs = obs, num = num_graph, graph = "bar.nb", 
              W = W, labvar = labvar, symbol = pch, couleurs = col)
  }
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function() {
    
    fig_save <- "fig_GeoXp.pdf"
    map_save <- "map_GeoXp.pdf"
    k <- 1
    while(file.exists(fig_save)) {
      fig_save <- paste("fig_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }
    pdf(fig_save)
    graphique(var1 = nb, obs = obs, num = dev.list()[length(dev.list())], graph = "bar.nb", 
              W = W, labvar = labvar, symbol = pch, couleurs = col)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = dev.list()[length(dev.list())],
          lablong = lablong, lablat = lablat, label = label, symbol = pch, 
          method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
          legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
          cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
          classe = card(object), cex.lab = cex.lab)
    dev.off()
    

    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = globalenv())
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure has been saved in", fig_save, "\n")

    assign("last.select", which(obs), envir = globalenv())
    
    dev.off(num_carte)
    dev.off(num_graph)
  }
  
  ####################################################
  # Representation des graphiques
  ####################################################

  carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
              lablong = lablong, lablat = lablat, label = label, symbol = pch, 
              method = "Neighbourplot1", W = W, axis = axes, legmap = legmap, 
              legends = legends, buble = buble, criteria = criteria, nointer = nointer, 
              cbuble = z, carte = carte, nocart = nocart, couleurs = col2, 
              classe = card(object), cex.lab = cex.lab)
        
  graphique(var1 = nb, obs = obs, num = num_graph, graph = "bar.nb", 
                  W = W, labvar = labvar, symbol = pch, couleurs = col)
        
  ####################################################
  # Boite de Dialogue
  ####################################################
  
  if(interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, 
                                 weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "barnbmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    
    point.but <- tkbutton(frame1a, text="Selection by point", command = pointfunc)
    poly.but <- tkbutton(frame1a, text="Selection by polygon ", command = polyfunc)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",
           fill = "x")
    
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the bar plot", font = "Times 12",
                   foreground = "darkred", background = "white"))
    barre.but <- tkbutton(frame1c, text = "Select a bar", command = barfunc)
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
                   foreground = "darkred", background = "white"),
           side = "left", fill="x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = fbubble)
    tkpack(nocou1.but, noint1.but, bubble.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
                   foreground = "blue", background = "white"))
    
    quit.but2 <- tkbutton(frame3, text = "Exit", command = quitfunc)
    
    tkpack(quit.but2, side = "left", expand = "TRUE",
           fill = "x")
    
    tkpack(frame3, expand = "TRUE", fill = "x")
  }
  ####################################################
  return(invisible())
  
}
