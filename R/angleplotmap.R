angleplotmap <- function(sf.obj, name.var, quantiles = TRUE, 
                         criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, 
                         pch = 16, col = "lightblue3", xlab = "angle", ylab = "absolute magnitude",
                         axes = FALSE, lablong = "", lablat = "") {
  
  ###################################################
  ########## COMMON to ALL FUNCTIONS in GeoXp
  
  envir <- globalenv()
  # Verification of the Spatial Object sf.obj
  class.obj <- class(sf.obj)[1]
  
  if(class.obj != "sf") 
    stop("sf.obj may be a sf object")
  
  if(!(name.var %in% names(sf.obj)))
    stop("name.var is not included in the sf object")
  
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
  
  var <- sf.obj[[name.var]]
  obs <- matrix(FALSE, nrow = length(long), ncol = length(long))
  names.slide <- "Alpha Quantile Value"
  
  # calcul des matrices theta et absvar
  long1 <- matrix(rep(t(long), length(long)), ncol = dim(t(long))[2],
                  byrow = FALSE)
  long2 <- matrix(rep(t(long), length(long)), ncol = dim(t(long))[2],
                  byrow = TRUE)
  lat1 <- matrix(rep(t(lat), length(lat)), ncol = dim(t(lat))[2],
                 byrow = FALSE)
  lat2 <- matrix(rep(t(lat), length(lat)), ncol = dim(t(lat))[2],
                 byrow = TRUE)
  theta <- matrix(0, nrow = length(long), ncol = length(long))
  numer <- lat2 - lat1
  denom <- long2 - long1
  
  for (i in 1:length(long)) {
    for (j in 1:length(long)) {
      if (denom[i, j] == 0) 
        theta[i, j] <- pi/2
      else 
        theta[i, j] <- atan(numer[i, j]/denom[i, j])
    }
  }
  
  theta[which(theta < 0)] <- theta[which(theta < 0)] + pi
  v1 <- matrix(rep(t(var), length(var)), ncol = dim(t(var))[2],
               byrow = FALSE)
  v2 <- matrix(rep(t(var), length(var)), ncol = dim(t(var))[2],
               byrow = TRUE)
  absvar <- abs(v1 - v2)
  
  ####################################################
  #choix des bornes 
  ####################################################
  borne1 <- 0.01
  borne2 <- 0.99
  alpha <- 0.5
  
  temp_data <- data.frame(
    var1 = sort(theta),
    var2 = absvar[order(theta)])
  # alpha1 <- qgam::qgam(var2 ~ s(var1, k = 20, bs = "ad"), 
  #                     data = temp_data, qu = alpha)
  alpha1 <- quantreg::rq(var2 ~ var1, data = temp_data, tau = alpha)
  
  ####################################################
  # selection d'un point sur l'angleplot
  ####################################################
  
  pointfunc <- function() {
    
    quit <- FALSE
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_graph)
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                  graph = "Angleplot", labvar = labvar, couleurs = col, 
                  symbol = pch, quantiles = quantiles, alpha1 = alpha1)
        next
      }
      
      obs <<- selectstat(var1 = theta, var2 = absvar, obs = obs, Xpoly = loc[1], Ypoly = loc[2],
                         method = "AnglePoint", long = long, lat = lat)
      
      diag(obs) <<- FALSE
      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                graph = "Angleplot", labvar = labvar, couleurs = col, 
                symbol = pch, quantiles = quantiles, alpha1 = alpha1)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, 
            criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label,
            cex.lab = cex.lab, symbol = pch, method = "Angleplot",
            axis = axes, legmap = legmap, legends = legends) 
    }
  }
  
  ####################################################
  # selection d'un polygone sur l'angleplot
  ####################################################
  
  polyfunc <- function() {
    quit <- FALSE
    polyX <- NULL
    polyY <- NULL
    
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_graph)
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        next
      }
      
      polyX <- c(polyX, loc[1])
      polyY <- c(polyY, loc[2])
      
      if (length(polyX) > 0)
        lines(polyX, polyY)
    }  
    
    polyX <- c(polyX, polyX[1])
    polyY <- c(polyY, polyY[1])
    
    if (length(polyX) > 0) {
      lines(polyX, polyY)
      for (i in 1:length(long)) {

        my_points <- st_as_sf(data.frame(x = as.numeric(theta[, i]), 
                                         y = as.numeric(absvar[, i])), 
                                 coords = c("x", "y"))
        polyg <- cbind(unlist(polyX), unlist(polyY))
        pol <- st_sfc(st_polygon(list(polyg)))
        def <- as.vector(st_intersects(my_points, pol, sparse = FALSE))
        obs[def, i] <<- !obs[def, i]    
      }
      
      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                graph = "Angleplot", labvar = labvar, couleurs = col, 
                symbol = pch, quantiles = quantiles, alpha1 = alpha1)
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, 
            criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label,
            cex.lab = cex.lab, symbol = pch, method = "Angleplot",
            axis = axes, legmap = legmap, legends = legends) 
    }
  }    
  
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() { 
    if (length(carte) != 0) {
      nocart <<- !nocart
        
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, 
            criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label,
            cex.lab = cex.lab, symbol = pch, method = "Angleplot",
            axis = axes, legmap = legmap, legends = legends) 
    } else {
      tkmessageBox(message = "Spatial contours have not been given",
                   icon = "warning", type = "ok")    
    }
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc<-function() {
    
    obs <<- matrix(FALSE, nrow = length(long), ncol = length(long))
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, 
          criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label,
          cex.lab = cex.lab, symbol = pch, method = "Angleplot",
          axis = axes, legmap = legmap, legends = legends) 
    
    graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
              graph = "Angleplot", labvar = labvar, couleurs = col, 
              symbol = pch, quantiles = quantiles, alpha1 = alpha1)
  }
  
  
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function() {
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    dev.off(num_graph)
    dev.off(num_carte)
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
    graphique(var1 = theta, var2 = absvar, obs = obs, num = dev.list()[length(dev.list())], 
              graph = "Angleplot", labvar = labvar, couleurs = col, 
              symbol = pch, quantiles = quantiles, alpha1 = alpha1)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())], buble = buble, 
          criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label,
          cex.lab = cex.lab, symbol = pch, method = "Angleplot",
          axis = axes, legmap = legmap, legends = legends) 
    dev.off()
    
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure has been saved in", fig_save, "\n")
    assign("last.select", which(obs), envir = envir)
    
    dev.off(num_carte)
    dev.off(num_graph)
  }
  
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer <- function() {
    if (length(criteria) != 0) {
      nointer <<- !nointer
    
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj,  num = num_carte, buble = buble, 
            criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label,
            cex.lab = cex.lab, symbol = pch, method = "Angleplot",
            axis = axes, legmap = legmap, legends = legends) 
    } else {
      tkmessageBox(message = "Criteria has not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble<-function()  {
    
    res2 <- choix.bubble(buble, listvar, listnomvar, legends, num_graph, num_carte)
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, 
          criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label,
          cex.lab = cex.lab, symbol = pch, method = "Angleplot",
          axis = axes, legmap = legmap, legends = legends) 
  }
  
  
  ####################################################
  # Pour le alpha 
  ####################################################
  
  refresh.code <- function(...) {
    
    alpha <<- slider1(names.slide = names.slide, no = 1)
  
    temp_data <- data.frame(
      var1 = sort(theta),
      var2 = absvar[order(theta)])
    
    # alpha1 <<- qgam::qgam(var2 ~ s(var1, k = 20, bs = "ad"), data = temp_data, qu = alpha)
    
    alpha1 <<- quantreg::rq(var2 ~ var1, data = temp_data, tau = alpha)
    
    graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
              graph = "Angleplot", labvar = labvar, couleurs = col, 
              symbol = pch, quantiles = quantiles, alpha1 = alpha1)
  }
  
  ####################################################
  # Representation graphique
  ####################################################
  
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, 
            criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label,
            cex.lab = cex.lab, symbol = pch, method = "Angleplot",
            axis = axes, legmap = legmap, legends = legends) 
      
      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                graph = "Angleplot", labvar = labvar, couleurs = col, 
                symbol = pch, quantiles = quantiles, alpha1 = alpha1)
      
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if(interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "angleplotmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the graph", font = "Times 12",
                   foreground = "darkred", background = "white"))
    
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pointfunc)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = polyfunc)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",fill = "x")
    
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
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
                   foreground = "darkred", background = "white"),side = "left", 
           fill = "x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = fbubble)
    tkpack(nocou1.but, noint1.but, bubble.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    if(quantiles) {
      frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
      
      slider1(frame1c, refresh.code, names.slide = names.slide,
              borne1, borne2, (borne2 - borne1)/100, alpha)
      
      tkpack(frame1c, expand = "TRUE", fill = "x")
    
    }
    
    frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
                   foreground = "blue", background = "white"))
    
    quit.but <- tkbutton(frame3, text = "Save results", command = quitfunc2)
    quit.but2 <- tkbutton(frame3, text = "Exit without saving", command = quitfunc)
    
    tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
           fill = "x")
    
    tkpack(frame3, expand = "TRUE", fill = "x")
  }
  ####################################################
  return(invisible())
  
  # return(obs)
}

