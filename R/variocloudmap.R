variocloudmap <- function(sp.obj, name.var, bin = NULL, quantiles = TRUE, names.attr = names(sp.obj), criteria = NULL,
                          carte = NULL, identify = FALSE, cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", 
                          axes = FALSE, lablong = "", lablat = "", xlim = NULL, ylim = NULL) {

  envir <- as.environment(1)
  
  # Verification of the Spatial Object sp.obj
  class.obj <- class(sp.obj)[1]
  
  if (substr(class.obj, 1, 7) != "Spatial")
    stop("sp.obj may be a Spatial object")
  
  if (substr(class.obj, nchar(class.obj) - 8, nchar(class.obj)) != "DataFrame")
    stop("sp.obj should contain a data.frame")
  
  if (!is.numeric(name.var) &
      is.na(match(as.character(name.var), names(sp.obj))))
    stop("name.var is not included in the data.frame of sp.obj")
  
  if (length(names.attr) != length(names(sp.obj)))
    stop("names.attr should be a vector of character with a length equal to the number of variable")
  
  # we propose to refind the same arguments used in first version of GeoXp
  long <- coordinates(sp.obj)[, 1]
  lat <- coordinates(sp.obj)[, 2]
  
  var <- sp.obj@data[, name.var]
  
  # verify the type of the main variable
  if (!(is.integer(var) ||
        is.double(var)))
    stop("the variable name.var should be a numeric variable")
  
  
  listvar <- sp.obj@data
  listnomvar <- names.attr
  
  # Code which was necessary in the previous version
  if (is.null(carte) &
      class.obj == "SpatialPolygonsDataFrame")
    carte <- spdf2list(sp.obj)$poly
  
  # for identifying the selected sites
  if (identify)
    label <- row.names(listvar)
  else
    label <- ""
  
  # initialisation
  nointer <- FALSE
  nocart <- FALSE
  buble <- FALSE
  legends <- list(FALSE, FALSE, "", "")
  z <- NULL
  legmap <- NULL
  inout <- NULL
  labvar <- c(xlab, ylab)
  
  opt1 <- 1
  opt2 <- 1
  
  angle <- 0
  names.slide <- "Alpha Quantile Value"
  obs <- matrix(FALSE, nrow = length(long), ncol = length(long))
  
  directionnel <- FALSE
  
  # Is there a Tk window already open ?
  if (interactive()) {
    if (!exists("GeoXp.open", envir = baseenv()) ||
        length(ls(envir = .TkRoot$env, all.names = TRUE)) == 2) {
      assign("GeoXp.open", TRUE, envir = baseenv())
    } else {
      if (get("GeoXp.open", envir = baseenv())) {
        stop(
          "A GeoXp function is already open. 
          Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")
      } else {
        assign("GeoXp.open", TRUE, envir = baseenv())
      }
    }
  }
  
  # Windows device
  if (length(dev.list()) == 0 & options("device") == "RStudioGD")
    dev.new()
  # for graphic
  dev.new(noRStudioGD = FALSE)
  num_graph <- dev.list()[length(dev.list())]
  # for map
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]

  # Transformation data.frame en matrix
  if ((length(listvar) > 0) &&
      (dim(as.matrix(listvar))[2] == 1))
    listvar <- as.matrix(listvar)
  
  
  ####################################################
  # calcul des matrices diff et dist
  ####################################################
  
  long1 <- matrix(rep(t(long), length(long)), ncol = dim(t(long))[2], byrow = FALSE)
  long2 <- matrix(rep(t(long), length(long)), ncol = dim(t(long))[2], byrow = TRUE)
  
  lat1 <- matrix(rep(t(lat), length(lat)), ncol = dim(t(lat))[2], byrow = FALSE)
  lat2 <- matrix(rep(t(lat), length(lat)), ncol = dim(t(lat))[2], byrow = TRUE)
  
  v1 <- matrix(rep(t(var), length(var)), ncol = dim(t(var))[2], byrow = FALSE)
  v2 <- matrix(rep(t(var), length(var)), ncol = dim(t(var))[2], byrow = TRUE)
  
  theta <- matrix(0, nrow = length(long), ncol = length(long))
  numer <- lat2 - lat1
  denom <- long2 - long1
  
  theta[which(denom == 0, arr.ind = TRUE)] <- pi / 2
  theta[which(denom != 0, arr.ind = TRUE)] <-
    atan(numer[which(denom != 0, arr.ind = TRUE)] / denom[which(denom != 0, arr.ind =
                                                                  TRUE)])
  
  theta[which(theta < 0)] <- theta[which(theta < 0)] + pi
  
  dist <- sqrt((long1 - long2) ^ 2 + (lat1 - lat2) ^ 2)
  dif <-  (v1 - v2) ^ 2 / 2
  dif2 <-  (abs(v1 - v2)) ^ (1 / 2)
  
  borne1 <- 0.01
  borne2 <- 0.99
  alpha <- 0.5
  
  ####################################################
  # selection d'un point sur le variocloud
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
        graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
                  graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
                  alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)
        next
      }
      obs <<- selectstat(var1 = dist, var2 = dif, obs = obs,Xpoly = loc[1], Ypoly = loc[2], 
                         method = "Variopoint", long = long, lat = lat)
      
      graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
                graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
                alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
            label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
            carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
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
      for (i in 1:length(long)) 
        obs[, i] <<- inout(cbind(dist[, i], dif[, i]), cbind(polyX, polyY), bound = TRUE)
      
      graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
                graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
                alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
            label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
            carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
            axis = axes, legmap = legmap, legends = legends)
    }
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function()  {
    if (length(carte) != 0) {
      nocart <<- !nocart
      carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
            label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
            carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
            axis = axes, legmap = legmap, legends = legends)
      
    } else {
      tkmessageBox(message = "Spatial contours have not been given", icon = "warning", type = "ok")    
    }
  }
  
  ####################################################
  # Pour le alpha 
  ####################################################
  refresh.code <- function(...) {
    alpha <<- slider1(names.slide = names.slide, no = 1)
    graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
              graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
              alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- matrix(FALSE, nrow = length(long), ncol = length(long))
    
    graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
              graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
              alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)
    
    carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
          label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
          carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
          axis = axes, legmap = legmap, legends = legends)
  }
  
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function() {
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = baseenv())
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
    graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = dev.list()[length(dev.list())], 
              graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
              alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, num = dev.list()[length(dev.list())], lablong = lablong, lablat = lablat, 
          label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
          carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
          axis = axes, legmap = legmap, legends = legends)
    dev.off()
    
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = baseenv())
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
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
            label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
            carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
            axis = axes, legmap = legmap, legends = legends)
    } else {
      tkmessageBox(message = "Criteria has not been given", icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function()   {
    res2 <- choix.bubble(buble,listvar,listnomvar,legends, num_graph, num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
          label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
          carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
          axis = axes, legmap = legmap, legends = legends)
  }
  
  ####################################################
  # Dessin du variogramme 
  ####################################################
  
  vari <- function() {
    
    opt1 <<- ifelse(opt1 == 1, 2, 1)
    
    graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
              graph = "Variocloud", labvar = labvar, symbol = pch, couleurs = col, quantiles = quantiles, 
              alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)   
  }
  
  vari2 <- function() {
    opt2 <<- ifelse(opt2 == 1, 2, 1)
    
    graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
              graph = "Variocloud", labvar = labvar, symbol = pch, couleurs = col, quantiles = quantiles, 
              alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)   
  }
  
  ####################################################
  # Choisir angle
  ####################################################
  
  choixangle <- function() {
    
    directionnel <- !directionnel
    
    if (directionnel) {
      SGfunc()
      tt1 <- tktoplevel()
      Name <- tclVar("0.5")
      entry.Name <-tkentry(tt1, width = "3", textvariable = Name)
      tkgrid(tklabel(tt1, text = "Please enter a decimal x between 0 and 1 (angle=x.Pi)"), entry.Name)
      
      OnOK <- function() {
        angle <<- tclvalue(Name)
        tkdestroy(tt1)
        
        if (is.na(as.numeric(angle)) || (as.numeric(angle) > 1) || (as.numeric(angle) < 0)) {
          tkmessageBox(message = "Sorry, but you have to choose a decimal number between 0 and 1 (exemple : 0.5)",
                       icon = "warning", type = "ok")
        } else {
        msg <- paste("You choose", angle, "pi")
        tkmessageBox(message = msg)
        
        dist <<- sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
        dif <<-  (v1 - v2)^2
        dif2 <<-  (abs(v1 - v2))^(1/2)
        
        if(as.numeric(angle) < 1/10) { 
          dist[which((theta > as.numeric(angle)*pi + pi/10) & (theta < 9*pi/10 + as.numeric(angle)*pi))] <<- -10 
          dif[which((theta > as.numeric(angle)*pi + pi/10) & (theta < 9*pi/10 + as.numeric(angle)*pi))] <<- -10 
          dif2[which((theta > as.numeric(angle)*pi + pi/10) & (theta < 9*pi/10 + as.numeric(angle)*pi))] <<- -10
        } else if (as.numeric(angle)>9/10) {
          dist[which((theta > -2*pi + (as.numeric(angle)*pi + 11*pi/10)) & (theta < as.numeric(angle)*pi - pi/10))] <<- -10
          dif[which((theta > -2*pi + (as.numeric(angle)*pi + 11*pi/10)) & (theta < as.numeric(angle)*pi - pi/10))] <<- -10 
          dif2[which((theta > -2*pi + (as.numeric(angle)*pi + 11*pi/10)) & (theta < as.numeric(angle)*pi - pi/10))] <<- -10
        } else {
          dist[which((theta > as.numeric(angle)*pi + pi/10))] <<- -10
          dist[which((theta < as.numeric(angle)*pi - pi/10))] <<- -10
          dif[which((theta > as.numeric(angle)*pi + pi/10))] <<- -10
          dif[which((theta < as.numeric(angle)*pi - pi/10))] <<- -10
          dif2[which((theta > as.numeric(angle)*pi + pi/10))] <<- -10
          dif2[which((theta < as.numeric(angle)*pi - pi/10))] <<- -10
        }
        
        graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
                  graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
                  alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)  
        }
      }
      
      OK.but <-tkbutton(tt1, text = "   OK   ", command = OnOK)
      tkgrid(OK.but)
      tkfocus(tt1)
    } else {
      SGfunc()
      
      dist <<- sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
      dif <<-  (v1 - v2)^2/2
      dif2 <<-  (abs(v1 - v2))^(1/2)
      
      graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
                graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
                alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)   
    }
    
  }
  
  ####################################################
  # Repreentation Graphique
  ####################################################

  graphique(var1 = dist, var2 = dif, var3 = dif2, obs = obs, opt1 = opt1, opt2 = opt2, num = num_graph, 
                graph = "Variocloud", labvar = labvar, symbol = pch, couleurs=col, quantiles = quantiles, 
                alpha1 = alpha, bin = bin, xlim = xlim, ylim = ylim)  
      
  carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, lablat = lablat, 
            label = label, buble = buble, criteria = criteria, nointer = nointer, cbuble = z,
            carte = carte, nocart = nocart, cex.lab = cex.lab, method = "Variocloud",
            axis = axes, legmap = legmap, legends = legends)

  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "variocloudmap")
    
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
                   foreground = "darkred", background = "white"), side = "left", fill="x",expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = fbubble)
    tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    frame2c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2c, text = "Classic empirical semi-variogram ", font = "Times 11",
                   foreground = "darkred", background = "white"),
           tklabel(frame2c, text = "Robust empirical semi-variogram  ", font = "Times 11", 
                   foreground = "darkred", background = "white"), side = "left", fill="x",expand = "TRUE")
    tkpack(frame2c, expand = "TRUE", fill = "x")
    
    frame2d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    vari.but <- tkbutton(frame2d, text="On/Off", command=vari)
    vari2.but <- tkbutton(frame2d, text="On/Off", command=vari2)
    tkpack(vari.but, vari2.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2d, expand = "TRUE", fill = "x")
    
    frame2e <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2e, text = "Directional Semi-Variogram Cloud", font = "Times 11",
                   foreground = "darkred", background = "white"),side = "left", fill="x", expand = "TRUE")
    tkpack(frame2e, expand = "TRUE", fill = "x")
    
    frame2f <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    angle.but <- tkbutton(frame2f, text = "On/Off", command = choixangle)
    tkpack(angle.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2f, expand = "TRUE", fill = "x")
    
    
    if (quantiles) {
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
  
}

