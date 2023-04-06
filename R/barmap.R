barmap <- function(sf.obj, name.var, type = c("count", "percent"), names.arg = "", 
                   criteria = NULL, carte = NULL, identify = NULL,
                   cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "",
                   axes = FALSE, lablong = "", lablat = "") {
  
 
  ###################################################
  ########## COMMON to ALL FUNCTIONS in GeoXp
  
  envir <- globalenv()
  # Verification of the Spatial Object sf.obj
  class.obj <- class(sf.obj)[1]
  
  if(class.obj != "sf") 
    stop("sf.obj may be a sf object")
  
  # verification on attributes
  listvar <- as.data.frame(st_drop_geometry(sf.obj))
  listnomvar <- colnames(listvar)
  
  if (is.numeric(name.var)) {
    if (name.var <= ncol(listvar)) 
      name.var <- listnomvar[name.var]
    else
      stop("Dimension of name.var is not good")
  }
  
  if(!(name.var %in% names(sf.obj)))
    stop("name.var is not included in the sf object")
  
  # we propose to refind the same arguments used in first version of GeoXp
  if (st_geometry_type(sf.obj, by_geometry = F) %in% c("POINT"))
    my_coords <- st_coordinates(st_geometry(sf.obj))
  else
    my_coords <- st_coordinates(st_point_on_surface(st_geometry(sf.obj)))
  long <- my_coords[, 1]
  lat <- my_coords[, 2]
  
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
  obs <- vector(mode = "logical", length = length(long))
  
  if (names.arg[1] == "")
    names.arg <- levels(as.factor(var))
  
  # for colors in map and new grahics
  if (length(col) == 1)
    col2 <- "blue"
  else
    col2 <- col
  
  col3 <- "lightblue3"
  ####################################################
  # selection d'un point
  ####################################################
  
  pointfunc <- function() {
    quit <- FALSE
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = 'red')
    title(sub = "To stop selection, click on the right button of the mouse or use ESC",
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_carte)
      
      if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
      }
      
      #selection des points
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
              cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
              lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
              method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
              legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
        next
      }
      
      if (nrow(sf.obj) > 100 | st_geometry_type(sf.obj, by_geometry = F) == "POINT")
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                          Xpoly = loc[1], Ypoly = loc[2], method = "point")
      else {
        my_points <- st_as_sf(data.frame(x = loc$x, y = loc$y), coords = c("x", "y"),
                              crs = st_crs(sf.obj))
        def <- as.vector(st_intersects(my_points, sf.obj, sparse = FALSE))
        obs[def] <<- !obs[def]
      }
      
      # graphiques
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC",
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2)) {
        graphique(var1 = listvar[, which(listnomvar == varChoice1)],
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch[1], 
                  labvar = c(varChoice1, varChoice2))
      }
    }
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunc <- function() {
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC",
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    points(long, lat, pch = 16, col = "royalblue")
    
    while (!quit) {
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
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
      
      #graphiques
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
      
      if ((graphChoice != "") &&
          (varChoice1 != "") && (length(dev.list()) > 2)) {
        graphique(var1 = listvar[, which(listnomvar == varChoice1)],
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch[1], 
                  labvar = c(varChoice1, varChoice2))
      }
    }
  }
  
  ####################################################
  # selection of a bar in the barplot
  ####################################################
  
  barfunc <- function() {
    SGfunc()
    quit <- FALSE
    
    while (!quit) {
      dev.set(num_graph)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC",
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                  symbol = pch, labmod = names.arg, couleurs = col)
        next
      }
      obs <<- selectstat(var1 = var, obs = obs, Xpoly = loc[1], Ypoly = loc[2], method = "Barplot")
      
      # graphiques
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC",
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      if ((graphChoice != "") &&
          (varChoice1 != "") && (length(dev.list()) > 2)) {
        graphique(var1 = listvar[, which(listnomvar == varChoice1)],
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch[1], 
                  labvar = c(varChoice1, varChoice2))
      }
    }
  }
  
  
  ####################################################
  # choix d'un autre graphique
  ####################################################
  
  graphfunc <- function() {
    if ((length(listvar) != 0) && (length(listnomvar) != 0)) {
      choix <<- selectgraph(listnomvar, listgraph)
      varChoice1 <<- choix$varChoice1
      varChoice2 <<- choix$varChoice2
      graphChoice <<- choix$graphChoice
      
      if ((length(graphChoice) != 0) && (length(varChoice1) != 0)) {

        if(is.na(num_supp)) {
          dev.new(noRStudioGD = FALSE)
          num_supp <<- dev.list()[length(dev.list())]
        }
        
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                  symbol = pch, labvar = c(varChoice1, varChoice2))
      } else {
        tkmessageBox(message = "You must choose a variable",
                     icon = "warning", type = "ok")
      }
    } else {
      tkmessageBox(message = "You must give variables (listvar) and their names (listnomvar)",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  
  cartfunc <- function()  {
    if (length(carte) != 0) {
      
      nocart <<- !nocart
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
    } else {
      tkmessageBox(message = "Spatial contours have not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- vector(mode = "logical", length = length(long))
    
    carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
          lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
          method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
          legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
    
    graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
              symbol = pch, labmod = names.arg, couleurs = col)
    
    if ((graphChoice != "") &&
        (varChoice1 != "") && (length(dev.list()) > 2)) {
      graphique(var1 = listvar[, which(listnomvar == varChoice1)],
                var2 = listvar[, which(listnomvar == varChoice2)],
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch[1], 
                labvar = c(varChoice1, varChoice2))
    }
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
    graphique(var1 = var, obs = obs, num = dev.list()[length(dev.list())], graph = "Barplot", bin = type, labvar = labvar,
              symbol = pch, labmod = names.arg, couleurs = col)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = dev.list()[length(dev.list())],
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
          lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
          method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
          legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
    dev.off()
    
    if(!is.na(num_supp)) {
      fig_supp <- "fig_supp_GeoXp.pdf"
      k <- 1
      while(file.exists(fig_supp)) {
        fig_supp <- paste("fig_supp_GeoXp", "_", k, ".pdf", sep = "")
        k <- k + 1
      }
      pdf(fig_supp)
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                var2 = listvar[, which(listnomvar == varChoice2)],
                obs = obs, num = dev.list()[length(dev.list())], graph = graphChoice, couleurs = col3, 
                symbol = pch, labvar = c(varChoice1, varChoice2))
      dev.off()
    }
    
    
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure has been saved in", fig_save, "\n")
    if(!is.na(num_supp))
      cat("Supplemental figure has been saved in", fig_supp, "\n")
    
    assign("last.select", which(obs), envir = envir)
    
    dev.off(num_carte)
    dev.off(num_graph)
    if(!is.na(num_supp))
      dev.off(num_supp)
  }
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer <- function()  {
    if (length(criteria) != 0) {
      nointer <<- !nointer
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
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
    
    carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
          lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
          method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
          legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
    
  }
  
  ####################################################
  # Representation des graphiques
  ####################################################
  
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
      
  
  ####################################################
  # creation de la boite de dialogue to create legens
  ####################################################
  
  if (interactive()) {
    OnOK <- function() {
      tkdestroy(tt1)
      msg <- paste("Click on the map to indicate the location of the upper left corner of the legend box")
      tkmessageBox(message = msg)
      
      dev.set(num_carte)
      loc <- locator(1)
      loc$name <- names(listvar[, name.var])
      legends <<- list(legends[[1]], TRUE, legends[[3]], loc)
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
    }
    
    OnOK2 <- function() {
      legends <<- list(legends[[1]], FALSE, legends[[3]], "")
      tkdestroy(tt1)
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
    }
    
    
    
    if (length(col) == length(levels(as.factor(var))) ||
        length(pch) == length(levels(as.factor(var)))) {
      
      tt1 <- tktoplevel()
      
      labelText12 <- tclVar("Do you want a legend for factors")
      label12 <-
        tklabel(tt1, justify = "center", wraplength = "3i",
                text = tclvalue(labelText12))
      tkconfigure(label12, textvariable = labelText12)
      tkgrid(label12, columnspan = 2)
      
      point.but <- tkbutton(tt1, text = "  Yes  ", command = OnOK)
      
      poly.but <- tkbutton(tt1, text = " No ", command = OnOK2)
      
      tkgrid(point.but, poly.but)
      tkgrid(tklabel(tt1, text = "    "))
      
      tkfocus(tt1)
    }
    
  }
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14,
                                 weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "barmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pointfunc)
    
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = polyfunc)
    
    tkpack(point.but, poly.but, side = "left", expand = "TRUE", fill = "x")
    
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the bar plot", font = "Times 12", 
                   foreground = "darkred", background = "white"))
    barre.but <- tkbutton(frame1c, text = "  Bar  ", command = barfunc)
    
    tkpack(barre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1c, expand = "TRUE", fill = "x")
    
    frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nettoy.but <- tkbutton(frame1b, text = "     Reset selection     " , command = SGfunc)
    
    tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1b, expand = "TRUE", fill = "x")
    
    frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2, text = "Options", font = "Times 14", 
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame2, text = "Spatial contours", font = "Times 11",
                   foreground = "darkred", background = "white"),
           tklabel(frame2, text = "Preselected sites", font = "Times 11",
                   foreground = "darkred", background = "white"),
           side = "left", fill = "x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    tkpack(nocou1.but, noint1.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    frame2c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2c, text = "Bubbles", font = "Times 11", foreground = "darkred",
                   background = "white"),
           tklabel(frame2c, text = "Additional graph", font = "Times 11", foreground = "darkred",
                   background = "white"),
           side = "left", fill = "x", expand = "TRUE")
    tkpack(frame2c, expand = "TRUE", fill = "x")
    
    frame2d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    bubble.but <- tkbutton(frame2d, text = "On/Off", command = fbubble)
    autre.but <- tkbutton(frame2d, text = "     OK     " , command = graphfunc)
    tkpack(bubble.but, autre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2d, expand = "TRUE", fill = "x")
    
    frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame3, text = "Exit", font = "Times 14", foreground = "blue",
                   background = "white"))
    
    quit.but <- tkbutton(frame3, text = "Save results", command = quitfunc2)
    
    quit.but2 <- tkbutton(frame3, text = "Exit without saving", command = quitfunc)
    
    tkpack(quit.but, quit.but2, side = "left", expand = "TRUE", fill = "x")
    
    tkpack(frame3, expand = "TRUE", fill = "x")
  }
  
  ####################################################
  # Fin
  ####################################################
  
  return(invisible())
}
