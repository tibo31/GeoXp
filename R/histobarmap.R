histobarmap <- function(sf.obj, names.var, nbcol = 10, type = "count",
                          names.arg = "", criteria = NULL, carte = NULL, identify = NULL,
                          cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = c("barplot", "histogram"), ylab = rep("count", 2),
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
  
  if (is.numeric(names.var)) {
    if (all(names.var <= ncol(listvar))) 
      names.var <- listnomvar[names.var]
    else
      stop("Dimension of names.var is not good")
  }
  
  if(!(all(names.var %in% names(sf.obj))))
    stop("At least one component of names.var is not included in the data.frame of sf.obj")
  
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
  num_graph1 <- dev.list()[length(dev.list())]
  # for map
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]
  # number of devices
  num_supp <- NA
  
  #####################################################
  ##### Arguments proper to each function 
  obs <- vector(mode = "logical", length = length(long))
  var1 <- sf.obj[[names.var[1]]]
  var2 <- sf.obj[[names.var[2]]]
  
  if(!(is.integer(var2) || is.double(var2))) 
    stop("Second element of names.var should be a numeric variable")
  
  # for colors in map and new grahics
  if(length(col) == 1)
    col2 <- "blue"
  else 
    col2 <- col
  
  col3 <- "lightblue3"
  labvar1 <- c(xlab[1], ylab[1])
  labvar2 <- c(xlab[2], ylab[2])
  if (names.arg[1] == "")
    names.arg <- levels(as.factor(var1))
  
  # for the 2nd density graph 
  dev.new(noRStudioGD = FALSE)
  num_graph2 <- dev.list()[length(dev.list())]
  
  ####################################################
  # selection d'un point
  ####################################################
  
  pointfunc <- function() {
    quit <- FALSE
    
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
        carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
              cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
              cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
              method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
              axis = axes, lablong = "", lablat = "")
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
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
                symbol = pch, labmod = names.arg, couleurs = col, bin = type)
      
      graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
                bin = type, labvar = labvar2, couleurs = col[1])
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2)) {
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                  symbol = pch, labvar = c(varChoice1, varChoice2))
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
      dev.set(num_carte)
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
      
      obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                        Xpoly = polyX, Ypoly = polyY, method = "poly")
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
      
      graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
                symbol = pch, labmod = names.arg, couleurs = col, bin = type)
      
      graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
                bin = type, labvar = labvar2, couleurs = col[1])
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                  symbol = pch, labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  # selection d'une barre sur le diagramme
  ####################################################
  
  bar1func <- function() {
    SGfunc()
    quit <- FALSE
    
    dev.set(num_graph1)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_graph1)
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
                  symbol = pch, labmod = names.arg, couleurs = col, bin = type)
        next
      }
      
      obs <<- selectstat(var1 = var1, obs = obs, Xpoly = loc[1], Ypoly = loc[2],
                         method = "Barplot")    
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
      
      graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
                symbol = pch, labmod = names.arg, couleurs = col, bin = type)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
                bin = type, labvar = labvar2, couleurs = col[1])
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                  symbol = pch, labvar = c(varChoice1, varChoice2))       
    }
  }
  
  ####################################################
  # selection d'une barre sur l'histogramme
  ####################################################
  
  bar2func <- function() {
    SGfunc()
    quit <- FALSE
    
    dev.set(num_graph2)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_graph2)
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
                  bin = type, labvar = labvar2, couleurs = col[1])
        next
      }
      
      obs <<- selectstat(var1 = var2, obs = obs, Xpoly = loc[1], 
                         Ypoly = loc[2], method = "Histogram", nbcol = nbcol)   
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
      
      graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
                symbol = pch, labmod = names.arg, couleurs = col, bin = type)
      
      graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
                bin = type, labvar = labvar2, couleurs = col[1])
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                  symbol = pch, labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  #  d'un autre graphique
  ####################################################
  
  graphfunc <- function() {
    if ((length(listvar) != 0) && (length(listnomvar) != 0)) {
      choix <<- selectgraph(listnomvar, listgraph)
      varChoice1 <<- choix$varChoice1
      varChoice2 <<- choix$varChoice2
      graphChoice <<- choix$graphChoice
      
      if ((graphChoice != "") && (varChoice1 != "")) {
        if(is.na(num_supp)) {
          dev.new(noRStudioGD = FALSE)
          num_supp <<- dev.list()[length(dev.list())]
        }
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                  symbol = pch, labvar = c(varChoice1, varChoice2))
      }
    } else {
      tkmessageBox(message = "Listvar and listnomvar must have been given in input",
                   icon = "warning", type = "ok")
    }
  }
  
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  
  cartfunc <- function() {  
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
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
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
          cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
          method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
          axis = axes, lablong = "", lablat = "")
    
    graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
              symbol = pch, labmod = names.arg, couleurs = col, bin = type)
    
    graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
              bin = type, labvar = labvar2, couleurs = col[1])
    
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                var2 = listvar[, which(listnomvar == varChoice2)],
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                symbol = pch, labvar = c(varChoice1, varChoice2))
  }
  
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function() {
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    dev.off(num_graph1)
    dev.off(num_graph2)
    dev.off(num_carte)
    if (!is.na(num_supp))
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
    graphique(var1 = var1, obs = obs, num = dev.list()[length(dev.list())], graph = "Barplot", labvar = labvar1, 
              symbol = pch, labmod = names.arg, couleurs = col, bin = type)
    dev.off()
    
    fig_save2 <- paste("fig_GeoXp", "_", k, ".pdf", sep = "")
    
    pdf(fig_save2)
    graphique(var1 = var2, obs = obs, num = dev.list()[length(dev.list())], graph = "Histogram", nbcol = nbcol, 
              bin = type, labvar = labvar2, couleurs = col[1])
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = dev.list()[length(dev.list())],
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
          cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
          method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
          axis = axes, lablong = "", lablat = "")
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
    cat("Figure 1 has been saved in", fig_save, "\n")
    cat("Figure 2 has been saved in", fig_save2, "\n")
    if(!is.na(num_supp))
      cat("Supplemental figure has been saved in", fig_supp, "\n")
    
    assign("last.select", which(obs), envir = envir)
    
    dev.off(num_carte)
    dev.off(num_graph1)
    dev.off(num_graph2)
    if(!is.na(num_supp))
      dev.off(num_supp)
  }
  
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer <- function() {

      if (length(criteria) != 0) {
      nointer <<- !nointer

      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "") 
    } else {
      tkmessageBox(message = "Criteria has not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function() {
    
    res2 <- choix.bubble(buble, listvar, listnomvar, legends, num_carte = num_carte)
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
          cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
          method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
          axis = axes, lablong = "", lablat = "")
    
  }
  
  ####################################################
  # Representation graphique
  ####################################################

  carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
      
  graphique(var1 = var1, obs = obs, num = num_graph1, graph = "Barplot", labvar = labvar1, 
                symbol = pch, labmod = names.arg, couleurs = col, bin = type)
      
  graphique(var1 = var2, obs = obs, num = num_graph2, graph = "Histogram", nbcol = nbcol, 
                bin = type, labvar = labvar2, couleurs = col[1])
   
  
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
      loc$name <- names(listvar[, names.var[1]])
      legends <<- list(legends[[1]], TRUE, legends[[3]], loc)
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")
      
    }
    
    OnOK2 <- function() { 
      legends <<- list(legends[[1]], FALSE, legends[[3]], "")
      tkdestroy(tt1)	
      
      carte(long = long, lat = lat, buble = buble, sf.obj = sf.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, label = label,
            cex.lab = cex.lab, symbol = pch, carte = carte, nocart = nocart, couleurs = col2,
            method = "Cluster", classe = var1, legmap = legmap, legends = legends, labmod = names.arg,
            axis = axes, lablong = "", lablat = "")      
    }
    
    
    if (length(col) == length(levels(as.factor(var1))) || length(pch) == length(levels(as.factor(var1)))) {
      tt1 <- tktoplevel()
      
      labelText12 <- tclVar("Do you want a legend for factors")
      label12 <- tklabel(tt1,justify = "center", wraplength = "3i", text = tclvalue(labelText12))
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
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "histobarmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    
    point.but <- tkbutton(frame1a, text="Selection by point", command = pointfunc)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = polyfunc)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE", fill = "x")
    
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Select bar(s)", font = "Times 12",
                   foreground = "darkred", background = "white"))
    barre1.but <- tkbutton(frame1c, text = "On barplot ", command = bar1func)
    barre2.but <- tkbutton(frame1c, text = "On histogram ", command = bar2func)
    tkpack(barre1.but,barre2.but, side = "left", expand = "TRUE", fill = "x")
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
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    tkpack(nocou1.but,noint1.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    frame2c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2c, text = "Bubbles", font = "Times 11",
                   foreground = "darkred", background = "white"),
           tklabel(frame2c, text = "Additional graph", font = "Times 11", 
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2c, expand = "TRUE", fill = "x")
    
    frame2d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    bubble.but <- tkbutton(frame2d, text = "On/Off", command = fbubble)
    autre.but <- tkbutton(frame2d, text = "     OK     " , command = graphfunc)
    tkpack(bubble.but,autre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2d, expand = "TRUE", fill = "x")
    
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

