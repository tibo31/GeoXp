densitymap <- function(sf.obj, name.var, kernel = "triweight",  
                       criteria = NULL, carte = NULL, identify = FALSE, cex.lab = 0.8,
                       pch = 16, col = "lightblue3", xlab = "", ylab = "", axes = FALSE,
                       lablong = "" , lablat = "") {
  
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
  obs <- vector(mode = "logical", length = length(long))
  graph <- "Densityplot2"
  
  # verify the type of the main variable
  if(!(is.integer(var) || is.double(var))) 
    stop("the variable name.var should be a numeric variable")
  
  # for the slider
  alpha1 <- 20
  names.slide <- "Alpha"
  polyX2 <- NULL
  
  # if add a graphic barplot
  labmod <- ""
  # if colors 
  col2 <- "blue"
  col3 <- col[1]
  pch2 <- pch[1]
  
  ####################################################
  # selection d'un point
  ####################################################
  
  pointfunc <- function() {  
    if (graph == "Densityplot2") {
    graph <<- "Densityplot1"
    
    graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
              graph = graph, labvar = labvar, symbol = pch, couleurs = col,
              kernel = kernel)
  }
    quit <- FALSE
    loc <- NULL
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
      points(long, lat, pch = 16, col = "royalblue")
    }
    
    while(!quit) {
      dev.set(num_carte)
      loc <- locator(1)
      
      if(is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
              buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
              label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
              legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
              lablat = lablat, symbol = pch2, couleurs = col2, method = method,
              classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
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
      if (length(var[obs]) >= 2)  {
        graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
                  graph = graph, labvar = labvar, symbol = pch, couleurs = col,
                  kernel = kernel)
      } else {
        dev.set(num_graph)
        title(sub = "You have to choose at least two sites to represent the sub-density", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
      }
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
            legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
            lablat = lablat, symbol = pch2, couleurs = col2, method = method,
            classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      
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
    if (graph == "Densityplot2") {
    graph <<- "Densityplot1"
    
    graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
              graph = graph, labvar = labvar, symbol = pch, couleurs = col,
              kernel = kernel)
  }
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
      if(is.null(loc)) {
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
      
      #graphique
      if (length(var[obs]) >= 2) {
        graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
                  graph = graph, labvar = labvar, symbol = pch, couleurs = col,
                  kernel = kernel)
      } else {
        dev.set(num_graph)
        title(sub = "You have to choose at least two sites to represent the sub-density", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
      }
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
            legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
            lablat = lablat, symbol = pch2, couleurs = col2, method = method,
            classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                  symbol = pch, labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  # selection d'une aire sous la courbe de densite
  ####################################################
  
  interfunc <- function() {
    
    if (graph == "Densityplot1") {
      SGfunc()
      graph <<- "Densityplot2"
    }
    
    n.inter <- length(polyX2)
    polyX <- NULL
    
    while (length(polyX) < 2) {
      dev.set(num_graph)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "Click two times to select an interval", cex.sub = 0.8, font.sub = 3, col.sub = "red")
      loc <- locator(1)
      polyX <- c(polyX, loc[1])
    }
    
    polyX2[[n.inter + 1]] <<- polyX
    
    obs <<- selectstat(var1 = var, obs = obs,
                       Xpoly = polyX[1], Ypoly = polyX[2],
                       method = "Densityplot")
    
    # graphiques
    graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
              graph = graph, Xpoly = polyX2, labvar = labvar, symbol = pch, couleurs = col,
              kernel = kernel)
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
          label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
          legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
          lablat = lablat, symbol = pch2, couleurs = col2, method = method,
          classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
    
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                var2 = listvar[, which(listnomvar == varChoice2)],
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                symbol = pch, labvar = c(varChoice1, varChoice2))
  }
  
  
  
  ####################################################
  # Choisir une valeur
  ####################################################
  
  choixvalue <- function() {
    if (graph == "Densityplot1") {
      SGfunc()
      graph <<- "Densityplot2"
    }
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    
    tt1 <- tktoplevel()
    Name <- tclVar("1st value")
    Name2 <- tclVar("2nd value")
    entry.Name <- tkentry(tt1, width = "8", textvariable = Name)
    entry.Name2 <- tkentry(tt1, width = "8", textvariable = Name2)
    tkgrid(tklabel(tt1, text = "Please enter values"), entry.Name, entry.Name2)
    
    OnOK <- function() {
      value1 <- tclvalue(Name)
      value2 <- tclvalue(Name2)
      n.inter <- length(polyX2)
      
      tkdestroy(tt1)
      if (is.na(as.numeric(value1)) || is.na(as.numeric(value2))) {
        tkmessageBox(message = "Sorry, but you have to choose decimal values",
                     icon = "warning", type = "ok")
      } else {
        polyX2[[n.inter + 1]] <<- c(as.numeric(value1), as.numeric(value2))
        obs <<- selectstat(var1 = var, obs = obs, Xpoly = as.numeric(value1), 
                           Ypoly = as.numeric(value2), method = "Densityplot")
        
        # graphiques
        graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
                  graph = graph, Xpoly = polyX2, labvar = labvar, symbol = pch, couleurs = col,
                  kernel = kernel)
        
        carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
              buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
              label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
              legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
              lablat = lablat, symbol = pch2, couleurs = col2, method = method,
              classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
        
        if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
          graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                    var2 = listvar[, which(listnomvar == varChoice2)],
                    obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, 
                    symbol = pch, labvar = c(varChoice1, varChoice2))
      }
    }
    
    OK.but <-tkbutton(tt1, text = "   OK   ", command = OnOK)
    tkgrid(OK.but)
    tkfocus(tt1)
  }
  
  ####################################################
  # modification du alpha pour l'estimateur de tous les points
  ####################################################
  
  refresh.code <- function(...) {
    
    alpha1 <<- slider1(names.slide = names.slide, no = 1)
    if (graph == "Densityplot1") { 
      if (length(var[obs]) > 1) {
        graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
                  graph = graph, labvar = labvar, symbol = pch, couleurs = col,
                  kernel = kernel)
        } else {
          dev.set(3)
        title(sub = "You have to choose at least two sites to represent the sub-density", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
      }
    } else {
      graphique(var1 = var, obs = obs, alpha1 = alpha1,  num = num_graph, 
                graph = "Densityplot2", Xpoly = polyX2, labvar = labvar, 
                symbol = pch, couleurs = col, kernel = kernel)
    }
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() {
    if (length(carte) != 0) {
      nocart <<- !nocart

      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
            legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
            lablat = lablat, symbol = pch2, couleurs = col2, method = method,
            classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
      } else {
      tkmessageBox(message = "Spatial contours have not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  
  ####################################################
  # choix d'un autre graphique
  ####################################################
  
  graphfunc <- function() {
    if ((length(listvar) != 0) && (length(listnomvar) != 0)) {
      choix <<- selectgraph(listnomvar,listgraph)
      varChoice1 <<- choix$varChoice1
      varChoice2 <<- choix$varChoice2
      graphChoice <<- choix$graphChoice
      
      if ((graphChoice != "") && (varChoice1 != "")) {
        if (((graphChoice == "Histogram")&&(!is.numeric(listvar[,which(listnomvar == varChoice1)]))) || 
            ((graphChoice == "Scatterplot") && ((!is.numeric(listvar[,which(listnomvar == varChoice1)])) || 
                                                (!is.numeric(listvar[, which(listnomvar == varChoice2)]))))) {
          tkmessageBox(message = "Variables choosed are not in a good format",
                       icon = "warning", type = "ok")
        } else {
          res1 <- choix.couleur(graphChoice, listvar, listnomvar, 
                                varChoice1, legends, col, pch, spdf = F,
                                num_graph, num_carte)
          method <<- res1$method
          col2 <<- res1$col2
          col3 <<- res1$col3
          pch2 <<- res1$pch2
          legends <<- res1$legends
          labmod <<- res1$labmod
          
          if(is.na(num_supp)) {
            dev.new(noRStudioGD = FALSE)
            num_supp <<- dev.list()[length(dev.list())]
          }
          
          graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                    var2 = listvar[, which(listnomvar == varChoice2)],
                    obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                    symbol = pch, labvar = c(varChoice1, varChoice2))
          
          if(is.na(num_supp)) {
            dev.new(noRStudioGD = FALSE)
            num_supp <<- dev.list()[length(dev.list())]
          }
          
          graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                    var2 = listvar[, which(listnomvar == varChoice2)],
                    obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                    symbol = pch, labvar = c(varChoice1, varChoice2))
          
          carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
                buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
                symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
                legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
                cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
          }
      }
    } else {
      tkmessageBox(message = "Variables (listvar) and their names (listnomvar) must have been given",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- vector(mode = "logical", length = length(long))
    graph <<- "Densityplot1"
    
    polyX2 <<- NULL
    
    graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, 
              graph = "Densityplot1", labvar = labvar, symbol = pch,
              couleurs = col, kernel = kernel)
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
          symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
          legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
          cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
    
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
    dev.off(num_graph)
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
    if (graph == "Densityplot1")
      graphique(var1 = var, obs = obs, alpha1 = alpha1, num = dev.list()[length(dev.list())], 
                graph = graph, labvar = labvar, symbol = pch, couleurs = col,
                kernel = kernel)
    else
    graphique(var1 = var, obs = obs, alpha1 = alpha1, num = dev.list()[length(dev.list())], 
              graph = graph, Xpoly = polyX2, labvar = labvar, symbol = pch, couleurs = col,
              kernel = kernel)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())],
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
          label = label, cex.lab = cex.lab, carte = carte, nocart = nocart, 
          legmap = legmap, legends = legends, axis = axes, lablong = lablong, 
          lablat = lablat, symbol = pch2, couleurs = col2, method = method,
          classe = listvar[, which(listnomvar == varChoice1)], labmod = labmod)
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
  
  fnointer <- function() {
    if (length(criteria) != 0) {
      
      nointer <<- !nointer
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
            symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
            legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
            cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
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
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
          symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
          legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
          cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
  }
  
  ####################################################
  # Representation graphique
  ####################################################

   graphique(var1 = var, obs = obs, alpha1 = alpha1, num = num_graph, graph = graph, 
                Xpoly = polyX2, labvar = labvar, symbol = pch, couleurs = col, kernel = kernel)
      
   carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
            symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
            legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
            cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14,
                                 weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "densitymap")
    
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
    tkpack(tklabel(frame1c, text = "Select an interval on the graphic", font = "Times 12",
                   foreground = "darkred", background = "white"))
    intervalle1.but <- tkbutton(frame1c, text = "by selecting on graph", command = interfunc)
    intervalle11.but <- tkbutton(frame1c, text = "by specifying bounds", command = choixvalue)
    tkpack(intervalle1.but,intervalle11.but, side = "left", expand = "TRUE", fill = "x")
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
                   foreground = "darkred", background = "white"), 
           side = "left", fill = "x", expand = "TRUE")
    tkpack(frame2c, expand = "TRUE", fill = "x")
    
    frame2d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    bubble.but <- tkbutton(frame2d, text = "On/Off", command = fbubble)
    autre.but <- tkbutton(frame2d, text = "     OK     " , command = graphfunc)
    tkpack(bubble.but,autre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2d, expand = "TRUE", fill = "x")
    
    frame2e <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    slider1(frame2e, refresh.code, names.slide, 3, 100, 1, alpha1)
    tkpack(frame2e, expand = "TRUE", fill = "x")
    
    frame3 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame3, text = "Exit", font = "Times 14",
                   foreground = "blue", background = "white"))
    
    quit.but <- tkbutton(frame3, text = "Save results", command = quitfunc2)
    quit.but2 <- tkbutton(frame3, text = "Exit without saving", command = quitfunc)
    
    tkpack(quit.but, quit.but2, side = "left", expand = "TRUE",
           fill = "x")
    
    tkpack(frame3, expand = "TRUE", fill = "x")
  }
  
  return(invisible())
}

