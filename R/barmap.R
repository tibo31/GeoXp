barmap <- function(sp.obj, name.var, type = c("count", "percent"), names.arg = "", 
                   names.attr = names(sp.obj), criteria = NULL, carte = NULL, identify = FALSE,
                   cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "",
                   axes = FALSE, lablong = "", lablat = "") {
  
  envir <- as.environment(1)
  # Verification of the Spatial Object sp.obj
  class.obj <- class(sp.obj)[1]
  spdf <- (class.obj == "SpatialPolygonsDataFrame")
  
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
  
  listvar <- sp.obj@data
  listnomvar <- names.attr
  
  # for colors in map and new grahics
  if (length(col) == 1)
    col2 <- "blue"
  else
    col2 <- col
  
  col3 <- "lightblue3"
  
  # for identifyng the selected sites
  if(identify)
    label <- row.names(listvar)
  else
    label <- ""
  
  #initialisation
  nointer <- FALSE
  nocart <- FALSE
  buble <- FALSE
  z <- NULL
  legmap <- NULL
  labvar <- c(xlab, ylab)
  legends <- list(FALSE, FALSE, "", "")
  if (names.arg[1] == "")
    names.arg <- levels(as.factor(var))
  obs <- vector(mode = "logical", length = length(long))
  fin <- tclVar(FALSE)
  
  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix <- ""
  listgraph <- c("Histogram", "Barplot", "Scatterplot")
  
  #transformation data.frame en matrice
  if ((length(listvar) > 0) &&
      (dim(as.matrix(listvar))[2] == 1))
    listvar <- as.matrix(listvar)
  
  # Windows device
  if(length(dev.list()) == 0 & options("device") == "RStudioGD")
    dev.new()
  # if(!(2%in%dev.list())) 
  dev.new(noRStudioGD = FALSE)
  num_graph <- dev.list()[length(dev.list())]
  # if(!(3%in%dev.list())) 
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]
  # number of devices
  num_supp <- NA
  
  
  
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
      if (spdf & nrow(sp.obj) > 75 & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
      }
      
      #selection des points
      loc <- locator(1)
      if (is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
              cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
              lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
              method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
              legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
        next
      }
      
      if (!spdf | nrow(sp.obj) > 75) {
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, Xpoly = loc[1],
                          Ypoly = loc[2], method = "point")
      } else {
        if (gContains(sp.obj, SpatialPoints(cbind(loc$x, loc$y), proj4string = CRS(proj4string(sp.obj))))) {
          for (i in 1:nrow(sp.obj)) {
            if (gContains(sp.obj[i, ], SpatialPoints(cbind(loc$x, loc$y), proj4string =
                                                     CRS(proj4string(sp.obj))))) {
              obs[i] <<- !obs[i]
              break
            }
          }
        }
      }
      
      # graphiques
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
    
    if (spdf)
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
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
      
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
      
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
    
    carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
    assign("GeoXp.open", FALSE, envir = baseenv())
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
    carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = dev.list()[length(dev.list())],
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
    assign("GeoXp.open", FALSE, envir = baseenv())
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
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
    
    carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
          cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
          lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
          method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
          legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
  }
  
  ####################################################
  # Representation des graphiques
  ####################################################
  
  # Is there a Tk window already open ?
  if (interactive()) {
    if (!exists("GeoXp.open", envir = baseenv()) ||
        length(ls(envir = .TkRoot$env, all.names = TRUE)) == 2) {
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
            cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
            lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
            method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
            legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
      
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                symbol = pch, labmod = names.arg, couleurs = col)
      
      assign("GeoXp.open", TRUE, envir = baseenv())
    } else {
      if (get("GeoXp.open", envir = baseenv())) {
        stop("Warning : a GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")
      } else {
        carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
              cbuble = z, criteria = criteria, nointer = nointer, obs = obs, lablong = lablong, 
              lablat = lablat, label = label, symbol = pch, carte = carte, nocart = nocart,
              method = "Cluster", classe = var, couleurs = col2, legmap = legmap,
              legends = legends, labmod = names.arg, axis = axes, cex.lab = cex.lab)
        
        graphique(var1 = var, obs = obs, num = num_graph, graph = "Barplot", bin = type, labvar = labvar,
                  symbol = pch, labmod = names.arg, couleurs = col)
        assign("GeoXp.open", TRUE, envir = baseenv())
      }
    }
  }
  
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
      loc$name <- names(sp.obj[, name.var])
      legends <<- list(legends[[1]], TRUE, legends[[3]], loc)
      
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
      
      carte(long = long, lat = lat, buble = buble, sp.obj = sp.obj, num = num_carte,
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
