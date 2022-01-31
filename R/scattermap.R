scattermap <- function(sp.obj, names.var, lin.reg = TRUE, quantiles = TRUE,
                         names.attr = names(sp.obj), criteria = NULL, carte = NULL, 
                         identify = FALSE, cex.lab = 0.8, pch = 16, col = "lightblue3",
                         xlab = "", ylab = "", axes = FALSE, lablong = "", lablat = "") {
  
  # Verification of the Spatial Object sp.obj
  class.obj <- class(sp.obj)[1]
  spdf <- (class.obj=="SpatialPolygonsDataFrame")
  
  if(substr(class.obj, 1, 7) != "Spatial") 
    stop("sp.obj may be a Spatial object")
  
  if(substr(class.obj,nchar(class.obj) - 8, nchar(class.obj)) != "DataFrame") 
    stop("sp.obj should contain a data.frame")
  
  if(!is.numeric(names.var) & length(match(names.var, names(sp.obj))) != length(names.var)) 
    stop("At least one component of names.var is not included in the data.frame of sp.obj")
  
  if(length(names.attr) != length(names(sp.obj))) 
    stop("names.attr should be a vector of character with a length equal to the number of variable")
  
  # we propose to refind the same arguments used in first version of GeoXp
  long <- coordinates(sp.obj)[, 1]
  lat <- coordinates(sp.obj)[, 2]
  
  var1 <- sp.obj@data[, names.var[1]]
  var2 <- sp.obj@data[, names.var[2]]
  
  listvar <- sp.obj@data
  listnomvar <- names.attr
  
  # for identifyng the selected sites
  if(identify)
    label <- row.names(listvar)
  else
    label <- ""
  
  # initialisation
  nointer <- FALSE
  nocart <- FALSE
  buble <- FALSE
  maptt <- FALSE
  z <- NULL
  legmap <- NULL
  legends <- list(FALSE, FALSE, "", "")
  labvar <- c(xlab, ylab)
  obs <- vector(mode = "logical", length = length(long))
  
  var1 <- as.matrix(var1)
  var2 <- as.matrix(var2)
  lat <- as.matrix(lat)
  long <- as.matrix(long)
  
  graphChoice <- ""
  varChoice1 <- ""
  varChoice2 <- ""
  choix <- ""
  listgraph <- c("Histogram","Barplot","Scatterplot")
  
  labmod <- ""
  col2 <- "blue"
  col3 <- col[1]
  method <- ""
  pch2 <- pch[1]
  
  # transformation data.frame en matrix
  if((length(listvar) > 0) && (dim(as.matrix(listvar))[2] == 1)) 
    listvar <- as.matrix(listvar)
  
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
  # if(!(2%in%dev.list())) 
  dev.new(noRStudioGD = FALSE)
  num_graph <- dev.list()[length(dev.list())]
  # if(!(3%in%dev.list())) 
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]
  # number of devices
  num_supp <- NA
  
  borne1 <- 0.01
  borne2 <- 0.99
  alpha <- 0.5
  temp_data <- data.frame(
    var1 = sort(var1),
    var2 = var2[order(var1)])
  alpha1 <- qgam(var2 ~ s(var1, k = 20, bs = "ad"), data = temp_data, qu = alpha)
  
  names.slide <- "Value of alpha-quantile"
  
  ####################################################
  # selection d'un point
  ####################################################
  
  pointfunc<-function() {
    quit <- FALSE
    
    if (maptt) { 
      dev.set(num_carte)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      if (spdf & nrow(sp.obj) > 85 & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
        }
    } else { 
      dev.set(num_graph)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
    }
    
    while (!quit) {
      if (maptt) {
        dev.set(num_carte)
        if(spdf & nrow(sp.obj) > 75 & !buble) {
          points(long, lat, pch = 16, col = "royalblue")
        }
        loc <- locator(1)
        if(is.null(loc)) {
          quit <- TRUE
          carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
                buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
                label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
                legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
                lablat = lablat, cex.lab = cex.lab, method = method, 
                classe = listvar[, which(listnomvar == varChoice1)]) 
          next
        }           
        if(!spdf | nrow(sp.obj) > 75) {
          obs <<- selectmap(var1 = long, var2 = lat, obs = obs, Xpoly = loc[1], 
                            Ypoly = loc[2], method = "point")
          } else {
            if (gContains(sp.obj, SpatialPoints(cbind(loc$x, loc$y), proj4string = CRS(proj4string(sp.obj))))) {
              for (i in 1:nrow(sp.obj)) {
                if (gContains(sp.obj[i, ], SpatialPoints(cbind(loc$x, loc$y), proj4string = CRS(proj4string(sp.obj))))) {
                  obs[i] <<- !obs[i]
                  break
                }  
              }
            } 
          }
      } else {
        dev.set(num_graph)
        loc <- locator(1)
        if(is.null(loc)) {
          quit<-TRUE
          graphique(var1 = var1, var2 = var2, obs = obs, num = num_graph, graph = "Scatterplot", 
                    labvar = labvar, symbol = pch, couleurs = col, opt1 = lin.reg,
                    quantiles = quantiles, alpha1 = alpha1)
          next
        }
        obs <<- selectmap(var1 = var1, var2 = var2, obs = obs, 
                          Xpoly = loc[1], Ypoly = loc[2], method = "point")
      }
      
      # graphiques     
      graphique(var1 = var1, var2 = var2, obs = obs, num = num_graph, graph = "Scatterplot", 
                labvar = labvar, symbol = pch, couleurs = col, opt1 = lin.reg, 
                quantiles = quantiles, alpha1 = alpha1)
      
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
            legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
            lablat = lablat, cex.lab = cex.lab, method = method, 
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
      if(maptt) {
        dev.set(num_carte)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
        title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
      } else {
        dev.set(num_graph)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
        title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
      }
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2)) {
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                  symbol = pch, labvar = c(varChoice1, varChoice2))
      }
    }
  }
  
  
  ####################################################
  # selection d'un point sur la carte
  ####################################################
  
  pt1func <- function() {
    maptt <<- TRUE 
    pointfunc()
  }
  
  ####################################################
  # selection d'un point sur le graphique
  ####################################################
  
  pt2func <- function() {
    maptt <<- FALSE
    pointfunc()
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunc<-function() {
    
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    if(maptt) {
      dev.set(num_carte)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      if(spdf)
        points(long, lat, pch = 16, col = "royalblue")
    } else { 
      dev.set(num_graph)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
    }
    
    while(!quit) {
      if(maptt) {
        dev.set(num_carte)
        loc <- locator(1)
        if (is.null(loc)) {
          quit <- TRUE
          next
        }
      } else {
        dev.set(num_graph)
        loc <- locator(1)
        if (is.null(loc)) {
          quit <- TRUE
          next
        }   
      }
      polyX <- c(polyX, loc[1])
      polyY <- c(polyY, loc[2])
      lines(polyX, polyY)
    }
    polyX <- c(polyX, polyX[1])
    polyY <- c(polyY, polyY[1])
    
    if (length(polyX)>0) {
      lines(polyX, polyY)
      
      if (maptt) {
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                          Xpoly = polyX, Ypoly = polyY, method = "poly")
      } else {
        obs <<- selectmap(var1 = var1, var2 = var2, obs = obs, Xpoly = polyX, 
                          Ypoly = polyY, method = "poly")
      }
      
      # graphiques
      graphique(var1 = var1, var2 = var2, obs = obs, num = num_graph, graph = "Scatterplot", 
                labvar = labvar, symbol = pch, couleurs = col, opt1 = lin.reg, 
                quantiles = quantiles, alpha1 = alpha1)
      
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
            legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
            lablat = lablat, cex.lab = cex.lab, method = method, 
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2)) {
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                  var2 = listvar[, which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                  symbol = pch, labvar = c(varChoice1, varChoice2))
      }
    }
  }
  
  ####################################################
  # selection d'un polygone sur la carte
  ####################################################
  
  poly1func <- function() {
    maptt <<- TRUE
    polyfunc()
  }
  
  ####################################################
  # selection d'un polygone sur le graphique
  ####################################################
  
  poly2func <- function() {
    maptt <<- FALSE
    polyfunc()
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() {  
    
    if (length(carte) != 0) {
      
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
            legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
            lablat = lablat, cex.lab = cex.lab, method = method, 
            classe = listvar[, which(listnomvar == varChoice1)]) 
    } else {
      tkmessageBox(message="Spatial contours have not been given",
                   icon = "warning", type = "ok")    
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
      
      if ((graphChoice != "") && (varChoice1 != "")) {
        if (((graphChoice == "Histogram") && (!is.numeric(listvar[, which(listnomvar == varChoice1)]))) || 
            ((graphChoice == "Scatterplot")&&((!is.numeric(listvar[, which(listnomvar == varChoice1)])) || 
                                              (!is.numeric(listvar[, which(listnomvar == varChoice2)]))))) {
          tkmessageBox(message = "Variables choosed are not in a good format",
                       icon = "warning", type = "ok")
        } else {
          res1 <- choix.couleur(graphChoice, listvar, listnomvar, 
                                varChoice1, legends, col, pch, spdf = spdf,
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
          
          carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
                buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
                label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
                legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
                lablat = lablat, cex.lab = cex.lab, method = method, 
                classe = listvar[, which(listnomvar == varChoice1)]) 
        }
      }   
    } else {
      tkmessageBox(message="Variables (listvar) and their names (listnomvar) must have been given",
                   icon = "warning", type = "ok")
    }  
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc<-function() {
    obs <<- vector(mode = "logical", length = length(long))
    
    # graphiques
    graphique(var1 = var1, var2 = var2, obs = obs, num = num_graph, graph = "Scatterplot", 
              labvar = labvar, symbol = pch, couleurs = col, opt1 = lin.reg, 
              quantiles = quantiles, alpha1 = alpha1)
    
    carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
          label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
          legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
          lablat = lablat, cex.lab = cex.lab, method = method, 
          classe = listvar[, which(listnomvar == varChoice1)]) 
    
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2)) {
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], 
                var2 = listvar[, which(listnomvar == varChoice2)],
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                symbol = pch, labvar = c(varChoice1, varChoice2))
    } 
  }
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function() {
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = globalenv())
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
    graphique(var1 = var1, var2 = var2, obs = obs, num = dev.list()[length(dev.list())], 
              graph = "Scatterplot", labvar = labvar, symbol = pch, couleurs = col, 
              opt1 = lin.reg, quantiles = quantiles, alpha1 = alpha1)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = dev.list()[length(dev.list())],
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
          label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
          legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
          lablat = lablat, cex.lab = cex.lab, method = method, 
          classe = listvar[, which(listnomvar == varChoice1)]) 
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
    assign("GeoXp.open", FALSE, envir = globalenv())
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure has been saved in", fig_save, "\n")
    if(!is.na(num_supp))
      cat("Supplemental figure has been saved in", fig_supp, "\n")
    
    assign("last.select", which(obs), envir = globalenv())
    
    dev.off(num_carte)
    dev.off(num_graph)
    if(!is.na(num_supp))
      dev.off(num_supp)
  }
  
  ####################################################
  # modification du alpha de l'estimateur
  ####################################################
  
  refresh.code <- function(...) {
    alpha <<- slider1(names.slide = names.slide, no = 1)
    
    temp_data <- data.frame(
      var1 = sort(var1),
      var2 = var2[order(var1)])
    
    alpha1 <<- qgam(var2 ~ s(var1, k = 20, bs = "ad"), data = temp_data, qu = alpha)
    
    graphique(var1 = var1, var2 = var2, obs = obs, num = num_graph, graph = "Scatterplot", 
              labvar = labvar, symbol = pch, couleurs = col, opt1 = lin.reg, 
              quantiles = quantiles, alpha1 = alpha1)
   
   
    #fitmax  <- rqss(var2 ~ qss(var1, constraint= "CD", lambda = etendue), 
    #                tau = alpha1, data = temp_data, control = sfn.control(warn.mesg = FALSE))    
    # lines(xSeq$var1, pred$fit, col = "blue") 
  }
  
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer<-function() {
    if (length(criteria) != 0) {
      nointer <<- !nointer
      carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
            legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
            lablat = lablat, cex.lab = cex.lab, method = method, 
            classe = listvar[, which(listnomvar == varChoice1)]) 
    } else {
      tkmessageBox(message="Criteria has not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function() {
    res2 <- choix.bubble(buble, listvar, listnomvar, 
                         legends, num_graph, num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
          buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
          label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
          legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
          lablat = lablat, cex.lab = cex.lab, method = method, 
          classe = listvar[, which(listnomvar == varChoice1)]) 
  }
  
  ####################################################
  # Representation des graphiques
  ####################################################
  graphique(var1 = var1, var2 = var2, obs = obs, num = num_graph, graph = "Scatterplot", 
                labvar = labvar, symbol = pch, couleurs = col, opt1 = lin.reg, 
                quantiles = quantiles, alpha1 = alpha1)
      
  carte(long = long, lat = lat, obs = obs, sp.obj = sp.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer,
            label = label, symbol = pch2, couleurs = col2, carte = carte, nocart = nocart,
            legmap = legmap, legends = legends, axis = axes, labmod = labmod, lablong = lablong,
            lablat = lablat, cex.lab = cex.lab, method = method, 
            classe = listvar[, which(listnomvar == varChoice1)]) 
   
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  if(interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "scattermap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text="Selection by point", command = pt1func)
    poly.but <- tkbutton(frame1a, text="Selection by polygon ", command = poly1func)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the graphic", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but2 <- tkbutton(frame1c, text = "Selection by point", command = pt2func)
    poly.but2 <- tkbutton(frame1c, text = "Selection by polygon ", command = poly2func)
    tkpack(point.but2, poly.but2, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1c, expand = "TRUE", fill = "x")
    
    frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nettoy.but <- tkbutton(frame1b, text = "     Reset selection     ", command=SGfunc)
    tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1b, expand = "TRUE", fill = "x")
    
    frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2, text = "Options", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame2, text = "Spatial contours", font = "Times 11",
                   foreground = "darkred", background = "white"), 
           tklabel(frame2, text = "Preselected sites", font = "Times 11",
                   foreground = "darkred", background = "white"), 
           side = "left", fill = "x",expand = "TRUE")
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
           side = "left", fill = "x",expand = "TRUE")
    tkpack(frame2c, expand = "TRUE", fill = "x")
    
    frame2d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    bubble.but <- tkbutton(frame2d, text = "On/Off", command = fbubble)
    autre.but <- tkbutton(frame2d, text = "     OK     " , command = graphfunc)
    tkpack(bubble.but,autre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2d, expand = "TRUE", fill = "x")
    
    if(quantiles) {
      frame2e <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
      
      slider1(frame2e, refresh.code, names.slide = names.slide,
              borne1, borne2, (borne2 - borne1)/100, alpha)
      
      tkpack(frame2e, expand = "TRUE", fill = "x")
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

