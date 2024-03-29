pcamap <- function(sf.obj, names.var, direct = c(1, 2), weight = rep(1/nrow(sf.obj), length = nrow(sf.obj)), 
                   metric = diag(length(names.var)), center = NULL, reduce = TRUE, qualproj = FALSE, 
                   criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16, 
                   col = "lightblue3", xlab = paste(direct[1]), ylab = paste(direct[2]), axes = FALSE, 
                   lablong = "", lablat = "") {

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
  # for graphic 1
  dev.new(noRStudioGD = FALSE)
  num_graph1 <- dev.list()[length(dev.list())]
  # for map
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]
  # number of devices
  num_supp <- NA
  
  #####################################################
  ##### Arguments proper to each function 
  
  if (is.numeric(names.var)) {
    if (all(names.var <= ncol(listvar))) 
      names.var <- listnomvar[names.var]
    else
      stop("Dimension of names.var is not good")
  }
  
  if(!(all(names.var %in% names(sf.obj))))
    stop("names.var is not included in the sf object")
  
  dataset <- listvar[, names.var]
  name.dataset <- names(dataset)
  noms <- name.dataset
  
  dataset <- as.matrix(dataset)
  obs <- vector(mode = "logical", length = length(long))
  
  # calcul de l'ACP et recuperation des resultats
  p <- genpca(dataset, w = weight, m = metric, center = center, reduc = reduce)
  casecoord <- p$casecoord
  varcoord <- p$varcoord
  
  if ((max(varcoord[, 1]) + min(varcoord[, 1])) < 0) {
    varcoord[, 1] <- -varcoord[, 1]
    casecoord[, 1] <- -casecoord[, 1]
  }
  
  casecoord1 <- casecoord[, direct[1]]
  casecoord2 <- casecoord[, direct[2]]
  
  inertia <- p$inertia
  
  inerpart <- inertia / sum(inertia)
  inertpartperc <- inerpart * 100
  
  casequal <- casecoord[, direct[1]] ^ 2 + casecoord[, direct[2]] ^ 2
  den <- colSums(t(casecoord) * t(casecoord))
  casequal <- sqrt(casequal / den)
  casequalperc <- casequal * 100
  
  varqual <- varcoord[, direct[1]] ^ 2 + varcoord[, direct[2]] ^ 2
  denv <- colSums(t(varcoord) * t(varcoord))
  varqual <- sqrt(varqual / denv)
  varqualperc <- varqual * 100
  
  # for graphic 2
  dev.new(noRStudioGD = FALSE)
  num_graph2 <- dev.list()[length(dev.list())]
  
  # for colors in map and new grahics
  if (length(col) == 1)
    col2 <- "blue"
  else
    col2 <- col
  
  col3 <- "lightblue3"
  method <- ""
  pch2 <- pch[1]
  labmod <- ""
  maptest <- FALSE
  
  ####################################################
  # selection d'un point
  ####################################################
  
  pointfunc <- function() {
    quit <- FALSE
    if (maptest) { 
      dev.set(num_carte)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
      }
    
    } else {
      dev.set(num_graph1)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
    }
    
    while (!quit) {
      if (maptest) {
        dev.set(num_carte)
        loc <- locator(1)
        if (is.null(loc)) {
          quit <- TRUE 
          carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
                criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
                carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
                labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
                classe = listvar[, which(listnomvar == varChoice1)]) 
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

      } else {
        dev.set(num_graph1)
        loc <- locator(1)
        if (is.null(loc)) {
          quit<-TRUE
          graphique(var1 = casecoord1, var2 = casecoord2, obs = obs, num = num_graph1, graph = "Acp1",
                    symbol = pch, labmod = casequalperc, direct = direct, inertie = inertpartperc, label = qualproj, cex.lab = cex.lab,
                    labvar = labvar, couleurs = col)   
          next
        }
        
        obs <<- selectmap(var1 = casecoord1, var2 = casecoord2, obs = obs,
                          Xpoly = loc[1], Ypoly = loc[2], method = "point")
      }
      
      # graphiques
      
      graphique(var1 = casecoord1, var2 = casecoord2, obs = obs, num = num_graph1, graph = "Acp1",
                symbol = pch, labmod = casequalperc, direct = direct, inertie = inertpartperc, label = qualproj, cex.lab = cex.lab,
                labvar = labvar, couleurs = col)   
      
      # Remarque : s'il y a tous ces If, c'est pour prevoir de rajoutter un barplot avec options de couleurs
      # sur la carte 
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
            labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
      if (maptest) {
        dev.set(num_carte)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
        title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
        
        points(long, lat, pch = 16,col = "royalblue")
      } else { 
        dev.set(num_graph1)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
        title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
      }
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch, labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  # selection d'un point sur la carte
  ####################################################
  
  pt1func <- function() {
    maptest <<- TRUE
    pointfunc()
  }
  
  ####################################################
  # selection d'un point sur le graphique
  ####################################################
  
  pt2func <- function() {
    maptest <<- FALSE
    pointfunc()
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunc<-function() {
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    if (maptest) {
      dev.set(num_carte)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      points(long, lat, pch = 16, col = "royalblue") 
    
    } else { 
      dev.set(num_graph1)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
    }
    
    while (!quit) {
      if (maptest) {
        dev.set(num_carte)
        loc<-locator(1)
        if (is.null(loc)) {
          quit<-TRUE
          next
        }
      } else {
        dev.set(num_graph1)
        loc <- locator(1)
        if (is.null(loc)) {
          quit<-TRUE
          next
        }   
      }
      polyX <- c(polyX, loc[1])
      polyY <- c(polyY, loc[2])
      lines(polyX, polyY)
    }
    polyX <- c(polyX, polyX[1])
    polyY <- c(polyY, polyY[1])
    lines(polyX, polyY)
    
    if (length(polyX) > 0) {
      lines(polyX, polyY)
      
      if (maptest) {
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                          Xpoly = polyX, Ypoly = polyY, method = "poly")
      } else {
        obs <<- selectmap(var1 = casecoord[, direct[1]], var2 = casecoord[, direct[2]], 
                          obs = obs, Xpoly = polyX, Ypoly = polyY, method = "poly")
      }
      
      graphique(var1=casecoord[, direct[1]], var2=casecoord[, direct[2]], obs = obs, num = num_graph1, graph = "Acp1",
                symbol = pch, labmod = casequalperc, direct = direct, inertie = inertpartperc, label = qualproj, cex.lab = cex.lab,
                labvar = labvar, couleurs = col)   
      
      
      # Remarque : s'il y a tous ces If, c'est pour prevoir de rajoutter un barplot avec options de couleurs
      # sur la carte 
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
            labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch, labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  # selection d'un polygone sur la carte
  ####################################################
  
  poly1func <- function() {
    SGfunc()
    maptest <<- TRUE
    polyfunc()
  }
  
  ####################################################
  # selection d'un polygone sur le graphique
  ####################################################
  
  poly2func <- function() {
    SGfunc()
    maptest <<- FALSE
    polyfunc()
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() {  
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
            labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
    } else {
      tkmessageBox(message = "Spatial contours have not been given", icon = "warning", type = "ok")    
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
        if (((graphChoice == "Histogram") && (!is.numeric(listvar[, which(listnomvar == varChoice1)]))) 
            || ((graphChoice == "Scatterplot") && ((!is.numeric(listvar[, which(listnomvar == varChoice1)])) || 
                                                   (!is.numeric(listvar[, which(listnomvar == varChoice2)]))))) {
          tkmessageBox(message = "Variables choosed are not in a good format", icon = "warning", type = "ok")
        } else {
          res1 <- choix.couleur(graphChoice, listvar, listnomvar, varChoice1, legends, col, pch, spdf = F, 
                                num_graph1, num_carte)
          
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
          
          graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[, which(listnomvar == varChoice2)],
                    obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch, labvar = c(varChoice1, varChoice2))  
          
          carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
                criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
                carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
                labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
                classe = listvar[, which(listnomvar == varChoice1)]) 
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
    
    # graphiques
    graphique(var1 = casecoord1, var2 = casecoord2, obs = obs, num = num_graph1, graph = "Acp1",
              symbol = pch, labmod = casequalperc, direct = direct, inertie = inertpartperc, label = qualproj, cex.lab = cex.lab,
              labvar = labvar, couleurs = col)   
    
    # Remarque : s'il y a tous ces If, c'est pour pr?voir de rajoutter un barplot avec options de couleurs
    # sur la carte 
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
          criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
          carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
          labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
          classe = listvar[, which(listnomvar == varChoice1)]) 
    
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)],
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col3, symbol = pch, labvar = c(varChoice1, varChoice2))
    
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
    graphique(var1 = casecoord1, var2 = casecoord2, obs = obs, num = dev.list()[length(dev.list())], graph = "Acp1",
              symbol = pch, labmod = casequalperc, direct = direct, inertie = inertpartperc, label = qualproj, cex.lab = cex.lab,
              labvar = labvar, couleurs = col)   
    dev.off()
    
    fig_save2 <- paste("fig_GeoXp", "_", k, ".pdf", sep = "")  
    pdf(fig_save2)
    if (length(name.dataset) != dim(dataset)[2]) {
      graphique(var1 = varcoord[, direct[1]], var2 = varcoord[, direct[2]], obs = obs, num = dev.list()[length(dev.list())], graph = "Acp2", labvar = labvar, 
                legmap = noms, cex.lab = cex.lab, labmod = varqualperc, direct = direct, inertie = inertpartperc, couleurs = col)
    } else {
      graphique(var1 = varcoord[, direct[1]], var2 = varcoord[, direct[2]], obs = obs, num = dev.list()[length(dev.list())], graph = "Acp2", labvar = labvar, 
                legmap = name.dataset, cex.lab = cex.lab, labmod = varqualperc, direct = direct, inertie = inertpartperc, couleurs = col)
    }
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())], buble = buble, cbuble = z,
          criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
          carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
          labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
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
    assign("GeoXp.open", FALSE, envir = envir)
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure 1 has been saved in", fig_save, "\n")
    cat("Figure 2 has been saved in", fig_save2, "\n")
    if(!is.na(num_supp))
      cat("Supplemental figure has been saved in", fig_supp, "\n")
    
    assign("last.select", which(obs), envir = envir)
    
    dev.off(num_graph1)
    dev.off(num_graph2)
    dev.off(num_carte)
    if (!is.na(num_supp))
      dev.off(num_supp)
  }
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer <- function() {
    if (length(criteria) != 0) {
      nointer <<- !nointer
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
            labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
    } else {
      tkmessageBox(message = "Criteria has not been given", icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function() {
    res2 <- choix.bubble(buble,listvar,listnomvar,legends, num_carte = num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
          criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
          carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
          labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
          classe = listvar[, which(listnomvar == varChoice1)]) 
    
  }
  
  ####################################################
  # Representation graphique
  ####################################################
  
  graphique(var1 = casecoord1, var2 = casecoord2, obs = obs, num = num_graph1, graph = "Acp1",
                symbol = pch, labmod = casequalperc, direct = direct, inertie = inertpartperc, label = qualproj, cex.lab = cex.lab,
                labvar = labvar, couleurs = col)   
      
  carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2, 
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, 
            labmod = labmod, lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
  if (length(name.dataset) != dim(dataset)[2]) {
        graphique(var1 = varcoord[, direct[1]], var2 = varcoord[, direct[2]], obs = obs, num = num_graph2, graph = "Acp2", labvar = labvar, 
                 legmap = noms, cex.lab = cex.lab, labmod = varqualperc, direct = direct, inertie = inertpartperc, couleurs = col)
        } else {
          graphique(var1 = varcoord[, direct[1]], var2 = varcoord[, direct[2]], obs = obs, num = num_graph2, graph = "Acp2", labvar = labvar, 
                    legmap = name.dataset, cex.lab = cex.lab, labmod = varqualperc, direct = direct, inertie = inertpartperc, couleurs = col)
          }
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "pcamap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pt1func)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = poly1func)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the graph", font = "Times 12",
                   foreground = "darkred", background = "white"))
    p2.but <- tkbutton(frame1c, text = "Selection by point", command = pt2func)
    pol2.but <- tkbutton(frame1c, text = "Selection by polygon ", command = poly2func)
    tkpack(p2.but,pol2.but, side = "left", expand = "TRUE", fill = "x")
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

