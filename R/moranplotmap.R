moranplotmap <- function(sf.obj, name.var, listw.obj, flower = FALSE, locmoran = FALSE, names.arg = c("H.-H.", "L.-H.", "L.-L.", "H.-L."),
                         criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16, col = "lightblue3",
                         xlab = "", ylab = "", axes = FALSE, lablong = "", lablat = "") {

  
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
  
  # verify the type of the main variable
  if(!(is.integer(var) || is.double(var))) 
    stop("the variable name.var should be a numeric variable")
  
  # We create a spatial weight matrix by using a matrix object
  n <- nrow(sf.obj)
  W <- matrix(0, n, n)
  W.sn <- listw2sn(listw.obj)
  W[as.matrix(W.sn[, 1:2])] <- W.sn[, 3]
  
  # Is W normalized ?
  is.norm <- all(apply(W, 1, sum) == rep(1, n))
  
  #initialisation
  if (xlab == "")
    xlab <- name.var
  if (ylab == "")
    ylab <- paste("spatially lagged", name.var)
  
  obsq <- rep(0, length(long))
  buble2 <- FALSE
  maptest <- FALSE
  classe <- rep(1, length(long))
  
  legends2 <- list(FALSE, FALSE, "", "")
  
  z2 <- NULL
  legmap2 <- NULL
  num <- NULL
  labvar <- c(xlab, ylab)
  quad <- FALSE
  
  # Option sur le moran
  method <- ifelse(flower, "Neighbourplot1", "")
  
  choix.col <- FALSE
  
  graph <- "Moran"
  col2 <- rep(col[1], 4)
  col3 <- rep("blue", 4)
  pch2 <- rep(pch[1], 4)
  
  wvar <- W %*% var
  
  stdvar <- var / sd(var)
  uns <- rep(1, length(var))
  result <- nonormmoran(stdvar, uns, W)
  MORAN <- result$morani
  prob.I <- pnorm(result$istat)
  rvar <- qr(var)
  beta.I <- qr.coef(rvar, wvar)
  
  # calcul de la variable obsq (pour les quadrants)
  obsq[which((var > mean(var)) & (wvar >= mean(wvar)))] <- 1
  obsq[which((var <= mean(var)) & (wvar > mean(wvar)))] <- 2
  obsq[which((var < mean(var)) & (wvar <= mean(wvar)))] <- 3
  obsq[which((var >= mean(var)) & (wvar < mean(wvar)))] <- 4
  
  # i de moran local
  x.centre <- (var - mean(var))
  wx.centre <- (W %*% x.centre)
  ilocal <- (x.centre / var(x.centre)) * (wx.centre)  
  
  ilocal <- matrix(spdep::localmoran(var, listw = listw.obj)[, 4], ncol = 1)
  
  ####################################################
  # selection d'un point
  ####################################################
  
  pointfunc <- function() {
    quit <- FALSE
    quad <<- FALSE
    
    if (maptest) { 
      dev.set(num_carte)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
      }
    } else { 
      dev.set(num_graph)
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
          carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
                carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
                W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
                nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
                label = label, cex.lab = cex.lab, labmod = names.arg)   
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
        dev.set(num_graph)
        loc <- locator(1)
        if (is.null(loc)) {
          quit <- TRUE
          graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
                    graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
                    locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
                    legmap = legmap2, legends = legends2, bin = is.norm)
          next
        }
        obs <<- selectmap(var1 = var, var2 = wvar, obs = obs,
                          Xpoly = loc[1], Ypoly = loc[2], method = "point")
      }
      
      #graphiques
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
            W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
            nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
            label = label, cex.lab = cex.lab, labmod = names.arg)    
      
      graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
                graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
                locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
                legmap = legmap2, legends = legends2, bin = is.norm)
      
      if (maptest) {
        dev.set(num_carte)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
        title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
        if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
          points(long, lat, pch = 16, col = "royalblue")
        }
      } else {
        dev.set(num_graph)
        title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
        title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
              cex.sub = 0.8, font.sub = 3, col.sub = "red")
        }
      
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)], 
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col[1], symbol = pch[1], labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  # selection d'un point sur la carte
  ####################################################
  
  pt1func <- function() {
    method <<- ifelse(flower, "Neighbourplot1", "Quadrant")
    graph <<- "Moran" 
    maptest <<- TRUE  
    pointfunc()
  }
  
  ####################################################
  # selection d'un point sur le graphique
  ####################################################
  
  pt2func <- function() {    
    method <<- ifelse(flower, "Neighbourplot1", "Quadrant")
    graph <<- "Moran" 
    maptest <<- FALSE  
    pointfunc()
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunc <- function() {
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    quad <<- FALSE
    
    if (maptest) { 
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
    
    while (!quit) {
      if (maptest) {
        dev.set(num_carte)
        points(long, lat, pch = 16, col = "royalblue")
        loc <- locator(1)
        if (is.null(loc)) {
          quit<-TRUE
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
    
    if (length(polyX) > 0) {
      lines(polyX, polyY)
      
      if (maptest) {
        obs <<- selectmap(var1 = long, var2 = lat, obs = obs, Xpoly = polyX, 
                          Ypoly = polyY, method = "poly")
      } else {
        obs <<- selectmap(var1 = var, var2 = wvar, obs = obs, 
                          Xpoly = polyX, Ypoly = polyY, method = "poly")
      }
      
      #graphiques
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
            W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
            nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
            label = label, cex.lab = cex.lab, labmod = names.arg)         
      
      graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
                graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
                locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
                legmap = legmap2, legends = legends2, bin = is.norm)
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)], 
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col[1], symbol = pch[1], labvar = c(varChoice1, varChoice2))
    }
  }
  
  ####################################################
  # selection d'un polygone sur la carte
  ####################################################
  
  poly1func <- function() {
    method <<- ifelse(flower, "Neighbourplot1", "")
    graph <<- "Moran" 
    
    if (quad) {
      SGfunc()
      quad <<- FALSE
    }
    
    maptest <<- TRUE
    polyfunc()
  }
  
  ####################################################
  # selection d'un polygone sur le graphique
  ####################################################
  
  poly2func <- function() {
    
    method <<- ifelse(flower, "Neighbourplot1", "")
    graph <<- "Moran" 
    
    if (quad) {
      SGfunc()
      quad <<- FALSE
    }
    maptest <<- FALSE
    polyfunc()
  }
  
  ####################################################
  # selection d'un quadrant
  ####################################################
  
  quadfunc <- function() {
    if (!quad) {
      SGfunc()
      quad <<- TRUE
    }
    
    obs[which(obsq == num)] <<- !obs[which(obsq == num)]
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
          carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
          W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
          nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
          label = label, cex.lab = cex.lab, labmod = names.arg)   
    
    graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
              graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
              locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
              legmap = legmap2, legends = legends2, bin = is.norm)
    
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) > 2))
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)], 
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col[1], symbol = pch[1], labvar = c(varChoice1, varChoice2))
  }
  
  ####################################################
  # selection d'un des 4 quadrants 
  ####################################################
  
  quad1func <- function() {
    num <<- 1 
    method <<- ifelse(flower, "Neighbourplot1", "Quadrant") 
    quadfunc()
  }
  
  quad2func <- function() {
    num <<- 2  
    method <<- ifelse(flower, "Neighbourplot1", "Quadrant") 
    quadfunc()
  }
  
  quad3func <- function() {
    num <<- 3
    method <<- ifelse(flower, "Neighbourplot1", "Quadrant") 
    quadfunc()
  }
  
  quad4func <- function() {
    num <<- 4
    method <<- ifelse(flower, "Neighbourplot1", "Quadrant") 
    quadfunc()
  }
  
  ####################################################
  # Differentes couleurs selon le quadrant
  ####################################################
  
  colfunc <- function() {
    if (!choix.col) {
      choix.col <<- TRUE 
      method <<- "Quadrant" 
      res1 <- choix.couleur("Moran", col = col[1], pch = pch[1], legends = legends, spdf = T, 
                            num_graph = num_graph, num_carte = num_carte)     
      
    if (length(res1$col2) == 4) {
      col2 <<- res1$col2
      col3 <<- col2 
      } else {
        col2 <<- rep(col[1], 4) 
        col3 <<- rep("blue", 4)
      }
      
    if (length(res1$pch2) == 4)  
      pch2 <<- res1$pch2
      else 
        pch2 <<- rep(pch[1], 4)
    
      legends <<- res1$legends
    } else {
      choix.col <<- FALSE 
      col2 <<- rep(col[1], 4) 
      col3 <<- rep("blue", 4) 
      pch2 <<- rep(pch[1], 4) 
      legends <<- list(legends[[1]], FALSE, legends[[3]], "") 
      }
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
          carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
          W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
          nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
          label = label, cex.lab = cex.lab, labmod = names.arg)     
    
    graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
              graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
              locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
              legmap = legmap2, legends = legends2, bin = is.norm)
    
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
        
        if(is.na(num_supp)) {
          dev.new(noRStudioGD = FALSE)
          num_supp <<- dev.list()[length(dev.list())]
        }
        
        graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)], 
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col[1], symbol = pch[1], labvar = c(varChoice1, varChoice2))
      }
    } else {
      tkmessageBox(message = "List of Variables and list of variables names must have been given",
                   icon ="warning", type = "ok")
    }
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  
  cartfunc <- function() {  
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
            W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
            nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
            label = label, cex.lab = cex.lab, labmod = names.arg)     
    } else {
      tkmessageBox(message = "Spatial contours have not been given",
                   icon = "warning", type = "ok")    
    }
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc<-function() {
    obs <<- vector(mode = "logical", length = length(long));
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
          carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
          W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
          nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
          label = label, cex.lab = cex.lab, labmod = names.arg)       
    
    graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
              graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
              locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
              legmap = legmap2, legends = legends2, bin = is.norm)
    
    if ((graphChoice != "") && (varChoice1 != ""))
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)], 
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col[1], symbol = pch[1], labvar = c(varChoice1, varChoice2))
    
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
    graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = dev.list()[length(dev.list())],
              graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
              locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
              legmap = legmap2, legends = legends2, bin = is.norm)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())],
          carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
          W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
          nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
          label = label, cex.lab = cex.lab, labmod = names.arg) 
    dev.off()
    
    if(!is.na(num_supp)) {
      fig_supp <- "fig_supp_GeoXp.pdf"
      k <- 1
      while(file.exists(fig_supp)) {
        fig_supp <- paste("fig_supp_GeoXp", "_", k, ".pdf", sep = "")
        k <- k + 1
      }
      pdf(fig_supp)
      graphique(var1 = listvar[, which(listnomvar == varChoice1)], var2 = listvar[,which(listnomvar == varChoice2)], 
                obs = obs, num = dev.list()[length(dev.list())], graph = graphChoice, 
                couleurs = col[1], symbol = pch[1], labvar = c(varChoice1, varChoice2))
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
  
  fnointer<-function() {
    if (length(criteria) != 0) {
      nointer <<- !nointer
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
            W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
            nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
            label = label, cex.lab = cex.lab, labmod = names.arg)      
    } else {
      tkmessageBox(message = "Criteria has not been given", icon = "warning", type = "ok")
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
          carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
          W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
          nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
          label = label, cex.lab = cex.lab, labmod = names.arg)   
    
  }
  
  ####################################################
  # Bubble  of Lisa
  ####################################################
  
  lisa <- function() {        
    res3 <- choix.bubble(buble2, abs(ilocal), "ilocal", legends2, num_graph, num_carte)
    
    buble2 <<- res3$buble
    legends2 <<- res3$legends
    z2 <<- res3$z
    legmap2 <<- res3$legmap
    
    graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
              graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
              locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
              legmap = legmap2, legends = legends2, bin = is.norm)
    
  }
  
  ####################################################
  # Permutation
  ####################################################
  permutation <- function() {
    tt1 <- tktoplevel()
    Name <- tclVar("n")
    entry.Name <- tkentry(tt1, width = "5", textvariable = Name)
    tkgrid(tklabel(tt1, text = "Number of simulations"), entry.Name)
    
    OnOK <- function() { 
      value1 <- tclvalue(Name)
      tkdestroy(tt1)
      if(is.na(as.integer(value1))) {
        tkmessageBox(message = "Sorry, but you have to choose decimal values",
                     icon = "warning", type = "ok")
      } else {
        n <- as.integer(value1)
        perm <- NULL
        for (i in 1:n) {
          sam <- sample(var, length(var))
          sam <- sam - mean(sam)
          
          epe <- sam %*% sam
          mi <- (sam %*% W %*% sam)/epe
          morani <- round(mi, 4)
          perm <- c(perm, morani)
        }
        msg <- paste("The p-value of the permutation test is :", 
                     1 - length(which(perm < rep(MORAN, n)))/n)
        tkmessageBox(message = msg, icon = "info", type = "ok") 
      }   
    }
    
    
    OK.but <-tkbutton(tt1, text = "   OK   ", command = OnOK)
    tkgrid(OK.but)
    tkfocus(tt1)
  }
  
  #################################################
  ###########   Representation
  

  carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            carte = carte, nocart = nocart, classe = obsq, couleurs = col3, symbol = pch2,
            W = W, method = method, buble = buble, cbuble = z, criteria = criteria,
            nointer = nointer, legmap = legmap, legends = legends, axis = axes, lablong = lablong, lablat = lablat,
            label = label, cex.lab = cex.lab, labmod = names.arg)   
      
  graphique(var1 = var, var2 = wvar, var3 = ilocal, obs = obs, num = num_graph,
                graph = graph, labvar = labvar, couleurs = col2, symbol = pch2, 
                locmoran = locmoran, obsq = obsq, cex.lab = cex.lab, buble = buble2, cbuble = z2, 
                legmap = legmap2, legends = legends2, bin = is.norm)
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "moranplotmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pt1func)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = poly1func)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",fill = "x")
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the graphic", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point2.but <- tkbutton(frame1c, text = "Selection by point", command = pt2func)
    poly2.but <- tkbutton(frame1c, text = "Selection by polygon ", command = poly2func)
    tkpack(point2.but, poly2.but, side = "left", expand = "TRUE",fill = "x")
    tkpack(frame1c, expand = "TRUE", fill = "x")
    
    frame1d <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1d, text = "Select a quadrant", font = "Times 12",
                   foreground = "darkred", background = "white"))
    HH.but <- tkbutton(frame1d, text = "H.-H.", command = quad1func) 
    LH.but <- tkbutton(frame1d, text = "H.-L.", command = quad4func)
    LL.but <- tkbutton(frame1d, text = "L.-L.", command = quad3func)
    HL.but <- tkbutton(frame1d, text = "L.-H.", command = quad2func)
    tkpack(HH.but, LH.but, LL.but, HL.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1d, expand = "TRUE", fill = "x")
    
    frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nettoy.but <- tkbutton(frame1b, text = "     Reset selection     " , command = SGfunc)
    tkpack(nettoy.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1b, expand = "TRUE", fill = "x")
    
    
    frame11a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame11a, text = "Moran test", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(frame11a, expand = "TRUE", fill = "x")
    
    frame11b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    msg <- paste("Moran index ",ifelse(is.norm,"(W normalized)","(W not normalized)"), ": ", 
                 MORAN, " - ", "p-value (Gaussian Test) : ", 
                 ifelse(round(1 - prob.I, 4) < 0.0001, "<0.0001", round(1 - prob.I, 4)))
    tkgrid(tklabel(frame11b, text = msg), columnspan = 2)
    tkpack(frame11b, expand = "TRUE", fill = "x")
    
    frame11c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    noint10.but <- tkbutton(frame11c, text = "Permutation Test", command = permutation)
    tkpack(noint10.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame11c, expand = "TRUE", fill = "x")
    
    frame2 <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2, text = "Options", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame2, text = "Spatial contours", font = "Times 10",
                   foreground = "darkred", background = "white"),
           tklabel(frame2, text = "Preselected sites", font = "Times 10",
                   foreground = "darkred", background = "white"),
           tklabel(frame2, text = "Bubbles", font = "Times 10",
                   foreground = "darkred", background = "white"),
           tklabel(frame2, text = "Additional graph", font = "Times 10",
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = fbubble)
    autre.but <- tkbutton(frame2b, text = "     OK     " , command = graphfunc)
    tkpack(nocou1.but,noint1.but,bubble.but,autre.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    
    
    frame2e <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2e, text = "Print different colors by quadrant  ", font = "Times 10",
                   foreground = "darkred", background = "white"),
           tklabel(frame2e, text = "Bubles of LISA", font = "Times 10",
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2e, expand = "TRUE", fill = "x")
    
    frame2f <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    colquad.but <- tkbutton(frame2f, text = "On/Off", command = colfunc)
    lisa.but <- tkbutton(frame2f, text = "     OK     " , command = lisa)
    tkpack(colquad.but,lisa.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2f, expand = "TRUE", fill = "x")
    
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

