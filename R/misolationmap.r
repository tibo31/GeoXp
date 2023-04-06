misolationmap <- function(sf.obj, nb.obj, names.var, propneighb = 0.4, chisqqu = 0.975,
                          criteria = NULL, carte = NULL, identify = NULL, cex.lab = 0.8, 
                          pch = 16, col = "lightblue3", xlab = "degree of isolation", 
                          ylab = "Pairwise Mahalanobis distances", 
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
  if (is.numeric(names.var)) {
    if (all(names.var <= ncol(listvar))) 
      names.var <- listnomvar[names.var]
    else
      stop("Dimension of names.var is not good")
  }
  
  if(!(all(names.var %in% names(sf.obj))))
    stop("names.var is not included in the sf object")
  
  dataset <- listvar[, names.var]
  
  obs <- matrix(FALSE, nrow = length(long), ncol = length(long))
  graf <- "Neighbourplot1"
  obsref <- obs
  
  outselect = TRUE
  outselect2 = TRUE
  orderselect = TRUE
  
  W <- nb2mat(nb.obj)
  Wref <- W
  
  # calcul des matrices theta et absvar
  n <- nrow(dataset)
  p <- ncol(dataset)
  
  covr <- covMcd(dataset, alpha = 0.75)
  cinv <- solve(covr$cov)
  MDglobal <- sqrt(mahalanobis(dataset, covr$center, cinv, inverted = TRUE))
  
  # TRUE/FALSE for non-outlying/outlying:
  qchi <- sqrt(qchisq(chisqqu, p))
  MDglobalTF <- (MDglobal < qchi)
  
  idx <- matrix(1:n, n, n)
  se <- as.vector(idx[lower.tri(idx)])
  hlp <- as.matrix(dataset[rep(1:(n - 1), seq((n - 1), 1)), ] - dataset[se, ])
  MDij <- sqrt(rowSums((hlp %*% cinv) * hlp))
  
  MDpair <- matrix(0, n, n)
  MDpair[lower.tri(MDpair)] <- MDij
  MDpair <- t(MDpair)
  MDpair[lower.tri(MDpair)] <- MDij
  
  MDpairN <- vector("list", n)
  
  # boundary that should include required proportion of neighbors:
  chibound <- rep(NA, n)
  theta <- matrix(0, n, n)
  absvar <- matrix(0, n, n)
  
  for (i in 1:n) {
    MDpairN[[i]] <- MDpair[i, nb.obj[[i]]]
    nn <- max(1, round(length(MDpairN[[i]]) * propneighb)) # number of neighbors in tolerance ellipse
    nval <- MDpairN[[i]][order(MDpairN[[i]])][nn]  # value of largest neighbor to be included
    chibound[i] <- pchisq(nval ^ 2, p, MDglobal[i] ^ 2)
    
    theta[i, nb.obj[[i]]] <- chibound[i]
    absvar[i, nb.obj[[i]]] <- MDpair[i, nb.obj[[i]]]
  }
  
  thetaref <- theta
  absvaref <- absvar
  
  thetaref2 <- thetaref
  absvaref2 <- absvaref
  
  # sort according to values of "chibound" - separately for outliers and non-outliers
  idx1 <- order(chibound[MDglobalTF])
  idx0 <- order(chibound[!MDglobalTF])
  
  idxg <- sort(chibound, index.return = TRUE)
  
  # calcul des distances de Mahalanobis par site
  rd <- sqrt(mahalanobis(dataset, center = covr$center, cov = covr$cov))
  
  pcrit <- ifelse(p <= 10, (0.24 - 0.003 * p) / sqrt(n), (0.252 - 0.0018 * p) / sqrt(n))
  delta <- qchisq(1 - 0.025, p)
  
  d2 <- mahalanobis(dataset, covr$center, covr$cov)
  d2ord <- sort(d2)
  dif <- pchisq(d2ord, p) - (0.5:n) / n
  i <- (d2ord >= delta) & (dif > 0)
  alfan <- ifelse(sum(i) == 0, 0, max(dif[i]))
  if (alfan < pcrit)
    alfan <- 0
  cn <- ifelse(alfan > 0, max(d2ord[n - ceiling(n * alfan)], delta), Inf)
  
  alphab <- ifelse(cn != Inf, sqrt(c(cn, qchisq(c(0.75, 0.5, 0.25), ncol(dataset)))), 
                   sqrt(qchisq(c(0.975, 0.75, 0.5, 0.25), ncol(dataset))))
  
  chi2.quant <- rep(0, n)
  lalpha <- length(alphab)
  for (j in 1:lalpha) {
    if (j == 1) {
      (chi2.quant[which(rd >= alphab[j])] <- lalpha)
    } else {
      chi2.quant[which((rd < alphab[j - 1]) &
                         (rd >= alphab[j]))] <- lalpha + 1 - j
      
    }
  }
  
  ####################################################
  # selection d'un point sur la carte
  ####################################################
  
  pointfunca <- function() {
    if (graf == "pairwise") SGfunc()
    graf <<- "Neighbourplot1"
    
    quit <- FALSE
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while (!quit) {
      dev.set(num_carte)
      loc <- locator(1)
      if (is.null(loc)) {
        quit<-TRUE
        carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
              buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
              nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
              symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
        next
      }
      obs2 <- selectmap(var1 = long, var2 = lat, obs = obs,
                        Xpoly = loc[1], Ypoly = loc[2], method = "point")
      
      obs <<- (W*obs2 > 0)
      
      if (!outselect)
        obs[MDglobalTF, ] <<- FALSE
      
      if (!outselect2)
        obs[!MDglobalTF, ] <<- FALSE
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
            symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, 
                symbol = pch, direct = propneighb)
    }
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunca <- function() {
    if (graf=="pairwise") 
      SGfunc()
    
    graf <<- "Neighbourplot1"
    
    polyX <- NULL
    polyY <- NULL
    quit <- FALSE
    
    dev.set(num_carte)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
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
      
      obs2 <- selectmap(var1 = long, var2 = lat, obs = obs, 
                        Xpoly = polyX, Ypoly = polyY, method = "poly")
      obs <<- (W*obs2 > 0)
      
      if (!outselect)
        obs[MDglobalTF, ] <<- FALSE
      
      if(!outselect2)
        obs[!MDglobalTF, ] <<- FALSE
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
            symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
      
      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, graph = "pairwise", labvar = labvar,
                couleurs = col, symbol = pch, direct = propneighb)
      
    }
  }
  
  
  ####################################################
  # selection d'un point sur l'angleplot
  ####################################################
  
  pointfunc <- function() {
    if (graf == "Neighbourplot1") 
      SGfunc()
    
    graf <<- "pairwise"
    
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
        graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, graph = "pairwise", labvar = labvar,
                  couleurs = col, symbol = pch, direct = propneighb)
        next
      }
      
      obs <<- selectstat(var1 = theta, var2 = absvar, obs = obs,Xpoly = loc[1], Ypoly = loc[2],
                         method = "AnglePoint",long = long, lat = lat)
      
      if (!outselect)
        obs[MDglobalTF, ] <<- FALSE
      
      if (!outselect2)
        obs[!MDglobalTF, ] <<- FALSE
  
      diag(obs) <<- FALSE
      
      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, 
                symbol = pch, direct = propneighb)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
            symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
    }
  }
  
  ####################################################
  # selection des outliers
  ####################################################
  
  outlier <- function() {
    
    outselect <<- !outselect
    labvar <- c(xlab,ylab)
    
    if (outselect) { 
      tkmessageBox(message = "You have remited non global outliers")
      theta <<- thetaref
      absvar <<- absvaref
      thetaref2 <<- thetaref
      absvaref2 <<- absvaref
      if (!orderselect) {
        labvar <- c("Rank",ylab)
        ind0 <- which(theta == 0, arr.ind = TRUE)
        theta[idxg$ix, ] <<- t(matrix(rep(1:n, each = length(idxg$ix)),
                                      length(idxg$ix), length(idxg$ix)))
        theta[ind0] <<- 0
        theta <<- theta
      }
    } else {
      tkmessageBox(message = "You have selected only global outliers")
      SGfunc()
      if (!outselect2) {
        theta <<- thetaref
        absvar <<- absvaref
        if (!orderselect) {
          labvar <- c("Rank", ylab)
          ind0 <- which(theta == 0, arr.ind = TRUE)
          theta[idxg$ix, ] <<- t(matrix(rep(1:n, each = length(idxg$ix)),
                                        length(idxg$ix), length(idxg$ix)))
          theta[ind0] <<- 0
          theta <<- theta
        }
      }

      theta[MDglobalTF, ] <<- 0
      absvar[MDglobalTF, ] <<- 0
      
      thetaref2 <<- theta
      absvaref2 <<- absvar
    }
    
    graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
              graph = "pairwise", labvar = labvar, couleurs = col, 
              symbol = pch, direct = propneighb)
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
          symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
  }
  
  ####################################################
  # selection des non outliers
  ####################################################
  
  outlier2 <- function() {
    
    outselect2 <<- !outselect2
    labvar <- c(xlab, ylab)
    
    if (!outselect) 
      tkmessageBox(message="You have selected non global outliers")
    
    if (outselect) 
      tkmessageBox(message="You have remited global outliers")
    
    if (outselect2) {
      theta <<- thetaref
      absvar <<- absvaref
      thetaref2 <<- thetaref
      absvaref2 <<- absvaref
      if (!orderselect) {
        labvar <- c("rank", ylab)
        ind0 <- which(theta == 0, arr.ind = TRUE)
        theta[idxg$ix, ] <<- t(matrix(rep(1:n, each = length(idxg$ix)),
                                      length(idxg$ix),length(idxg$ix)))
        theta[ind0] <<- 0
        theta <<- theta
      }
    } else {
      SGfunc()

      if (!outselect) {
        theta <<- thetaref
        absvar <<- absvaref
        if(!orderselect) {
          labvar <- c("Rank", ylab)
          ind0 <- which(theta == 0, arr.ind = TRUE)
          theta[idxg$ix, ] <<- t(matrix(rep(1:n, each = length(idxg$ix)),
                                        length(idxg$ix),length(idxg$ix)))
          theta[ind0] <<- 0
          theta <<- theta
        }
      }
      theta[!MDglobalTF, ] <<- 0
      absvar[!MDglobalTF, ] <<- 0
      
      thetaref2 <<- theta
      absvaref2 <<- absvar
    }
    
    graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
              graph = "pairwise", labvar = labvar, couleurs = col, 
              symbol = pch, direct = propneighb)
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
          symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
  }
  
  ####################################################
  # trier les observations
  ####################################################
  
  ordering <- function() {
    orderselect <<- !orderselect
    
    if (orderselect) 
      tkmessageBox(message="Degree of isolation on the x-axis")
    
    if (!orderselect) 
      tkmessageBox(message="Rank on the x-axis")
    
    if (orderselect) {
      theta <<- thetaref
      absvar <<- absvaref
    } else {
      labvar <- c("Rank", ylab)
      ind0 <- which(theta == 0, arr.ind = TRUE)
      theta[idxg$ix, ] <<- t(matrix(rep(1:n, each = length(idxg$ix)),
                                    length(idxg$ix), length(idxg$ix)))
      
      theta[ind0] <<- 0
      theta <<- theta
    }
    
    graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
              graph = "pairwise", labvar = labvar, couleurs = col, 
              symbol = pch, direct = propneighb)
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
          symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
  }
  
  ####################################################
  # selection d'un polygone sur l'angleplot
  ####################################################
  
  polyfunc <- function() {
    if (graf == "Neighbourplot1") 
      SGfunc()
    graf <<- "pairwise"
    
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
        
        my_points <- st_as_sf(data.frame(x = theta[, i], y = absvar[, i]), coords = c("x", "y"))
        polyg <- cbind(unlist(polyX), unlist(polyY))
        pol <- st_sfc(st_polygon(list(polyg)))
        def <- as.vector(st_intersects(my_points, pol, sparse = FALSE))
        obs[def, i] <<- !obs[def, i]    
      }
      
      if (!outselect)
        obs[MDglobalTF, ] <<- FALSE
      
      if (!outselect2)
        obs[!MDglobalTF,]<<-FALSE

      graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, 
                symbol = pch, direct = propneighb)
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
            symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
    }
  }    
  
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() { 
    if (length(carte) != 0) {
      nocart <<- !nocart

      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
            symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
    } else {
      tkmessageBox(message = "Spatial contours have not been given",
                   icon = "warning", type = "ok")    
    }
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- matrix(FALSE, nrow = length(long), ncol = length(long))
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
          symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
    
    graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
              graph = "pairwise", labvar = labvar, couleurs = col, 
              symbol = pch, direct = propneighb)
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
              graph = "pairwise", labvar = labvar, couleurs = col, 
              symbol = pch, direct = propneighb)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())], 
          buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
          symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
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
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
            buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
            nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
            symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends) 
    } else {
      tkmessageBox(message = "Criteria has not been given",
                   icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function() {
    res2<-choix.bubble(buble, cbind(chi2.quant, listvar), c("chi2.quant",listnomvar), legends, num_graph, num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    if (legends[[1]]) {
      if ((legmap[length(legmap)] == "chi2.quant")) {
        legmap <<- c(legmap, paste(">", round(alphab[1], 2)),
                     paste(round(alphab[2], 2), "-", round(alphab[1], 2)),
                     paste(round(alphab[3], 2), "-", round(alphab[2], 2)),
                     paste(round(alphab[4], 2), "-", round(alphab[3], 2)), 
                     paste("<", round(alphab[4], 2)), "Mahalanobis")
      }
    }
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, 
          buble = buble, criteria = criteria, nointer = nointer, cbuble = z, carte = carte,
          nocart = nocart, lablong = lablong, lablat = lablat, label = label, cex.lab = cex.lab,
          symbol = pch, method = "pairwise", axis = axes, legmap = legmap, legends = legends)
  }

  ####################################################
  # Representation graphique
  ####################################################
  
  carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, lablong = lablong, lablat = lablat, 
        label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
        axis = axes, legends = legends)
  
  graphique(var1 = theta, var2 = absvar, obs = obs, num = num_graph, 
            graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
            direct = propneighb)
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "misolationmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pointfunca)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon", command = polyfunca)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",
           fill = "x")
    
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the graphic", font = "Times 12",
                   foreground = "darkred", background = "white"))
    intervalle1.but <- tkbutton(frame1c, text = "Selection by point", command = pointfunc)
    intervalle11.but <- tkbutton(frame1c, text = "Selection by polygon", command = polyfunc)
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
                   foreground = "darkred", background = "white"), 
           tklabel(frame2, text = "Bubbles", font = "Times 11",
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = fbubble)
    tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE",
           fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    frame2c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame2c, text = "Only global outliers ", font = "Times 11",
                   foreground = "darkred", background = "white"),
           tklabel(frame2c, text = "Without global outliers ", font = "Times 11",
                   foreground = "darkred", background = "white"), 
           tklabel(frame2c, text = "Rank on the x-axis", font = "Times 11",
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2c, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = outlier)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = outlier2)
    bubble.but <- tkbutton(frame2b, text = "On/Off", command = ordering)
    tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE",
           fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
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

