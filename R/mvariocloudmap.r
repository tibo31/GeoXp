mvariocloudmap <- function(sp.obj, nb.obj, names.var, quantiles = TRUE, names.attr = names(sp.obj), 
                           criteria = NULL, carte = NULL, identify = FALSE, cex.lab = 0.8, pch = 16, col = "lightblue3", 
                           xlab = "Pairwise spatial distances", ylab = "Pairwise Mahalanobis distances", axes = FALSE,
                           lablong = "", lablat = "") {
  
  # Verification of the Spatial Object sp.obj
  class.obj <- class(sp.obj)[1]
  
  if (substr(class.obj, 1, 7) != "Spatial")
    stop("sp.obj may be a Spatial object")
  
  if (substr(class.obj, nchar(class.obj) - 8, nchar(class.obj)) != "DataFrame")
    stop("sp.obj should contain a data.frame")
  
  if (!is.numeric(names.var) &
      length(match(names.var, names(sp.obj))) != length(names.var))
    stop("At least one component of names.var is not included in the data.frame of sp.obj")
  
  if (length(names.attr) != length(names(sp.obj)))
    stop("names.attr should be a vector of character with a length equal to the number of variable")
  
  # Is there a Tk window already open ?
  if (interactive()) {
    if (!exists("GeoXp.open", envir = globalenv())) {
      assign("GeoXp.open", TRUE, envir = globalenv())
    } else {
      if (get("GeoXp.open", envir = globalenv())) {
        stop(
          "A GeoXp function is already open. Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics"
        )
      } else
        assign("GeoXp.open", TRUE, envir = globalenv())
    }
  }
  
  # we propose to refind the same arguments used in first version of GeoXp
  long <- coordinates(sp.obj)[, 1]
  lat <- coordinates(sp.obj)[, 2]
  
  dataset <- sp.obj@data[, names.var]
  
  listvar <- sp.obj@data
  listnomvar <- names.attr
  
  # Code which was necessary in the previous version
  if (is.null(carte) &
      class.obj == "SpatialPolygonsDataFrame")
    carte <- spdf2list(sp.obj)$poly
  
  # for identifying the selected sites
  if(identify)
    label <- row.names(listvar)
  else
    label <- ""
  
  # Spatial weight matrix
  W <- nb2mat(nb.obj)
  
  # initialisation
  xy <- cbind(long, lat)
  nointer <- FALSE
  nocart <- FALSE
  buble <- FALSE
  legends <- list(FALSE, FALSE, "", "")
  z <- NULL
  legmap <- NULL
  inout <- NULL
  labvar <- c(xlab, ylab)
  obs <- matrix(FALSE, nrow = length(long), ncol = length(long))
  obs2 <- matrix(FALSE, nrow = length(long), ncol = length(long))
  graf <- "Neighbourplot1"
  names.slide <- c("Quantile smooth spline parameter")
  
  # Transformation d'un data.frame en matrix
  if ((length(listvar) > 0) &&
      (dim(as.matrix(listvar))[2] == 1))
    listvar <- as.matrix(listvar)
  if ((length(dataset) > 0) &&
      (dim(as.matrix(dataset))[2] == 1))
    dataset <- as.matrix(dataset)
  
  # Windows device
  if(length(dev.list()) == 0 & options("device") == "RStudioGD")
    dev.new()
  # for graphic
  dev.new(noRStudioGD = FALSE)
  num_graph <- dev.list()[length(dev.list())]
  # for map
  dev.new(noRStudioGD = FALSE)
  num_carte <- dev.list()[length(dev.list())]

  # calcul des matrices theta et absvar
  n <- nrow(dataset)
  p <- ncol(dataset)
  covr <- covMcd(dataset, alpha = 0.75)
  cinv <- solve(covr$cov)
  
  idx <- matrix(1:n, n, n)
  se <- as.vector(idx[lower.tri(idx)])
  dij <- sqrt((rep(xy[-n, 1], seq(n - 1, 1)) - xy[se, 1]) ^ 2 + 
                (rep(xy[-n, 2], seq(n - 1, 1)) - xy[se, 2]) ^ 2)
  hlp <- as.matrix(dataset[rep(1:(n - 1), seq((n - 1), 1)), ] - dataset[se, ])
  MDij <- sqrt(rowSums((hlp %*% cinv) * hlp))
  
  indij <- cbind(rep(1:(n - 1), seq(n - 1, 1)), se)
  
  theta <- matrix(0, n, n)
  theta[indij] <- dij
  theta <- theta + t(theta)
  
  absvar <- matrix(0, n, n)
  absvar[indij] <- MDij
  absvar <- absvar + t(absvar)
  
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
  
  borne1 <- 0.01
  borne2 <- 0.99
  alpha <- 0.5
  
  temp_data <- data.frame(
    var1 = sort(theta),
    var2 = absvar[order(theta)])
  alpha1 <- qgam(var2 ~ s(var1, k = 20, bs = "ad"), data = temp_data, qu = alpha)
  
  ####################################################
  # selection d'un point sur la carte
  ####################################################
  
  pointfunca <- function() {
    if (graf == "pairwise") 
      SGfunc()
    
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
        quit <- TRUE
        carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
              nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
              lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
              axis = axes, legmap = legmap, legends = legends)       
        next
      }
      obs2 <<- selectmap(var1 = long, var2 = lat, obs = obs2,
                         Xpoly = loc[1], Ypoly = loc[2], method = "point")
      
      obs <<- (diag(n)*obs2 > 0)
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
            lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
            axis = axes, legmap = legmap, legends = legends)    
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
                quantiles = quantiles, alpha1 = alpha1)
    }
  }
  
  ####################################################
  # selection d'un polygone
  ####################################################
  
  polyfunca <- function() {
    if (graf=="pairwise") SGfunc()
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
        quit<-TRUE
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
      
      obs2 <<- selectmap(var1 = long, var2 = lat, obs = obs2, 
                         Xpoly = polyX, Ypoly = polyY, method = "poly")
      obs <<- (diag(n)*obs2 > 0)
      
      # graphiques
      carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
            lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
            axis = axes, legmap = legmap, legends = legends)    
      
      graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
                quantiles = quantiles, alpha1 = alpha1)
    }
  }
  
  ####################################################
  # selection d'un point sur l'angleplot
  ####################################################
  
  pointfunc <- function() {
    if (graf == "Neighbourplot1") SGfunc()
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
        graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
                  graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
                  quantiles = quantiles, alpha1 = alpha1)
        next
      }
      
      obs <<- selectstat(var1 = theta, var2 = absvar, obs = obs,Xpoly = loc[1], Ypoly = loc[2],
                         method = "AnglePoint",long = long, lat = lat)
      
      diag(obs) <<- FALSE
      
      obs2 <<- obs
      graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
                quantiles = quantiles, alpha1 = alpha1)
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
            lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
            axis = axes, legmap = legmap, legends = legends)    
    }
  }
  
  ####################################################
  # selection d'un polygone sur l'angleplot
  ####################################################
  
  polyfunc <- function() {
    if (graf=="Neighbourplot1") 
      SGfunc()
    graf<<-"pairwise"
    
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
        def <- inout(cbind(theta[, i], absvar[, i]), 
                     cbind(polyX, polyY), bound = TRUE)       
        obs[def, i] <<- !obs[def, i]    
      }
      
      obs2 <<- obs
      graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
                graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
                quantiles = quantiles, alpha1 = alpha1)
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
            lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
            axis = axes, legmap = legmap, legends = legends)    
    }
  }    
  
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() { 
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
            lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
            axis = axes, legmap = legmap, legends = legends)    
    } else {
      tkmessageBox(message = "Spatial contours have not been given", icon = "warning", type = "ok")    
    }
  }
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- matrix(FALSE, nrow = length(long), ncol = length(long))
    obs2 <<- matrix(FALSE, nrow = length(long), ncol = length(long))
    
    carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
          nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
          lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
          axis = axes, legmap = legmap, legends = legends)    
    
    graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
              graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
              quantiles = quantiles, alpha1 = alpha1)
  }
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function() {
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = globalenv())
    dev.off(num_graph)
    dev.off(num_carte)
  }
  
  quitfunc2<-function(){
    
    fig_save <- "fig_GeoXp.pdf"
    map_save <- "map_GeoXp.pdf"
    k <- 1
    while(file.exists(fig_save)) {
      fig_save <- paste("fig_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }
    pdf(fig_save)
    graphique(var1 = theta, var2 = absvar, obs = obs2, num = dev.list()[length(dev.list())], 
              graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
              quantiles = quantiles, alpha1 = alpha1)
    dev.off()
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, num = dev.list()[length(dev.list())], buble = buble, criteria = criteria,
          nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
          lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
          axis = axes, legmap = legmap, legends = legends)    
    dev.off()

    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = globalenv())
    cat("Results have been saved in last.select object \n")
    cat("Map has been saved in", map_save, "\n")
    cat("Figure has been saved in", fig_save, "\n")
    
    assign("last.select", which(obs), envir = globalenv())
    
    dev.off(num_carte)
    dev.off(num_graph)
  }
  
  
  ####################################################
  # Open a no interactive selection
  ####################################################
  
  fnointer <- function() {
    if (length(criteria) != 0) {
      nointer <<- !nointer
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
            nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
            lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
            axis = axes, legmap = legmap, legends = legends)    
    } else {
      tkmessageBox(message = "Criteria has not been given", icon = "warning", type = "ok")
    }
    
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble <- function() {
    res2 <- choix.bubble(buble, cbind(chi2.quant,listvar), c("chi2.quant",listnomvar), legends, num_graph, num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    
    if (legends[[1]]) {
      if ((legmap[length(legmap)] == "chi2.quant")) {
        legmap <<- c(legmap, paste(">", round(alphab[1], 2)), paste(round(alphab[2], 2), "-", round(alphab[1], 2)), 
                     paste(round(alphab[3], 2), "-", round(alphab[2], 2)), paste(round(alphab[4], 2), "-", round(alphab[3],2)),
                     paste("<", round(alphab[4], 2)), "Mahalanobis")
      }
    }
  
    carte(long = long, lat = lat, obs = obs, num = num_carte, buble = buble, criteria = criteria,
          nointer = nointer, cbuble = z, carte = carte, nocart = nocart, lablong = lablong,
          lablat = lablat, label = label, cex.lab = cex.lab, symbol = pch, method = "pairwise",
          axis = axes, legmap = legmap, legends = legends)    
  }
  
  ####################################################
  # Pour le alpha 
  ####################################################
  
  refresh.code <- function(...)  {

    alpha <<- slider1(names.slide = names.slide, no = 1)
    
    temp_data <- data.frame(
      var1 = sort(theta),
      var2 = absvar[order(theta)])
    
    alpha1 <<- qgam(var2 ~ s(var1, k = 20, bs = "ad"), data = temp_data, qu = alpha)
    
    graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
              graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
              quantiles = quantiles, alpha1 = alpha1)
  }
  
  ####################################################
  # Representation graphique
  ####################################################
  
  carte(long = long, lat = lat, obs = obs, num = num_carte, lablong = lablong, 
        lablat = lablat, label = label, cex.lab = cex.lab, 
        symbol = pch, carte = carte, method = "pairwise",
        axis = axes, legends = legends)
  
  graphique(var1 = theta, var2 = absvar, obs = obs2, num = num_graph, 
            graph = "pairwise", labvar = labvar, couleurs = col, symbol = pch, 
            quantiles = quantiles, alpha1 = alpha1)
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "mvariocloudmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pointfunca)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = polyfunca)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",fill = "x")
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1c, text = "Work on the graph", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but2 <- tkbutton(frame1c, text = "Selection by point", command = pointfunc)
    poly.but2 <- tkbutton(frame1c, text = "Selection by polygon ", command = polyfunc)
    tkpack(point.but2, poly.but2, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1c, expand = "TRUE", fill = "x")
    
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
                   foreground = "darkred", background = "white"), side = "left", fill="x", expand = "TRUE")
    tkpack(frame2, expand = "TRUE", fill = "x")
    
    frame2b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nocou1.but <- tkbutton(frame2b, text = "On/Off", command = cartfunc)
    noint1.but <- tkbutton(frame2b, text = "On/Off", command = fnointer)
    bubble.but <- tkbutton(frame2b, text="On/Off", command=fbubble)
    tkpack(nocou1.but,noint1.but,bubble.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame2b, expand = "TRUE", fill = "x")
    
    if (length(quantiles) != 0) {
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

