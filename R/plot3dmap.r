plot3dmap <- function(sf.obj, names.var, box = TRUE, criteria = NULL, carte = NULL, 
                      identify = NULL, cex.lab = 0.8, pch = 16, col = "lightblue3", xlab = "", ylab = "", zlab = "",
                      axes = FALSE, lablong = "", lablat = "") {
  
  ###################################################
  ########## COMMON to ALL FUNCTIONS in GeoXp
  
  envir <- globalenv()
  # Verification of the Spatial Object sf.obj
  class.obj <- class(sf.obj)[1]
  
  if(class.obj != "sf") 
    stop("sf.obj may be a sf object")
  
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
  
  if (is.numeric(names.var)) {
    if (all(names.var <= ncol(listvar))) 
      names.var <- listnomvar[names.var]
    else
      stop("Dimension of names.var is not good")
  }
  
  if(!(all(names.var %in% names(sf.obj))))
    stop("names.var is not included in the sf object")
  
  var1 <- sf.obj[[names.var[1]]] 
  var2 <- sf.obj[[names.var[2]]] 
  var3 <- sf.obj[[names.var[3]]] 
  
  var1 <- as.matrix(var1)
  var2 <- as.matrix(var2)
  var3 <- as.matrix(var3)
  lat <- as.matrix(lat)
  long <- as.matrix(long)
  obs <- vector(mode = "logical", length = length(long))
  
  # for colors in map and new grahics
  col2 <- "blue"
  col3 <- col[1]
  pch2 <- pch[1]
  labmod <- ""
  
  # Windows device
  if(length(dev.list()) == 0 & options("device") == "RStudioGD")
    dev.new()
  # for map
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
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while(!quit) {
      #selection des points
      
      dev.set(num_carte)
      
      if (nrow(sf.obj) > 100 & st_geometry_type(sf.obj, by_geometry = F) == "POLYGON" & !buble) {
        points(long, lat, pch = 16, col = "royalblue")
      }
      
      loc <- locator(1)
      
      if (is.null(loc)) {
        quit <- TRUE
        carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
              criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
              carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
              lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
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
      
      # graphiques
      if (length(which(obs == TRUE)) == 0 || length(which(obs == TRUE)) == length(obs)) {
        if (length(which(obs == TRUE)) == 0)
        plot3d(var1[!obs], var2[!obs], var3[!obs], 
               xlab = xlab, ylab = ylab, zlab = zlab,
               col = col, type = "p", size = 5, box = box)
        else
          plot3d(var1[obs], var2[obs], var3[obs], 
                 xlab = xlab, ylab = ylab, zlab = zlab,
                 col = "red", type = "p", size = 6, box = box)    
      } else { 
        col.prov <- rep(col, length(var1))
        col.prov[obs] <- "red"

        plot3d(var1, var2, var3, 
               xlab = xlab, ylab = ylab, zlab = zlab,
               col = col.prov, type = "p", size = 5, box = box)
      }
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, sf.obj = sf.obj, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
            lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) >= 2))
        graphique(var1 = listvar[,which(listnomvar == varChoice1)], 
                  var2 = listvar[,which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                  symbol = pch, labvar = c(varChoice1, varChoice2))
      
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
      
      obs <<- selectmap(var1 = long, var2 = lat, obs = obs, 
                        Xpoly = polyX, Ypoly = polyY, method = "poly")
      
      # graphiques
      if (length(which(obs == TRUE)) == 0 || length(which(obs == TRUE)) == length(obs)) {
        if (length(which(obs == TRUE)) == 0)
        plot3d(var1[!obs], var2[!obs], var3[!obs],
               xlab = xlab, ylab = ylab, zlab = zlab,
               col = col, type = "p", size = 5, box = box)
        else
          plot3d(var1[obs], var2[obs], var3[obs],
                 xlab = xlab, ylab = ylab, zlab = zlab,
                 col = "red", type = "p", size = 6, box = box)    
      } else { 
        col.prov <- rep(col, length(var1))
        col.prov[obs] <- "red"
        plot3d(var1, var2, var3,
               xlab = xlab, ylab = ylab, zlab = zlab, col = col.prov,
               type = "p", size = 5, box = box)
      }
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, sf.obj = sf.obj, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
            lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)]) 
      
      if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) >= 2))
        graphique(var1 = listvar[,which(listnomvar == varChoice1)], 
                  var2 = listvar[,which(listnomvar == varChoice2)],
                  obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                  symbol = pch, labvar = c(varChoice1, varChoice2))
    }
  }
  
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() { 
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, num = num_carte, sf.obj = sf.obj, buble = buble, cbuble = z,
            criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
            lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
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
        if (((graphChoice == "Histogram") && (!is.numeric(listvar[,which(listnomvar == varChoice1)]))) || 
            ((graphChoice == "Scatterplot") && ((!is.numeric(listvar[, which(listnomvar == varChoice1)])) || 
                                                (!is.numeric(listvar[, which(listnomvar == varChoice2)]))))) {
          tkmessageBox(message = "Variables choosed are not in a good format", icon = "warning", type = "ok")
        } else {
          res1 <- choix.couleur(graphChoice, listvar, listnomvar, varChoice1, legends, col, pch,
                                num_carte = num_carte)
          
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
          
          layout(1)
          par(mar = c(5.1, 4.1, 4.1, 2.1))
          graphique(var1 = listvar[,which(listnomvar == varChoice1)], 
                    var2 = listvar[,which(listnomvar == varChoice2)],
                    obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                    symbol = pch, labvar = c(varChoice1, varChoice2))  
          
          carte(long = long, lat = lat, obs = obs, num = num_carte, sf.obj = sf.obj, buble = buble, cbuble = z,
                criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
                carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
                lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
                classe = listvar[, which(listnomvar == varChoice1)]) 
        }
      }   
    } else {
      tkmessageBox(message = "Variables (listvar) and their names (listnomvar) must have been given", icon = "warning", type = "ok")
    }  
  }
  
  
  ####################################################
  # rafraichissement des graphiques
  ####################################################
  
  SGfunc <- function() {
    obs <<- vector(mode = "logical", length = length(long))
    
    # graphiques
    plot3d(var1, var2, var3, 
           xlab = xlab, ylab = ylab, zlab = zlab,
           col = col, type = "p", size = 5, box = box) 
    
    carte(long = long, lat = lat, obs = obs, num = num_carte, sf.obj = sf.obj, buble = buble, cbuble = z,
          criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
          carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
          lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
          classe = listvar[, which(listnomvar == varChoice1)]) 
    
    if ((graphChoice != "") && (varChoice1 != "") && (length(dev.list()) >= 2))
      graphique(var1 = listvar[,which(listnomvar == varChoice1)], 
                var2 = listvar[,which(listnomvar == varChoice2)],
                obs = obs, num = num_supp, graph = graphChoice, couleurs = col3,
                symbol = pch, labvar = c(varChoice1, varChoice2))
    
  }
  
  ####################################################
  # quitter l'application
  ####################################################
  
  quitfunc <- function(){
    tkdestroy(tt)
    assign("GeoXp.open", FALSE, envir = envir)
    close3d()
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
    rgl.postscript(fig_save,"pdf") 
    
    k <- 1
    while(file.exists(map_save)) {
      map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
      k <- k + 1
    }    
    
    pdf(map_save)
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())], buble = buble, cbuble = z,
          criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
          carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
          lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
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
    cat("Figure has been saved in", fig_save, "\n")
    if(!is.na(num_supp))
      cat("Supplemental figure has been saved in", fig_supp, "\n")
    
    assign("last.select", which(obs), envir = envir)
    
    dev.off(num_carte)
    close3d()
    if(!is.na(num_supp))
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
            carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
            lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
            classe = listvar[, which(listnomvar == varChoice1)])  
    } else {
      tkmessageBox(message = "Criteria has not been given", icon = "warning", type = "ok")
    }
  }
  
  ####################################################
  # Bubble
  ####################################################
  
  fbubble<-function() {
    res2<-choix.bubble(buble,listvar,listnomvar,legends, num_carte = num_carte)
    
    buble <<- res2$buble
    legends <<- res2$legends
    z <<- res2$z
    legmap <<- res2$legmap
    
    carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
          criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
          carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
          lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
          classe = listvar[, which(listnomvar == varChoice1)])  
    
  }
  
  
  ####################################################
  # graphiques 
  ####################################################
  
  carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte, buble = buble, cbuble = z,
        criteria = criteria, nointer = nointer, label = label, symbol = pch2, couleurs = col2,
        carte = carte, nocart = nocart, legmap = legmap, legends = legends, axis = axes, labmod = labmod,
        lablong = lablong, lablat = lablat, cex.lab = cex.lab, method = method,
        classe = listvar[, which(listnomvar == varChoice1)]) 
  
  par3d(windowRect = c(0, 32, 512, 544), userMatrix = rotationMatrix(5*pi/180, 0, 1, 0) %*% par3d("userMatrix"), family = "sans")
  bg3d(color = "white")
  
  plot3d(var1, var2, var3, 
         xlab = xlab, ylab = ylab, zlab = zlab,
         xlim = range(var1), ylim = range(var2), zlim = range(var3), 
         col = col, size = 5, type = "n", box = box) 
  
  plot3d(var1, var2, var3, 
         xlab = xlab, ylab = ylab, zlab = zlab,
         xlim = range(var1), ylim = range(var2), zlim = range(var3),
         col = col, size = 5, type = "p", box = FALSE, add = TRUE)
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "plot3dmap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Work on the map", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text = "Selection by point", command = pointfunc)
    poly.but <- tkbutton(frame1a, text = "Selection by polygon ", command = polyfunc)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE", fill = "x")
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    
    frame1b <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    nettoy.but <- tkbutton(frame1b, text = "     Reset selection     " , command=  SGfunc)
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
    nocou1.but <- tkbutton(frame2b, text="On/Off", command=cartfunc)
    noint1.but <- tkbutton(frame2b, text="On/Off", command=fnointer)
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

