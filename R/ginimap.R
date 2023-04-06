ginimap <-function(sf.obj, name.var, criteria = NULL,
                     carte = NULL, identify = NULL, cex.lab = 0.8, pch = 16, col = "lightblue3",
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
  
  # Parametres sur Gini
  graph <- "Lorentz"
  angle <- NULL
  ptX <- NULL
  result <- gini(-var)
  F <- result$F
  G <- result$G
  GINI <- result$gini
  
  # verify the type of the main variable
  if(!(is.integer(var) || is.double(var))) 
    stop("the variable name.var should be a numeric variable")
  # if add a graphic barplot
  labmod <- ""
  # if colors 
  col2 <- "blue"
  col3 <- col[1]
  pch2 <- pch[1]
  
  ####################################################
  # selection d'un point sur la courbe de Lorentz
  ####################################################
  
  ginifunc <- function() {
    SGfunc()
    quit <- FALSE
    ptX <- NULL
    
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
    title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
          cex.sub = 0.8, font.sub = 3, col.sub = "red")
    
    while(!quit) {
      
      dev.set(num_graph)
      loc <- locator(1)
      
      if (is.null(loc)) {
        quit <- TRUE
        graphique(var1 = var, obs = obs, num = num_graph, graph = "Lorentz", 
                  Xpoly = ptX, labvar = labvar, symbol = pch, couleurs = col, F = F, G = G)
        ptX <<- ptX
        graph <<- "Lorentz"
        next
      }
      ptX <- loc[1]
      
      obs <<- selectstat(var1 = var, obs = obs, Xpoly = ptX$x, method = "Lorentz", F = F)
      
      # graphiques
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Lorentz", 
                Xpoly = ptX, labvar = labvar, symbol = pch, couleurs = col, F = F, G = G)
      
      title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")
      title(sub = "To stop selection, click on the right button of the mouse or use ESC", 
            cex.sub = 0.8, font.sub = 3, col.sub = "red")
      
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
  }
  
  ####################################################
  # contour des unites spatiales
  ####################################################
  cartfunc <- function() {
    if (length(carte) != 0) {
      nocart <<- !nocart
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
            symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
            legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
            cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
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
      choix <<- selectgraph(listnomvar, listgraph)
      varChoice1 <<- choix$varChoice1
      varChoice2 <<- choix$varChoice2
      graphChoice <<- choix$graphChoice
      
      if ((graphChoice != "") && (varChoice1 != "")) {
        if (((graphChoice == "Histogram") && (!is.numeric(listvar[, which(listnomvar == varChoice1)]))) || 
            ((graphChoice == "Scatterplot") && ((!is.numeric(listvar[, which(listnomvar == varChoice1)])) || 
                                                (!is.numeric(listvar[,which(listnomvar == varChoice2)]))))) {
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
    
    # graphiques
    
    graphique(var1 = var, obs = obs, num = num_graph, 
              graph = "Lorentz", labvar = labvar, symbol = pch,
              couleurs = col, F = F, G = G)
    
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
  
  quitfunc2<-function() {
    
  fig_save <- "fig_GeoXp.pdf"
  map_save <- "map_GeoXp.pdf"
  k <- 1
  while(file.exists(fig_save)) {
    fig_save <- paste("fig_GeoXp", "_", k, ".pdf", sep = "")
    k <- k + 1
  }
  pdf(fig_save)
  if (graph == "Lorentz")
    graphique(var1 = var, obs = obs, num = dev.list()[length(dev.list())], graph = "Lorentz", 
            Xpoly = ptX, labvar = labvar, symbol = pch, couleurs = col, F = F, G = G)
  else
    graphique(var1 = var, obs = obs, num = dev.list()[length(dev.list())], graph = "VLorentz", Xpoly = angle, labvar = labvar, 
              symbol = pch, couleurs = col, F = F, G = G)
  dev.off()
  
  k <- 1
  while(file.exists(map_save)) {
    map_save <- paste("map_GeoXp", "_", k, ".pdf", sep = "")
    k <- k + 1
  }    
  
  pdf(map_save)
  carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = dev.list()[length(dev.list())],
        buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
        symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
        legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
        cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
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
  # Choisir une valeur
  ####################################################
  
  choixvalue <- function() {
    dev.set(num_graph)
    title("ACTIVE DEVICE", cex.main = 0.8, font.main = 3, col.main = "red")

    tt1 <- tktoplevel()
    Name <- tclVar("value")
    entry.Name <- tkentry(tt1, width = "8", textvariable = Name)
    tkgrid(tklabel(tt1, text = paste("Please enter a value between",
                                     min(var), "and", max(var))),
           entry.Name)
    
    OnOK <- function() {
      angle <- tclvalue(Name)
      tkdestroy(tt1)
      if (is.na(as.numeric(angle)))
        tkmessageBox(message = "Sorry, but you have to choose a decimal value",
                     icon = "warning", type = "ok")
      else {
        msg <- paste("You choose", angle) 
        tkmessageBox(message = msg)
      
      angle <<- as.numeric(angle)
      obs <<- (var <= as.numeric(angle))
      graph <<- "VLorentz"
      graphique(var1 = var, obs = obs, num = num_graph, graph = "VLorentz", Xpoly = as.numeric(angle), labvar = labvar, 
                symbol = pch, couleurs = col, F = F, G = G)
      
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
    }
    
    OK.but <-tkbutton(tt1, text = "   OK   ", command = OnOK)
    tkgrid(OK.but)
    tkfocus(tt1)
  }
  
  ############################
  #### First representation
  
  # Is there a Tk window already open ?
  if (interactive()) {
    if(!exists("GeoXp.open", envir = envir) || length(ls(envir = .TkRoot$env, all.names = TRUE)) == 2) {
      graphique(var1 = var, obs = obs, num = num_graph, graph = "Lorentz", 
                labvar = labvar, symbol = pch, couleurs = col, F = F, G = G)
      
      carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
            buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
            symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
            legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
            cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
      assign("GeoXp.open", TRUE, envir = envir)
    } else {
      if(get("GeoXp.open", envir = envir)) {
        stop("Warning : a GeoXp function is already open. 
             Please, close Tk window before calling a new GeoXp function to avoid conflict between graphics")
        } else {
          graphique(var1 = var, obs = obs, num = num_graph, graph = "Lorentz", 
                    labvar = labvar, symbol = pch, couleurs = col, F = F, G = G)
         
        carte(long = long, lat = lat, obs = obs, sf.obj = sf.obj, num = num_carte,
              buble = buble, cbuble = z, criteria = criteria, nointer = nointer, label = label,
              symbol = pch2, couleurs = col2, carte = carte, nocart = nocart, legmap = legmap,
              legends = legends, axis = axes, labmod = labmod, lablong = lablong, lablat = lablat,
              cex.lab = cex.lab, method = method, classe = listvar[, which(listnomvar == varChoice1)])
        assign("GeoXp.open", TRUE, envir = envir)}
    }
  }
  
  ####################################################
  # creation de la boite de dialogue
  ####################################################
  if (interactive()) {
    fontheading <- tkfont.create(family = "times", size = 14, weight = "bold")
    
    tt <- tktoplevel()
    tkwm.title(tt, "ginimap")
    
    frame1a <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    tkpack(tklabel(frame1a, text = "Interactive selection", font = "Times 14",
                   foreground = "blue", background = "white"))
    tkpack(tklabel(frame1a, text = "Select a point on the Lorentz curve", font = "Times 12",
                   foreground = "darkred", background = "white"))
    point.but <- tkbutton(frame1a, text = "by clicking on graph", command = ginifunc)
    poly.but <- tkbutton(frame1a, text = "by choosing a value", command = choixvalue)
    tkpack(point.but, poly.but, side = "left", expand = "TRUE",
           fill = "x")
    tkpack(frame1a, expand = "TRUE", fill = "x")
    
    frame1c <- tkframe(tt, relief = "groove", borderwidth = 2, background = "white")
    msg <- paste("Gini Index : ", round(abs(GINI), 4))
    tkgrid(tklabel(frame1c, text = msg), columnspan = 2)
    tkgrid(tklabel(frame1c, text = "    "))
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

