genpca <- function(data, w = rep(1 / nrow(data), length = nrow(data)),
           m = diag(ncol(data)), center = NULL, reduc = TRUE) {
    
  #initialisation
    x <- as.matrix(data)
    nr <- nrow(x)
    nc <- ncol(x)
    
    # calcul des differents outputs:
    # inertie, coordonnees des variables et des observations par rapport aux axes de l'ACP
    
    w <- w / sum(w)
    W <- diag(w)
    F1 <- matrix(rep(w, nc), ncol = nc, byrow = FALSE)
    xbar <- colSums(F1 * x)
    
    if (class(center)[1] == "NULL")
      xc <- x - rep(1, nr) %*% t(xbar)
    else {
      if (length(center) == nc)
        xc <- x - rep(1, nr) %*% t(center)
      else
        stop("The vector has not the good dimension")
    }
    
    sig2 <- colSums(F1 * (data ^ 2)) - (xbar ^ 2)
    sigmademi <- diag(sig2 ^ (-1 / 2))
    
    if (reduc)
      xcr <- xc %*% sigmademi
    else
      xcr <- xc
    
    cov <- t(xcr) %*% W %*% xcr
    l <- chol(m)
    covn <- l %*% cov %*% t(l)
    
    res <- eigen(covn)
    U <- res$vectors
    Dv <- res$values
    
    V <- t(l) %*% U
    casecoord <- xcr %*% V
    
    Dvdemi <- Dv ^ (-1 / 2)
    Ddemi <- diag(Dvdemi)
    varcoord <- cov %*% V %*% Ddemi
    inertia <- Dv
    
    return(list(
      inertia = inertia,
      varcoord = varcoord,
      casecoord = casecoord
    ))
  }
