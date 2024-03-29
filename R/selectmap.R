selectmap <- function(var1, var2, obs, Xpoly, Ypoly, method = "") {
  
####################################################
  # selection of a point 
####################################################

  if (method == "point") {
    diff <- abs(var1 - as.numeric(Xpoly)) * (max(var2) - min(var2)) + abs(var2 - as.numeric(Ypoly)) * (max(var1) - min(var1))
    if (min(diff[diff == min(diff)]/(max(var2) - min(var2))/(max(var1) - min(var1))) < 0.01) {
      if (length(obs) == length(var1))
        obs[diff == min(diff)] <- !obs[diff == min(diff)]  
      else
        obs[diff == min(diff), ] <- !obs[diff == min(diff), ]   
    }
    return(obs)
  }
  
####################################################
  # Selection d'un polygone
####################################################

 if (method == "poly") {
   polyg <- cbind(unlist(Xpoly), unlist(Ypoly))
   pol <- st_sfc(st_polygon(list(polyg)))
   my_points <- st_as_sf(data.frame(x = var1, y = var2), coords = c("x", "y"))
   def <- as.vector(st_intersects(my_points, pol, sparse = FALSE))
   if (length(obs) == length(var1))
     obs[def] <- !obs[def]
   else
     obs[def, ] <- !obs[def, ]
   
  return(obs)
 }
  
}

