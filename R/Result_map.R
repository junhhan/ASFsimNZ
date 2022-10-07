#' Draw a figure of ASFV presence on the map
#'
#' @description This function generates a map of the probability of ASFV presence over the simulated area at a certain timepoint.
#'
#' @param res An output vector (for single core) / list (for multiple cores) of \code{\link{ASF_simulation}}.
#' @param year Year of the chosen timepoint. Input value should be between 1 and 5 which are the minimum and maximum simulation periods in year, respectively.
#' @param week Week of the chosen timepoint. Input value should be between 1 and 52.
#'
#' @details In order to generate a map, end-users should designate the timepoint that they want to inspect.
#' For example, the very first week of the simulation period can be chosen by setting both the year and week as 1.
#' @export

Result_map <- function (res, year, week) {
  run <- 1
  
  v.total_iter <- res$param[1]
  v.ncore <- res$param[2]
  v.ncell <- res$param[3]
  v.simyear <- res$param[4]
  v.pop <- res$param[5]
  v.scale <- res$param[6]
  v.island <- res$param[7]
  v.introloc <- res$param[8]
  v.control_index <- res$param[9]
  v.foi <- res$param[10]
  v.year <- year
  v.week <- week
  
  if (v.year > 5 | v.year < 1) {
    print("Error: Timepoint of interest is out of bound. Check input year again.");
    run <- 0
  }
  
  if (v.week > 52 | v.week < 1) {
    print("Error: Timepoint of interest is out of bound. Check input week again.");
    run <- 0
  }
  
  if (v.island == 0) {
    map <- map_NI
    if (v.scale == 0) {
      ras <- r_NI4km
    } else {
      ras <- r_NI10km
    }
  } else {
    map <- map_SI
    if (v.scale == 0) {
      ras <- r_SI4km
    } else {
      ras <- r_SI10km
    }
  }
  
  if (run == 1) {
    if (v.ncore == 1) {
      v <- matrix(res[[2]]$map, nrow= v.ncell, byrow= F)
    } else {
      v <- res[[2]][[1]][[1]]$map
      for (i in 2:v.ncore) {
        v <- v + res[[2]][[i]][[1]]$map
      }
      v <- matrix(v, nrow= v.ncell, byrow= F)
    }
    k <- 52 * (v.year - 1) + v.week
    tmp <- ras
    tmp@data@values <- 100 * v[,k] / v.total_iter
    tmp@data@values[is.na(ras@data@values)] <- NA
    tmp <- as.data.frame(as(tmp, "SpatialPixelsDataFrame"))
    dat_m <- data.frame("x"= tmp$x, "y"= tmp$y, "val"= tmp$layer)
    
    ggplot(data= dat_m) +
      geom_tile(aes(x= x, y= y, fill= val), alpha=0.8) +
      geom_polygon(data= map, aes(x= long, y= lat, group= group), colour= "grey", fill= NA, size= 0.1) +
      scale_fill_gradientn(
        colours= c("white", "yellow", "red"), values= rescale(c(0, 25, 100)),
        breaks= c(0, 25, 50, 75, 100),
        limits= c(0, 100), labels= c("0%", "25%", "50%", "75%", "100%")) +
      coord_equal() +
      labs(x= NULL, y= NULL, fill= NULL) +
      theme(
        axis.text= element_blank(),
        axis.ticks= element_blank(),
        panel.background= element_blank(),
        strip.text.y= element_blank())
  }
}
