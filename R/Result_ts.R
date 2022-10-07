#' Draw a time-series of ASFV presence
#'
#' @description This function generates a time-series of the total area that ASFV present over the simulated area.
#'
#' @param res An output vector (for single core) / list (for multiple cores) of \code{\link{ASF_simulation}}.
#'
#' @return A figure of time-series of the total area that ASF virus present. The red line (grey area) indicates the median (interquartile) of the affected area.
#' \itemize{
#'   \item val Median.
#'   \item low 25th percentile.
#'   \item high 75th percentile.
#' }
#' @export

Result_ts <- function (res) {
  run <- 1
  
  if (!exists("asn.env")) {
    print("Error: Environment storing the parameters is missing.")
    run <- 0
  } else {
    v.simyear <- asn.env$simyear
    v.ncore <- asn.env$ncore
    v.ncell <- asn.env$ncell
    v.total_iter <- asn.env$total_iter
    v.scale <- asn.env$scale
  }
  
  if (run == 1) {
    if (v.ncore == 1) {
      v <- matrix(res$ts, nrow= v.total_iter / v.ncore, byrow= F)
    } else {
      v <- matrix(res[[1]][[1]]$ts, nrow= v.total_iter / v.ncore, byrow= F)
      for (i in 2:v.ncore) {
        v <- rbind(v, matrix(res[[i]][[1]]$ts, nrow= v.total_iter / v.ncore, byrow= F))
      }
    }
    
    simperiod <- v.simyear * 52
    if (v.scale == 0) {
      hs <- 4
    } else {
      hs <- 10
    }
    
    dat_g <- data.frame("t"= 0, "val"= 0, "low"= 0, "high"= 0)
    for (i in c(1:simperiod)) {
      dat_g <- rbind(dat_g, rbind(data.frame("t"= i, "val"= hs * quantile(v[,i], 0.5),
                                             "low"= hs * quantile(v[,i], 0.25),
                                             "high"= hs * quantile(v[,i], 0.75))))
    }
    
    yl <- 1.2 * max(dat_g$high)
    digit <- floor(log10(yl)) + 1
    if (digit == 1) {
      ylim <- 10
    } else {
      ylim <- round(yl, -(digit - 2))
    }
    brks <- seq(0, ylim, length.out= 4)
    labs <- paste(round(seq(0, ylim, length.out= 4),0))
    print(dat_g)
    
    ggplot(data= dat_g, aes(x= t, y= val)) +
      geom_line(colour= "red") +
      geom_ribbon(aes(x= t, ymin= low, ymax= high), fill= "grey", alpha= 0.2) +
      scale_y_continuous(expand= c(0,0), limits= c(0, ylim), breaks= brks, labels= labs) +
      scale_x_continuous(expand= c(0,0), breaks= 52 * (c(1:simyear) - 1) + 26.5,
                         labels= paste("Year", c(1:simyear), sep= " ")) +
      labs(x= "Simulation year", y= bquote('Areas with ASF virus ' ~(km^2)), colour= "black") +
      theme(
        axis.line= element_line(colour= "black"),
        legend.position= "none",
        panel.grid.major.y= element_line(colour= "gray"),
        panel.grid.minor.x= element_line(colour= "gray"),
        panel.background= element_blank())
  }
}
