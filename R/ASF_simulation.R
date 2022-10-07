#' Simulation of ASF transmission in New Zealand feral pigs
#'
#' @param total_iter Number of total iterations. Recommend either 100 or 200.
#' @param ncore Number of CPU threads available for parallel computation.
#' @param simyear Simulation period in year. Input value should be an integer between one and five.
#' @param pop Total population size of feral pig. Input value should be either;
#' \itemize{
#'    \item 0 (= small),
#'    \item 1 (= intermediate),
#'    \item 2 (= large).
#' }
#' @param scale Average habitat size of feral pigs. Input value should be either;
#' \itemize{
#'    \item 0 (= 4km2),
#'    \item 1 (= 10km2).
#' }
#' @param island Input value should be either;
#' \itemize{
#'    \item 0 (= North Island),
#'    \item 1 (= South Island).
#' }
#' @param introloc ASF entry point. Input value should be either;
#' \itemize{
#'    \item 0 (= Auckland or Picton for the North and South Island, respectively),
#'    \item 1 (= Wellington or Dunedin for the North and South Island, respectively).
#' }
#' @param control_index Control scenarios. Input value should be either;
#' \itemize{
#'    \item 0 (= no control),
#'    \item 1 (= Baseline control),
#'    \item 2 (= Increase control range),
#'    \item 3 (= Add carcass removal in buffer zone),
#'    \item 4 (= Add passive surveillance in surveillance zone).
#' }
#' Default scenario is no control.
#' @param foi (optional) A way to calculate force of infection. Input value should be either;
#' \itemize{
#'    \item 0 (= Han et al (2021)),
#'    \item 1 (= Pepin et al (2022)),
#'    \item 2 (= Thulke et al (2018)).
#' }
#' Default method is Han et al (2021).
#'
#' @return This function returns;
#' \itemize{
#'    \item the number of iterations having ASFV for each cell throughout the simulation period, and
#'    \item the number of habitat cells having ASFV for each week of the simulation period.
#'}
#' @export
#'
#' @examples ASF_simulation(total_iter= 200, ncore= 10, simyear= 5, pop= 1, scale= 0, island= 0, introloc= 0, control_index= 0)

ASF_simulation <- function (total_iter, ncore, simyear, pop, scale, island, introloc, control_index, foi) {
  run <- 0 # An index for validation of input values
  if (missing(control_index)) {
    control_index <- 0;
  }
  if (missing(foi)) {
    foi <- 0;
  }
  
  #if(!exists("asn.env")) {
  #  asn.env <<- new.env()
  #}
  
  v.total_iter <- total_iter
  v.ncore <- ncore
  v.simyear <- simyear
  v.pop <- pop
  v.scale <- scale
  v.island <- island
  v.introloc <- introloc
  v.control_index <- control_index
  v.foi <- foi
  
  while (run == 0) {
    
    if (v.ncore < 1) {
      print("Error: Number of core should be >= 1. Using at least 4 cores is recommended.");
      break;
    }
    if (v.ncore == 1) {
      print("Warning: Single core is used. Running time can be ~3 hours if the model is iterated for 200 times.")
    }
    if (v.total_iter %% v.ncore != 0) {
      print("Warning: Simulation cannot be evenly distributed to each core. Please adjust either the number of iterations or cores.")
      break;
    }
    
    if (v.simyear > 5) {
      print("Error: Maximum simulation period is five years.");
      break;
    } else if (v.simyear < 1) {
      print("Error: Simulation period cannot be negative.");
      break;
    }
    
    if (v.pop > 2 | v.pop < 0) {
      print("Error: Total population should be either low (= 0), intermediate (= 1), or high (= 2).");
      break;
    }
    
    if (v.scale > 1 | v.scale < 0) {
      print("Error: Habitat size should be either small (= 0), or large (= 1).");
      break;
    }
    
    if (v.island > 1 | v.island < 0) {
      print("Error: Selected island should be either the North (= 0) or South (= 1) Island.");
      break;
    }
    
    if (v.introloc > 1 | v.introloc < 0) {
      print("Error: ASF should be introduced via Auckland (= 0) or Wellington (=1) for the North Island,");
      print("or Picton (= 0) or Dunedin (= 1) for the South Island.");
      break;
    }
    
    if (v.control_index > 4 | v.control_index < 0) {
      print("Error: Please choose the proper control scenario:")
      print(" No control (= 0)");
      print(" Baseline control (= 1)");
      print(" Increase control range (= 2)");
      print(" Add carcass removal in buffer zone (= 3)");
      print(" Add passive surveillance in surveillance zone (= 4)");
      break;
    }
    
    if (v.foi > 3 | v.foi < 0) {
      print("Error: Please choose the proper way of calculating the force of infection:")
      print(" Han et al (2021) (= 0)");
      print(" Pepin et al (2022) (= 1)");
      print(" Lange et al (2018) (= 2)");
      print("Han et al (2021) is selected as a default.")
      v.foi <- 0
    }
    
    run <- run + 1;
  }
  
  if (run == 1) {
    if (v.scale == 0) {
      if (v.island == 0) {
        v.ncell <- 105324
      } else {
        v.ncell <- 109652
      }
    } else {
      if (v.island == 0) {
        v.ncell <- 42330
      } else {
        v.ncell <- 44000
      }
    }
    #asn.env$ncell <- v.ncell
    #assign("ncell", v.ncell, envir= .GlobalEnv)
    v.rep <- v.total_iter / v.ncore;
    #asn.env$rep <- v.rep
    
    if (v.ncore == 1) {
      res <- simpercore(v.rep, v.simyear, v.pop, v.scale, v.island, v.introloc, v.control_index, v.foi)
      res <- list("param"= c(v.total_iter, v.ncore, v.ncell, v.simyear,
                             v.pop, v.scale, v.island, v.introloc, v.control_index, v.foi), res)
    } else {
      cl <- makeCluster(v.ncore)
      clusterExport(cl, c('v.rep', 'v.simyear', 'v.pop', 'v.scale', 'v.island', 'v.introloc', 'v.control_index', 'v.foi'),
                    envir= environment())
      clusterApply(cl, seq_along(cl), function(i) coreid <<- i - 1)
      res <-
        clusterEvalQ(cl, {
          library(ASFsimNZ)
          l <- simpercore(v.rep, v.simyear, v.pop, v.scale, v.island, v.introloc, v.control_index, v.foi)
          list(l)
        })
      res <- list("param"= c(v.total_iter, v.ncore, v.ncell, v.simyear,
                             v.pop, v.scale, v.island, v.introloc, v.control_index, v.foi), res)
      stopCluster(cl)
    }
    gc()
    return(res);
  }
}
