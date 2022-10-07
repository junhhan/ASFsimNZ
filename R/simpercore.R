#' Allocating ASF simulations to each core
#'
#' @param rep Number of iterations per each core.
#' @param simyear Simulation period in year. Input value should be an integer between one and five.
#' @param pop Total population size of feral pig. Input value should be either 0 (= small), 1 (= intermediate), or 2 (= large).
#' @param scale Average habitat size of feral pigs. Input value should be either 0 (= 4km2) or 1 (= 10km2).
#' @param island Input value should be either 0 (= North Island) or 1 (= South Island).
#' @param introloc ASF entry point. Input value should be either 0 (= Auckland or Picton for the North and South Island, respectively) or 1 (= Wellington or Dunedin for the North and South Island, respectively).
#' @param control_index Control scenarios. Input value should be either 0, 1, 2, 3, or 4. (0: No control. 1: Baseline control. 2: Increase control range. 3: Add carcass removal in buffer zone. 4: Add passive surveillance in surveillance zone).
#' @param foi (optional) A way to calculate force of infection. Input value should be either 0, 1, or 2. (0: Han et al (2021). 1: Pepin et al (2022). 2: Thulke et al (2018). Default method is Han et al (2021).
#'
#' @return This function allocates and runs the simulation model to each core.
#' @export

simpercore <- function(rep, simyear, pop, scale, island, introloc, control_index, foi) {
  
  v.rep <- rep
  v.simyear <- simyear
  v.pop <- pop
  v.scale <- scale
  v.island <- island
  v.introloc <- introloc
  v.control_index <- control_index
  v.foi <- foi
  
  # Prepare datasets based on the input values
  if (v.scale == 0) {
    if (v.island == 0) {
      v.ncell <- 105324
      input_info <- as.numeric(as.matrix(NI4km))
      input_road <- as.numeric(as.matrix(road_NI4km))
      input_river <- as.numeric(as.matrix(river_NI4km))
      input_proxycid <- as.numeric(as.matrix(proxycid_NI4km))
      input_proxydist <- as.numeric(as.matrix(proxydist_NI4km))
      input_proxyroad <- as.numeric(as.matrix(proxyroad_NI4km))
      input_proxyriver <- as.numeric(as.matrix(proxyriver_NI4km))
    } else {
      v.ncell <- 109652
      input_info <- as.numeric(as.matrix(SI4km))
      input_road <- as.numeric(as.matrix(road_SI4km))
      input_river <- as.numeric(as.matrix(river_SI4km))
      input_proxycid <- as.numeric(as.matrix(proxycid_SI4km))
      input_proxydist <- as.numeric(as.matrix(proxydist_SI4km))
      input_proxyroad <- as.numeric(as.matrix(proxyroad_SI4km))
      input_proxyriver <- as.numeric(as.matrix(proxyriver_SI4km))
    }
  } else {
    if (v.island == 0) {
      v.ncell <- 42330
      input_info <- as.numeric(as.matrix(NI10km))
      input_road <- as.numeric(as.matrix(road_NI10km))
      input_river <- as.numeric(as.matrix(river_NI10km))
      input_proxycid <- as.numeric(as.matrix(proxycid_NI10km))
      input_proxydist <- as.numeric(as.matrix(proxydist_NI10km))
      input_proxyroad <- as.numeric(as.matrix(proxyroad_NI10km))
      input_proxyriver <- as.numeric(as.matrix(proxyriver_NI10km))
    } else {
      v.ncell <- 44000
      input_info <- as.numeric(as.matrix(SI10km))
      input_road <- as.numeric(as.matrix(road_SI10km))
      input_river <- as.numeric(as.matrix(river_SI10km))
      input_proxycid <- as.numeric(as.matrix(proxycid_SI10km))
      input_proxydist <- as.numeric(as.matrix(proxydist_SI10km))
      input_proxyroad <- as.numeric(as.matrix(proxyroad_SI10km))
      input_proxyriver <- as.numeric(as.matrix(proxyriver_SI10km))
    }
  }
  tmp_map <- rep(0, v.ncell * v.simyear * 52)
  tmp_ts <- rep(0, v.rep * v.simyear * 52)
  
  # Run the simulation model
  if (!exists("coreid")) {
    coreid <- 0
  }
  val <- coresim(
    input_info, input_road, input_river, input_proxycid, input_proxydist, input_proxyroad, input_proxyriver,
    coreid, v.rep, v.simyear, v.pop, v.scale, v.island, v.introloc, v.control_index, v.foi, tmp_map, tmp_ts)
  l <- list("map"= val[[17]], "ts"= val[[18]])
  gc()
  
  return(l);
}
