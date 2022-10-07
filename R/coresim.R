#' Core ASF simulation function
#'
#' @param input_info (internally generated) contains information of cells.
#' @param input_road (internally generated) contains the presence of roads between neighbouring cells.
#' @param input_river (internally generated) contains the presence of rivers between neighbouring cells.
#' @param input_proxycid (internally generated) contains the cell id of proximal cells.
#' @param input_proxydist (internally generated) contains the distance to proximal cells.
#' @param input_proxyroad (internally generated) contains the presence of roads between proximal cells.
#' @param input_proxyriver (internally generated) contains the presence of rivers between proximal cells.
#' @param coreid (internally generated) indicates CPU ID for parallel computation.
#' @param rep (internally generated) indicates the number of iterations for each core.
#' @param simyear Simulation period in year. Input value should be an integer between one and five.
#' @param pop Total population size of feral pig. Input value should be either 0 (= small), 1 (= intermediate), or 2 (= large).
#' @param scale Average habitat size of feral pigs. Input value should be either 0 (= 4km2) or 1 (= 10km2).
#' @param island Input value should be either 0 (= North Island) or 1 (= South Island).
#' @param introloc ASF entry point. Input value should be either 0 (= Auckland or Picton for the North and South Island, respectively) or 1 (= Wellington or Dunedin for the North and South Island, respectively).
#' @param control_index Control scenarios. Input value should be either 0, 1, 2, 3, or 4. (0: No control. 1: Baseline control. 2: Increase control range. 3: Add carcass removal in buffer zone. 4: Add passive surveillance in surveillance zone).
#' @param tmp_map (internally generated) is the temporary file to store the spatial result.
#' @param tmp_ts (internally generated) is the temporary file to store the temporal result.
#' @param foi (optional) A way to calculate force of infection. Input value should be either 0, 1, or 2. (0: Han et al (2021). 1: Pepin et al (2022). 2: Thulke et al (2018). Default method is Han et al (2021).
#'
#' @return This function returns the number of iterations that each cell having ASFV over the simulation period.
#' @export
#' @useDynLib ASFsimNZ
#'

coresim <- function (
    input_info, input_road, input_river, input_proxycid, input_proxydist, input_proxyroad, input_proxyriver,
    coreid, rep, simyear, pop, scale, island, introloc, control_index, foi, tmp_map, tmp_ts) {
  
  v.rep <- rep
  v.simyear <- simyear
  v.pop <- pop
  v.scale <- scale
  v.island <- island
  v.introloc <- introloc
  v.control_index <- control_index
  v.foi <- foi
  
  .C("_ASF_simulation",
     as.integer(input_info), as.integer(input_road), as.integer(input_river),
     as.integer(input_proxycid), as.double(input_proxydist),
     as.integer(input_proxyroad), as.integer(input_proxyriver),
     as.integer(coreid), as.integer(v.rep), as.integer(v.simyear), as.integer(v.pop), as.integer(v.scale),
     as.integer(v.island), as.integer(v.introloc), as.integer(v.control_index), as.integer(v.foi),
     as.integer(tmp_map), as.integer(tmp_ts))
}
