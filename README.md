---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
```

# ASFsimNZ
A user guide for *ASFsimNZ*: an R package for simulating ASF spread in New Zealand feral pigs.

## Introduction
*ASFsimNZ* is an R package to simulate African swine fever (ASF) transmission in New Zealand feral pigs. Specifically, it is designed to (1) simulate the spread of ASF under different assumptions about the demography of the animals, and (2) examine the efficiency of different ASF control scenarios. The model was developed based on the contract between the Ministry for Primary Industries (MPI) and EpiCentre, Massey University (hereafter, we) for the project “Understanding the role of feral pigs in exotic disease incursion – African swine fever disease model development (Agreement number: C0033724)”.

## Model description
To simulate the spread of ASF, we developed a stochastic grid-based simulation model. In the model, New Zealand is converted into a grid with each cell representing either;  
  
* a habitat for a sounder (i.e. a family group of feral pigs), 
* an accessible territory (i.e. locations where feral pigs could trespass but are not suitable to inhabit), or
* a barrier (i.e. locations where feral pigs cannot or do not approach)
  
depending on the landscape of the cell.  
  
For each habitat, feral pigs experienced demographic events, such as breeding, farrowing, and weaning. Also, feral pigs could migrate to a different habitat when they reached to a certain age.  
  
Feral pigs could be infected with ASF if there were any infectious live pigs or carcasses in the same habitat. ASF infection could also occur by the presence of infectious pigs or carcasses in habitats nearby with the strength of transmission exponentially reducing over the distance between habitats. Moreover, ASF could spread via the migration of feral pigs.  
  
Given the uncertainty around the total population and average habitat size of feral pigs in New Zealand, the model was developed with different total populations (i.e. small, medium, and large population sizes which are equivalent to approximately 160,000, 480,000, and 960,000 of non-piglet feral pigs, respectively, across the country) and habitat scales (i.e. small and large habitat scales which are 4km^2^ and 10km^2^, respectively). Considering the geographic isolation as well as the computational burden of the simulation, the feral pigs in the North and South Islands were modelled separately. For each island, two potential ASF entry points were considered (e.g. Auckland and Wellington for the North Island, and Picton and Dunedin for the South Island). This resulted in the development of 24 sub-models (i.e. 2 islands × 3 population sizes × 2 habitat scales × 2 entry points), one of which should be selected by users to run the model. The temporal unit of the model was one week, and the transmission can be simulated **up to five years** from the introduction of ASF virus. A more detailed description of the model development and specification is written in the Final Report of the project.  
  
## Installation
### System requirement
Considering the characteristic of stochastic models, it is recommended to iterate the simulation multiple times (e.g. 100 ~ 200 times at least). In order to maximise the efficiency of model running, we therefore designed *ASFsimNZ* to use multiple CPU threads (i.e. parallel computation). Using a modern desktop/laptop, it takes up to 200 megabytes of memory (i.e. ram) with the computation time of approximately one minute for a single iteration using a single CPU thread. Therefore, a system with **> 10 CPU threads** and **> 2 gigabytes of memory** is required to run the model 200 times within an acceptable computation time (i.e. 20 minutes or so).

### Package install
*ASFsimNZ* should be manually installed on your local machine. Before the installation, we recommend using RStudio to make the installation process easier. Also, your system may require 'Rtools' to install the package. For Windows users, Rtools can be downloaded from https://cran.r-project.org/bin/windows/Rtools/.  
Once your system is ready, the package can be installed as below;

1. Open RStudio. On RStudio, choose "Package" tab on the mid-right section of the screen.
2. Click "Install" button. It generates a pop-up menu "Install packages".
3. On the pop-up menu, select "Package Archive File (.zip; .tar.gz)" option for "Install from:". This generates another pop-up screen to browse the package archive file. 
4. Browse the package archive file. The file name should be like "ASFsimNZ_1.1.4.tar.gz". Please note that the version of the package (e.g. 1.1.4) would vary depending on the latest updates.
5. If the pop-up screen for browsing in step 3 is not appearing, you can manually navigate the file with "Browse..." button.
6. Click the "Install" button.  

### Loading the package
In order to run the simulation model and execute the function to export the model result, *ASFsimNZ* requires other packages below.
```{r, echo= T, eval= T, message= F}
library(parallel)
library(rgdal)
library(sp)
library(raster)
library(ggplot2)
library(scales)
```
If missing, any packages can be installed with the code below.
```{r, echo= T, eval= F}
install.packages("rgdal") # Install "rgdal" package
```
Please make sure that all the packages above are installed and loaded on your system. Once the packages are loaded, your system is ready to load *ASFsimNZ*.
```{r, echo= T, eval= F}
library(ASFsimNZ)
```

## Running the model
As mentioned above, users should select the model assumptions and settings. For the purpose of demonstration, we present the case of simulating ASF spread in the North Island when the medium population size and small habitat scale are assumed. Also, ASF is assumed to be introduced via Auckland. The model is simulated for five years from the introduction, and is iterated for 200 times.  
``` {r, echo= T, eval= F}
total_iter <- 200 # Number of iterations
simyear <- 5 # Simulation period in years
pop <- 1 # Indicator for the population size (0: small, 1: medium, 2: high)
scale <- 0 # Indicator for the habitat scale (0: small, 1: large)
island <- 0 # Indicator for the Island (0: North, 1: South)
introloc <- 0 # Indicator for ASF introduction (0: Auckland/Picton, 1: Wellington/Dunedin)
```
Users should designate the number of CPU threads for model running. In case someone is not sure how many CPU threads his/her system can use, the number of threads can be checked with the code below.
```{r, echo= T, eval= F}
detectCores()
```
Now, the system is ready to run the model.  
  
The model can be run with the code below.
```{r, ehco= T, eval= F}
res <- ASF_simulation(total_iter= total_iter, ncore= ncore, simyear= simyear, 
                      pop= pop, scale= scale, island= island, introloc= introloc)
```
The result of the simulation is stored in the vector "res".  
  
## Model results
The model result can be inspected in three ways.  
  
1. Time-series of total habitat areas that carrying ASF virus over the simulation period.
2. Probability of habitats carrying ASF virus at a certain time-point. 
3. List of farms located near the habitats carrying ASF virus at a certain time-point.  
  
### Time-series of affected area
The time-series is presented with the code below.
```{r, echo= T, eval= F}
Result_ts(res)
```
  
The code generates a time-series plot of the total habitat areas that carrying ASF virus. The red line (grey area) indicates the median (interquartile range) of the predicted area.  
  
Users can also examine the data that the plot is based on.
```{r, echo= T, eval= F}
v <- Result_ts(res)
head(v$data)
```
  
In the data, column "t" indicates the simulation week, and "val", "low", "high" indicate the median, 1^st^, 3^rd^ quartile of the affected habitat area. 

### Spatial extent of ASF spread
The spatial extent of ASF spread can be visualised. Given the spatial scale of ASF spread is time-varying, users should specify the timepoint for the visualisation. For example, if someone wants to examine 26 weeks after the introduction of ASF virus, the timepoint should be labelled as "year 1" and "week 26" in the code below.
```{r, echo= T, eval= F, warning= F}
Result_map(res, 1, 26) # 26 weeks after the introduction
```
  
The value of the map is the probability that a habitat carrying at least one infectious feral pig or carcass at the chosen week. If someone wants to examine 120 weeks after the introduction, the timepoint is then "year 3" and "week 16" (120 = 52 * 2 + 16, so week 16 of year 3).  
```{r, echo= T, eval= F, warning= F}
Result_map(res, 3, 16) # 120 weeks after the introduction
```
  
## ASF control
### Description
In *ASFsimNZ*, four control scenarios can be examined. Each control scenario is a combination of different control measures (e.g. increasing passive surveillance, increasing carcass removal, increasing hunting pressure) on different control zones, including the infected, buffer, and surveillance zones. The infected zone is an area that immediately surrounds a habitat cell where an ASF infectious carcass is found via passive surveillance. The buffer and surveillance zones are the areas surrounding an infected and buffer zone, respectively. The radii for the infected, buffer, and surveillance zones are 4km, 8km, and 15km, respectively.  
  
### Control scenarios
The control scenarios are;  
  
1. Scenario 1 (baseline scenario): In this scenario, there is an increased level of passive surveillance, carcass removal, and hunting pressure inside the infected zones. Also, increased passive surveillance is implemented in the buffer zones. No control measures are conducted inside the surveillance zones.
2. Scenario 2: Compared with Scenario 1, the radii of infected and buffer zones are expanded to 20km and 28km, respectively.
3. Scenario 3: On top of Scenario 1, carcass removal is added to the buffer zone.
4. Scenario 4: Over the baseline scenario, increased passive surveillance is implemented in the surveillance zone.  
  
### Efficiency of scenarios
For the demonstration, we choose Scenario 1. To examine the efficiency, users should run the model while indicating the scenario type. If "control_index" is missing or its value is 0, no control scenario is implemented in the simulation.
```{r, ehco= T, eval= F}
control_index <- 1 # Indicator for control scenario (1, 2, 3, or 4)
cnt <- ASF_simulation(total_iter= total_iter, ncore= ncore, simyear= simyear, 
                      pop= pop, scale= scale, island= island, introloc= introloc, 
                      control_index= control_index)
```
With the scenario, the spatial and temporal extent of ASF spread can be examined and compared.
```{r, echo= T, eval= F}
Result_ts(cnt)
Result_map(cnt, 3, 16)
```
  
## Additional feature
One of the most important parts of disease simulation modelling is how to mathematically formulate the transmission, which is often calculated using the force of infection (FOI). In **'ASF_simulation()'**, the transmission of ASF in susceptible animals is based on the method used in Han et al., (2021). However, to compare the way of calculating FOI on the scale of ASF spread, we also implemented different FOI calculation methods in **'ASF_simulation()'**. Specifically, two additional methods are considered;  

1. Based on Pepin et al., (2022)
2. Based on Lange et al., (2018)  
  
```{r, ehco= T, eval= F}
foi <- 1 # Indicator for FOI calculation 1 (1: Pepin(2022), 2: Lange(2018))
pm <- ASF_simulation(total_iter= total_iter, ncore= ncore, simyear= simyear, 
                     pop= pop, scale= scale, island= island, introloc= introloc, 
                     foi= foi)
```
 
