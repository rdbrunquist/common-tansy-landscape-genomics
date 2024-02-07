#!/usr/bin/env Rscript

# This script is a sample of the scripts run on a University cluster - prior to running the script 
# Julia and the CircuitScape models need to be set-up and the link to Julia through R needs to be 
# established - See ResistanceGA documentation for more information.

paste("Loading Libraries")
#Libraries
library(ResistanceGA)

paste("Reading in coordinate data")
#MN tansy coordinate
xy_MN <- read.csv("/circuitscape/data/xy_MN.csv", header = TRUE)
coordinates(xy_MN) <- ~ lon + lat
proj4string(xy_MN) <- CRS("+proj=longlat +ellps=WGS84 +no_defs")

paste("Reading in genetic distance matrix")
# Tansy genetic Chord Distance
gen_mat <- read.csv("/circuitscape/data/MN_tansy_gen_chord_dist_matrix.csv", header = TRUE)

paste("Reading in environmental data stack")
#Environmental raster stack
MN_env_stack <- stack("/circuitscape/data/MN_centric_environmental_rasters/MN_env_varbs.grd")

paste("Isolating layer to be analyzed: Mean Temp of the Warmest Quarter")
warm_q <- raster::subset(MN_env_stack, 2)


paste("Connecting R to Julia")
Sys.time()
JULIA_HOME <- "/.conda/envs/ResistanceGA-env/bin"
JuliaCall::julia_setup(JULIA_HOME)

paste("Preparing inputs for Resistance GA")
Sys.time()
GA.inputs <- GA.prep(ASCII.dir = warm_q,
                     Results.dir = "/circuitscape/",
                     select.trans = list("A"),
                     method = "LL",
                     parallel = 4)

paste("Preparing data to be included in Julia Circuitscape")
Sys.time()
jl.inputs <- jl.prep(n.Pops = length(xy_MN),
                     response = lower(gen_mat),
                     CS_Point.File = xy_MN,
                     cholmod = TRUE,
                     JULIA_HOME = JULIA_HOME)
                     
paste("Running Circuitscape")
Sys.time()
jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

paste("Saving objects")
save(GA.inputs, jl.inputs, jl.optim, file = "mean_temp_warm_q_circuitscape_results.rda")

Sys.time()

paste("Fin")
                   
Sys.time()


