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
MN_env_stack <- stack("/circuitscape/data/MN_centric_environmental_rasters/MN_env_run_varbs.grd")

temp_wq <- raster::subset(MN_env_stack, 1)
precip_wq <- raster::subset(MN_env_stack, 3)
temp_cm <- raster::subset(MN_env_stack, 4)
lc <- raster::subset(MN_env_stack, 6)
imp <- raster::subset(MN_env_stack, 7)
drain <- stack("/circuitscape/data/MN_centric_environmental_rasters/drainage_bb.grd")
soilcomp1 <- stack("/circuitscape/data/MN_centric_environmental_rasters/soil_comp_pca1_transpose_bb.grd")
soilcomp1 <- resample(soilcomp1, drain)
soilcomp2 <- stack("/circuitscape/data/MN_centric_environmental_rasters/soil_comp_pca2_transpose_bb.grd")
soilcomp2 <- resample(soilcomp2, drain)
soilchem1 <- stack("/circuitscape/data/MN_centric_environmental_rasters/soil_min_pca1_transpose_bb.grd")
soilchem1 <- resample(soilchem1, drain)
soilchem2 <- stack("/circuitscape/data/MN_centric_environmental_rasters/soil_min_pca2_transpose_bb.grd")
soilchem2 <- resample(soilchem2, drain)

MN_env_stack_1 <- stack(temp_wq, precip_wq, temp_cm, lc, imp, soilcomp1, soilcomp2, 
soilchem1, soilchem2, drain)

paste("Connecting R to Julia")
Sys.time()
JULIA_HOME <- "/.conda/envs/ResistanceGA-env/bin"
JuliaCall::julia_setup(JULIA_HOME)

paste("Preparing inputs for Resistance GA")
Sys.time()
GA.inputs <- GA.prep(ASCII.dir = MN_env_stack_1,
                     Results.dir = "/circuitscape/",
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

jl.ms.optim <- MS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

paste("Saving objects")
save(GA.inputs, jl.inputs, jl.ms.optim, file = "Multisurface_circuitscape_results.rda")

Sys.time()

paste("Fin")
                   
Sys.time()


