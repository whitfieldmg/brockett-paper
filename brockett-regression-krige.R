library(rgdal)
library(geoR)
library(gstat)
library(sf)
library(dplyr)

## Function from Stackoverflow to remove rows containing missing data from spatial objects
# FUNCTION TO REMOVE NA's IN sp DataFrame OBJECT
#   x           sp spatial DataFrame object
#   margin      Remove rows (1) or columns (2) 
sp.na.omit <- function(x, margin=1) {
  if (!inherits(x, "SpatialPointsDataFrame") & !inherits(x, "SpatialPolygonsDataFrame")) 
    stop("MUST BE sp SpatialPointsDataFrame OR SpatialPolygonsDataFrame CLASS OBJECT") 
  na.index <- unique(as.data.frame(which(is.na(x@data),arr.ind=TRUE))[,margin])
  if(margin == 1) {  
    cat("DELETING ROWS: ", na.index, "\n") 
    return( x[-na.index,]  ) 
  }
  if(margin == 2) {  
    cat("DELETING COLUMNS: ", na.index, "\n") 
    return( x[,-na.index]  ) 
  }
}

## extract_krige_stats
# Extracts summary statistics from kriging cross-validation
extract_krige_stats <- function(krige_model) {
  stat_table <- data.frame(me = NA,
                           mspe = NA,
                           msne = NA,
                           cor_obspred = NA,
                           cor_predres = NA,
                           rmse = NA)
  stat_table[1, "me"] <- mean(krige_model$residual)
  stat_table[1, "mspe"] <- mean(krige_model$residual^2)
  stat_table[1, "msne"] <- mean(krige_model$zscore^2)
  stat_table[1, "cor_obspred"] <- cor(krige_model$observed, krige_model$var1.pred)
  stat_table[1, "cor_predres"] <- cor(krige_model$var1.pred, krige_model$residual)
  stat_table[1, "rmse"] <- sqrt(mean((krige_model$observed - krige_model$var1.pred)^2))
  return(stat_table)
}

# Set paths to geodatabases
birkhowe.gdb <- "/Users/mikewhitfield/Brockett-paper/MikeBirkhowe.gdb"
hollins.gdb <- "/Users/mikewhitfield/Brockett-paper/MikeHollins.gdb"
lowsnab.gdb <- "/Users/mikewhitfield/Brockett-paper/MikeLowsnab.gdb"

# Examine Feature Classes present in geodatabases
birkhowe_fclist <- ogrListLayers(birkhowe.gdb)
birkhowe_fclist

hollins_fclist <- ogrListLayers(hollins.gdb)
hollins_fclist

lowsnab_fclist <- ogrListLayers(lowsnab.gdb)
lowsnab_fclist

# Read Feature Classes from geodatabases
# One per depth, three depths per site
birkhowe_d1 <- readOGR(dsn = birkhowe.gdb, layer = "birkhowe_soild1_pt_v5", pointDropZ = TRUE)
birkhowe_d2 <- readOGR(dsn = birkhowe.gdb, layer = "birkhowe_soild2_pt_v5", pointDropZ = TRUE)
birkhowe_d3 <- readOGR(dsn = birkhowe.gdb, layer = "birkhowe_soild3_pt_v5", pointDropZ = TRUE)
# birkhowe_man <- readOGR(dsn = birkhowe.gdb, layer = "birkhowe_managementareas_po_2")
hollins_d1 <- readOGR(dsn = hollins.gdb, layer = "hollins_soild1_pt_v3", pointDropZ = TRUE)
hollins_d2 <- readOGR(dsn = hollins.gdb, layer = "hollins_soild2_pt_v3", pointDropZ = TRUE)
hollins_d3 <- readOGR(dsn = hollins.gdb, layer = "hollins_soild3_pt_v3", pointDropZ = TRUE)
lowsnab_d1 <- readOGR(dsn = lowsnab.gdb, layer = "lowsnab_soild1_pt_v5", pointDropZ = TRUE)
lowsnab_d2 <- readOGR(dsn = lowsnab.gdb, layer = "lowsnab_soild2_pt_v5", pointDropZ = TRUE)
lowsnab_d3 <- readOGR(dsn = lowsnab.gdb, layer = "lowsnab_soild3_pt_v5", pointDropZ = TRUE)

# Read elevation data for Hollins (missing from geodatabase)
hollins_elevation <- read.csv("/Users/mikewhitfield/brockett-soilc/hollins_elevation.csv")

# Remove empty elevation column from spatial data
hollins_d1@data <- select(hollins_d1@data, -elevation)
hollins_d2@data <- select(hollins_d2@data, -elevation)
hollins_d3@data <- select(hollins_d3@data, -elevation)

hollins_d1@data <- left_join(hollins_d1@data, hollins_elevation, by = "gps_id")
hollins_d2@data <- left_join(hollins_d2@data, hollins_elevation, by = "gps_id")
hollins_d3@data <- left_join(hollins_d3@data, hollins_elevation, by = "gps_id")

# Subset necessary columns from data
sub_cols <- c("hls_plot", "elevation", "moisture", "totC_mean_mass_vol")
birkhowe_d1_sub <- birkhowe_d1[, sub_cols]
birkhowe_d2_sub <- birkhowe_d2[, sub_cols]
birkhowe_d3_sub <- birkhowe_d3[, sub_cols]

hollins_d1_sub <- hollins_d1[, sub_cols]
hollins_d2_sub <- hollins_d2[, sub_cols]
hollins_d3_sub <- hollins_d3[, sub_cols]

lowsnab_d1_sub <- lowsnab_d1[, sub_cols]
lowsnab_d2_sub <- lowsnab_d2[, sub_cols]
lowsnab_d3_sub <- lowsnab_d3[, sub_cols]

# Subset using dplyr select for sf. Currently doesn't work.
# birkhowe_d1_sub2 <- birkhowe_d1 %>%
#                       select.sf(hls_plot, elevation, moisture, totC_mean_mass_vol)

# Filter rows containing NAs
# Overwrite _sub object
birkhowe_d1_sub <- sp.na.omit(birkhowe_d1_sub)
birkhowe_d2_sub <- sp.na.omit(birkhowe_d2_sub)
birkhowe_d3_sub <- sp.na.omit(birkhowe_d3_sub)
hollins_d1_sub <- sp.na.omit(hollins_d1_sub)
hollins_d2_sub <- sp.na.omit(hollins_d2_sub)
hollins_d3_sub <- sp.na.omit(hollins_d3_sub)
# lowsnab_d1_sub <- sp.na.omit(lowsnab_d1_sub) # No empty cells. Function returns empty data
lowsnab_d2_sub <- sp.na.omit(lowsnab_d2_sub)
lowsnab_d3_sub <- sp.na.omit(lowsnab_d3_sub)

# Interactive variogram fitting procedure
# Run these lines once
# Save variogram model as geoR vgm object to disk

## Birkhowe
# birkhowe_d1.ifit <- eyefit(variog(as.geodata(birkhowe_d1["totC_mean_mass_vol"])))

# cov.model          sigmasq phi            tausq kappa kappa2   practicalRange
# 1 exponential 614.433374732558 800 76.8041718415697  <NA>   <NA> 2396.58581878812

# birkhowe_d1.vgm <- as.vgm.variomodel(birkhowe_d1.ifit[[1]])

# save(birkhowe_d1.vgm, file = "/Users/mikewhitfield/Brockett-paper/birkhowe_d1_vgm.RData")

# birkhowe_d2.ifit <- eyefit(variog(as.geodata(birkhowe_d2_sub["totC_mean_mass_vol"])))
# 
# # cov.model sigmasq    phi  tausq kappa kappa2   practicalRange
# # 1 exponential     385 2147.5 123.97  <NA>   <NA> 6433.33505730783
# 
# birkhowe_d2.vgm <- as.vgm.variomodel(birkhowe_d2.ifit[[1]])
# 
# save(birkhowe_d2.vgm, file = "/Users/mikewhitfield/Brockett-paper/birkhowe_d2_vgm.Rdata")

# birkhowe_d3.ifit <- eyefit(variog(as.geodata(birkhowe_d3_sub["totC_mean_mass_vol"])))
# 
# # cov.model sigmasq phi            tausq kappa kappa2   practicalRange
# # 1 exponential    5250 200 556.361701560112  <NA>   <NA> 599.146454697692
# 
# birkhowe_d3.vgm <- as.vgm.variomodel(birkhowe_d3.ifit[[1]])
# 
# save(birkhowe_d3.vgm, file = "/Users/mikewhitfield/Brockett-paper/birkhowe_d3_vgm.Rdata")

# Kriging with prediction. Doesn't work without a grid to predict over.
# birkhowe_d1.rk <- krige(formula = totC_mean_mass_vol ~ elevation + moisture + as.factor(hls_plot),
#                         birkhowe_d1,
#                         model = birkhowe_d1.vgm)

## Lowsnab
# lowsnab_d1.ifit <- eyefit(variog(as.geodata(lowsnab_d1_sub["totC_mean_mass_vol"])))
# 
# # cov.model sigmasq    phi            tausq kappa kappa2   practicalRange
# # 1 exponential  772.08 318.04 105.803498263631  <NA>   <NA> 952.762692259754
# # 2 exponential  772.08 318.04 105.803498263631  <NA>   <NA> 952.762692259754
# 
# lowsnab_d1.vgm <- as.vgm.variomodel(lowsnab_d1.ifit[[1]])
# 
# save(lowsnab_d1.vgm, file = "/Users/mikewhitfield/Brockett-paper/lowsnab_d1_vgm.Rdata")

# lowsnab_d2.ifit <- eyefit(variog(as.geodata(lowsnab_d2_sub["totC_mean_mass_vol"])))
# 
# # cov.model sigmasq   phi            tausq kappa kappa2   practicalRange
# # 1 exponential 1965.53 53.01 269.350236476898  <NA>   <NA> 158.803767818213
# 
# lowsnab_d2.vgm <- as.vgm.variomodel(lowsnab_d2.ifit[[1]])
# 
# save(lowsnab_d2.vgm, file = "/Users/mikewhitfield/Brockett-paper/lowsnab_d2_vgm.Rdata")

lowsnab_d3.ifit <- eyefit(variog(as.geodata(lowsnab_d3_sub["totC_mean_mass_vol"])))

# cov.model sigmasq   phi            tausq kappa kappa2   practicalRange
# 1 exponential    2100 53.01 295.711917421923  <NA>   <NA> 158.803767818213

lowsnab_d3.vgm <- as.vgm.variomodel(lowsnab_d3.ifit[[1]])

save(lowsnab_d3.vgm, file = "/Users/mikewhitfield/Brockett-paper/lowsnab_d3_vgm.Rdata")

# Load saved variograms (vgm objects)
load("/Users/mikewhitfield/Brockett-paper/birkhowe_d1_vgm.RData")
load("/Users/mikewhitfield/Brockett-paper/birkhowe_d2_vgm.Rdata")
load("/Users/mikewhitfield/Brockett-paper/birkhowe_d3_vgm.Rdata")

load("/Users/mikewhitfield/Brockett-paper/lowsnab_d1_vgm.Rdata")
load("/Users/mikewhitfield/Brockett-paper/lowsnab_d2_vgm.Rdata")
load("/Users/mikewhitfield/Brockett-paper/lowsnab_d3_vgm.Rdata")

# Kriging with cross-validation
birkhowe_d1.cv <- krige.cv(formula = totC_mean_mass_vol ~ elevation + moisture + hls_plot,
                           birkhowe_d1_sub,
                           model = birkhowe_d1.vgm)

birkhowe_d1.stats <- extract_krige_stats(birkhowe_d1.cv)

birkhowe_d2.cv <- krige.cv(formula = totC_mean_mass_vol ~ elevation + moisture + hls_plot,
                           birkhowe_d2_sub,
                           model = birkhowe_d2.vgm)

birkhowe_d2.stats <- extract_krige_stats(birkhowe_d2.cv)

birkhowe_d3.cv <- krige.cv(formula = totC_mean_mass_vol ~ elevation + moisture + hls_plot,
                           birkhowe_d3_sub,
                           model = birkhowe_d3.vgm)

birkhowe_d3.stats <- extract_krige_stats(birkhowe_d3.cv)

lowsnab_d1.cv <- krige.cv(formula = totC_mean_mass_vol ~ elevation + moisture + hls_plot,
                          lowsnab_d1_sub,
                          model = lowsnab_d1.vgm)

lowsnab_d1.stats <- extract_krige_stats(lowsnab_d1.cv)

lowsnab_d2.cv <- krige.cv(formula = totC_mean_mass_vol ~ elevation + moisture + hls_plot,
                          lowsnab_d2_sub,
                          model = lowsnab_d2.vgm)

lowsnab_d2.stats <- extract_krige_stats(lowsnab_d2.cv)

lowsnab_d3.cv <- krige.cv(formula = totC_mean_mass_vol ~ elevation+ moisture + hls_plot,
                          lowsnab_d3_sub,
                          model = lowsnab_d3.vgm)

lowsnab_d3.stats <- extract_krige_stats(lowsnab_d3.cv)
