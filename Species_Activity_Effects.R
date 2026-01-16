################################################################################
### Project: Thesis Ch 1. Evaluating which eye metrics best explain variation among species in their responses to light pollution
### Updated: 1/13/26
### Author(s): Clinton D. Francis (main) & Eva M. Moylan
################################################################################

### OBJECTIVES ####
# Calculate species-specific effect sizes for all species considered in Pease and Gilbert (2025) Science.
# Perform analysis of species-specific chances in activity in the morning and evenining in response to light pollution. 
# The model uses the same structure as originally used in Pease and Gilbert, but reduces the nested random effects to just the "sp_cell5_week" variable. 

# Load packages
library(tidyverse)
library(here)
library(janitor)
library(glmmTMB)
library(gt)
library(performance)
library(broom.mixed)
library(ggeffects)
library(patchwork)
library(clootl)
library(phylolm)
library(ape)
library(geiger)
library(picante)
library(phytools)
library(phyr) 
library(ggtree)
#devtools::install_github("xiangpin/ggtreeExtra")
library(ggtreeExtra)
library(ggnewscale)
library(colorspace)
library(treeio)

# key for joining up birdweather, EltonTraits, etc.
key <- readr::read_csv(here("Data/ActivityData/species_keys","birdweather_elton_botw_name_key.csv"))

# AvoNET trait database
avo <- readr::read_csv(here("Data/ActivityData/traits", "avonet.csv")) |> 
  dplyr::mutate(lat_range = Max.Latitude - Min.Latitude)

# table indicating whether each species nests in cavities/burrows (1) or not (0)
cavity <- readr::read_csv(here("Data/ActivityData/traits","cavity.csv")) |> 
  dplyr::rename(sci_name = sci_name_bw)

# EltonTraits database
elton <- read.delim(here::here("Data/ActivityData/traits/elton.txt")) |> 
  janitor::clean_names() |> 
  dplyr::select( family = bl_family_latin, 
                 sci_name_elton = scientific,
                 mass = body_mass_value,
                 starts_with("diet"),
                 starts_with("for_strat")) |> 
  dplyr::select(-diet_source, -diet_certainty, -diet_entered_by, 
                -for_strat_source, -for_strat_spec_level, -for_strat_entered_by) |> 
  dplyr::right_join( key ) |> 
  tibble::as_tibble() |> 
  dplyr::select(sci_name = sci_name_bw, 
                family, 
                mass,
                starts_with("diet"), 
                starts_with("for_strat")) |> 
  dplyr::distinct()


# import filtered and processed birdweather data
setwd(here::here("Data/ActivityData/vocalization_activity"))
load(here::here("Data/ActivityData/vocalization_activity/cessation_data_conf_0.75_det_100_grid_10.RData")) # cessation of evening vocalization
load(here::here("Data/ActivityData/vocalization_activity/onset_data_conf_0.75_det_100_grid_10.RData")) # onset of vocalization in the morning

# updated name key to correct for some idiosyncracies between AvoNet and BirdLife names
key2 <- key |> 
  dplyr::mutate(sci_name_avo = sci_name_botw) |> 
  dplyr::mutate(sci_name_avo = ifelse(sci_name_avo == "Curruca communis", 
                                      "Sylvia communis", sci_name_avo)) |> 
  dplyr::mutate(sci_name_avo = ifelse(sci_name_avo == "Curruca curruca", "Sylvia curruca",
                                      sci_name_avo)) |> 
  dplyr::mutate(sci_name_avo = ifelse(sci_name_avo == "Ortygornis sephaena",
                                      "Dendroperdix sephaena", sci_name_avo)) |> 
  dplyr::mutate(sci_name_avo = ifelse(sci_name_avo == "Grus canadensis", "Antigone canadensis",
                                      sci_name_avo)) |> 
  dplyr::mutate(sci_name_avo = ifelse(sci_name_avo  == "Tetrastes bonasia",
                                      "Bonasa bonasia", sci_name_avo))

load(here::here("Data/ActivityData/vocalization_activity/tot_voc.RData"))

eye <- readr::read_csv(here::here("Data/ActivityData/traits", "ritland_clean.csv"))

eye.sp <- key2 |> 
  dplyr::select(sci_name = sci_name_bw, sci_name_elton) |> 
  dplyr::distinct() |> 
  dplyr::left_join(
    eye |> 
      dplyr::rename(sci_name_elton = species_jetz)) |> 
  dplyr::select(sci_name, cd1) |> 
  dplyr::distinct()
# join warning message occurs here. It can be ignored

eye.fam <- eye |> 
  group_by(family_jetz) |> 
  dplyr::summarise(cd1.fam = mean(cd1),
                   count = n()) |> 
  dplyr::rename(family = family_jetz)


##################################################################################
########################### EVENING DATA PREP ####################################
##################################################################################
# final formatted dataset: evening
d_e <- final_cess |> 
  dplyr::ungroup() |>
  dplyr::rename(sci_name_bw = sci_name) |> 
  # join with avonet trait database
  dplyr::left_join(
    avo |> 
      janitor::clean_names() |> 
      dplyr::rename(sci_name_avo = species1) |> 
      dplyr::right_join(key2) |> 
      dplyr::group_by(sci_name_bw) |> 
      # have to account for some species that have multiple entries (splits/lumps)
      dplyr::summarise(hand_wing_index = mean(hand_wing_index), 
                       mass_avo = mean(mass), # body mass from avonet
                       habitat = unique(habitat), 
                       hd = unique(habitat_density),
                       migration = unique(migration), 
                       lat_range = mean(lat_range),
                       range_size = mean(range_size)) |> 
      dplyr::distinct() |> 
      dplyr::slice(1) |> 
      dplyr::ungroup()) |> 
  # covariate transformations
  dplyr::mutate(mass_avo_sc = as.numeric(scale(log(mass_avo))), # scaled and centered ln of body mass from avonet
                alan_sc = as.numeric(scale(log1p(avg_rad))), # scale ln of radiance + 1
                lat_sc = as.numeric(scale(abs(lat))), # scale absolute latitude
                value_hr = value / 60, # convert response variable from units of minutes to hours
                range_size_sc = as.numeric(scale(log(range_size))),
                lat_range_sc = as.numeric(scale(log(lat_range)))) |>  # scale ln of range size
  dplyr::rename(sci_name = sci_name_bw) |> 
  dplyr::group_by(sci_name, week) |> 
  # create group IDs for random effect groupings
  dplyr::mutate(sp_week = dplyr::cur_group_id()) |>  # species x week
  dplyr::group_by(sci_name, grid_ID_cell_5, week) |> 
  dplyr::mutate(sp_cell5_week = dplyr::cur_group_id()) |>  # species x 5 deg grid cell x week
  dplyr::group_by(sci_name, grid_ID_cell_10, week) |>  
  dplyr::mutate(sp_cell10_week = dplyr::cur_group_id()) |> # species x 10 deg grid cell x week 
  dplyr::group_by(sci_name, grid_ID_cell_15, week) |> 
  dplyr::mutate(sp_cell15_week = dplyr::cur_group_id()) |>  # species x 15 deg grid cell x week
  dplyr::ungroup() |> 
  dplyr::left_join(cavity) |> 
  dplyr::mutate(cavity = factor(cavity)) |> 
  dplyr::left_join( elton ) |> 
  # transformations of various trait columns
  dplyr::mutate(across(starts_with("for_strat"), function(x) x / 100)) |> 
  dplyr::mutate( inv_sc = as.numeric(scale(diet_inv)), 
                 ground_sc = as.numeric(scale(for_strat_ground)), 
                 under_sc = as.numeric(scale(for_strat_understory)), 
                 low_sc = as.numeric(scale (for_strat_ground + for_strat_understory)),
                 mass_sc = as.numeric(scale(log(mass))),
                 migration = factor(migration),
                 hd = factor(hd)) |> 
  dplyr::left_join( tot.voc ) |> 
  dplyr::mutate( 
    shan_raw = shan,
    shan = as.numeric(scale(shan)),
    sr_raw = sr,
    sr = as.numeric(scale(log(sr))),
    tot_raw = tot,
    tot = as.numeric(scale(log(tot)))) |> 
  dplyr::select(lat, lon, grid_ID_cell_5, date, week, sci_name, family, value_hr, avg_rad, alan_sc, lat_sc, mass_avo, mass_avo_sc, mass_sc,
                range_size, range_size_sc, lat_range, lat_range_sc, migration, hd, for_strat_ground, ground_sc, cavity, tot, tot_raw, sr, sr_raw, shan, shan_raw, sp_cell5_week, sp_week) |> 
  left_join(eye.sp) |> 
  left_join(eye.fam) |> 
  dplyr::mutate(
    cd = ifelse(!is.na(cd1), cd1, cd1.fam)) |> 
  dplyr::mutate(
    cd_raw = cd, 
    cd = as.numeric(scale(log(cd)))) # scaled and centered ln of CD (eye size)

##################################################################################
########################### MORNING DATA PREP ####################################
##################################################################################
# same data processing for morning onset
d_onset <- final |> 
  dplyr::ungroup() |>
  dplyr::rename(sci_name_bw = sci_name) |> 
  dplyr::left_join(
    avo |> 
      janitor::clean_names() |> 
      dplyr::rename(sci_name_avo = species1) |> 
      dplyr::right_join(key2) |> 
      dplyr::group_by(sci_name_bw) |> 
      dplyr::summarise(hand_wing_index = mean(hand_wing_index), 
                       mass_avo = mean(mass), # body mass from avonet
                       habitat = unique(habitat), 
                       migration = unique(migration),
                       hd = unique(habitat_density),
                       range_size = mean(range_size),
                       lat_range = mean(lat_range, na.rm = TRUE)) |> 
      dplyr::distinct() |> 
      dplyr::slice(1) |> 
      dplyr::ungroup()) |> 
  dplyr::mutate(mass_avo_sc = as.numeric(scale(mass_avo)), # scale and center ln of body mass
                alan_sc = as.numeric(scale(log1p(avg_rad))),
                lat_sc = as.numeric(scale(abs(lat))),
                value_hr = value / 60,
                range_size_sc = as.numeric(scale(log(range_size))),
                lat_range_sc = as.numeric(scale(log(lat_range))),
                migration = factor(migration),
                hd = factor(hd)) |> 
  dplyr::rename(sci_name = sci_name_bw) |> 
  dplyr::group_by(sci_name, week) |> 
  dplyr::mutate(sp_week = dplyr::cur_group_id()) |> 
  dplyr::group_by(sci_name, grid_ID_cell_5, week) |> 
  dplyr::mutate(sp_cell5_week = dplyr::cur_group_id()) |> 
  dplyr::group_by(sci_name, grid_ID_cell_10, week) |> 
  dplyr::mutate(sp_cell10_week = dplyr::cur_group_id()) |> 
  dplyr::group_by(sci_name, grid_ID_cell_15, week) |> 
  dplyr::mutate(sp_cell15_week = dplyr::cur_group_id()) |> 
  dplyr::ungroup() |> 
  dplyr::left_join(cavity) |> 
  dplyr::mutate(cavity = factor(cavity)) |> 
  dplyr::left_join( elton ) |> 
  dplyr::mutate(across(starts_with("for_strat"), function(x) x / 100)) |> 
  dplyr::mutate( inv_sc = as.numeric(scale(diet_inv)), 
                 ground_sc = as.numeric(scale(for_strat_ground)), 
                 under_sc = as.numeric(scale(for_strat_understory)), 
                 low_sc = as.numeric(scale (for_strat_ground + for_strat_understory)),
                 mass_sc = as.numeric(scale(log(mass)))) |> 
  dplyr::left_join( tot.voc ) |> 
  dplyr::mutate( 
    shan_raw = shan,
    shan = as.numeric(scale(shan)),
    sr_raw = sr,
    sr = as.numeric(scale(log(sr))),
    tot_raw = tot,
    tot = as.numeric(scale(log(tot)))) |> 
  dplyr::select(lat, lon, grid_ID_cell_5, date, week, sci_name, family, value_hr, avg_rad, alan_sc, lat_sc, mass_avo, mass_avo_sc, mass_sc,
                range_size, range_size_sc, lat_range, lat_range_sc, migration, hd, for_strat_ground, ground_sc, cavity, tot, tot_raw, sr, sr_raw, shan, shan_raw, sp_cell5_week, sp_week) |> 
  left_join(eye.sp) |> 
  left_join(eye.fam) |> 
  dplyr::mutate(
    cd = ifelse(!is.na(cd1), cd1, cd1.fam)) |> 
  dplyr::mutate(
    cd_raw = cd,
    cd = as.numeric(scale(log(cd)))) # scale and center ln of CD (eye size)


##################################################################################
########################### ANALYSES #############################################
### Function will run model for each species for both morning and evening activity
### Compiles responses into two data.frames
### Prints model diagnostics for each model to file (this might need to be able to identify which pdf is for evening and morning per species)
##################################################################################
##################################################################################

run_model_by_species_perf <- function(dat,
                                 db_path = NULL,
                                 table_name = "species_model_results",
                                 append_db = !is.null(db_path),
                                 save_plots = TRUE,
                                 plot_dir = "check_model_plots") {
  # Make sure glmmTMB is available
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    stop("Package 'glmmTMB' is required but not installed.")
  }
  
  # Optional model-check plots
  if (save_plots) {
    if (!requireNamespace("performance", quietly = TRUE)) {
      stop("Package 'performance' is required to save check_model plots.")
    }
    if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  }
  
  
  # Split data by species
  dat_split <- split(dat, dat$sci_name)
  
  # Function to fit model and extract stats for one species
  fit_one_species <- function(df_sp, sp_name) {
    n_sp <- nrow(df_sp)
    
    # Helper for NA row
    na_row <- function() {
      data.frame(
        sci_name        = sp_name,
        intercept_est   = NA_real_,
        intercept_se    = NA_real_,
        intercept_z     = NA_real_,
        alan_sc_est     = NA_real_,
        alan_sc_se      = NA_real_,
        alan_sc_z       = NA_real_,
        n               = n_sp,
        stringsAsFactors = FALSE
      )
    }
    
    # If too few rows to fit random effect, return NA
    if (n_sp < 2) return(na_row())
    
    # Fit model; catch errors gracefully
    fit <- try(
      glmmTMB::glmmTMB(
        value_hr ~ 1 + alan_sc + (1 | sp_cell5_week),
        data = df_sp
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) return(na_row())
    
    # Save check_model() plot(s) to a PDF named by sci_name (optional)
    if (save_plots) {
      # Make a filesystem-safe filename
      safe_name <- gsub("[^A-Za-z0-9_-]+", "_", sp_name)
      pdf_file  <- file.path(plot_dir, paste0(safe_name, ".pdf"))
      
      # Wrap plotting in try so model results still return even if plotting fails
      try({
        grDevices::pdf(pdf_file, width = 10, height = 8, onefile = TRUE)
        on.exit(grDevices::dev.off(), add = TRUE)
        
        p <- performance::check_model(fit)
        
        # check_model may return a patchwork/ggplot object, list, or print method
        # Just printing is the most robust way to ensure it draws to the device.
        print(p)
      }, silent = TRUE)
    }
    
    # Extract coefficient table (conditional model)
    s <- summary(fit)
    coef_tab <- s$coefficients$cond
    
    # Pull intercept
    intercept_est <- if ("(Intercept)" %in% rownames(coef_tab)) coef_tab["(Intercept)", "Estimate"] else NA_real_
    intercept_se  <- if ("(Intercept)" %in% rownames(coef_tab)) coef_tab["(Intercept)", "Std. Error"] else NA_real_
    intercept_z   <- if ("(Intercept)" %in% rownames(coef_tab)) coef_tab["(Intercept)", "z value"] else NA_real_
    
    # Pull alan_sc
    alan_sc_est <- if ("alan_sc" %in% rownames(coef_tab)) coef_tab["alan_sc", "Estimate"] else NA_real_
    alan_sc_se  <- if ("alan_sc" %in% rownames(coef_tab)) coef_tab["alan_sc", "Std. Error"] else NA_real_
    alan_sc_z   <- if ("alan_sc" %in% rownames(coef_tab)) coef_tab["alan_sc", "z value"] else NA_real_
    
    data.frame(
      sci_name        = sp_name,
      intercept_est   = unname(intercept_est),
      intercept_se    = unname(intercept_se),
      intercept_z     = unname(intercept_z),
      alan_sc_est     = unname(alan_sc_est),
      alan_sc_se      = unname(alan_sc_se),
      alan_sc_z       = unname(alan_sc_z),
      n               = n_sp,
      stringsAsFactors = FALSE
    )
  }
  
  # Apply to each species and row-bind
  res_list <- mapply(
    FUN     = fit_one_species,
    df_sp   = dat_split,
    sp_name = names(dat_split),
    SIMPLIFY = FALSE
  )
  
  results <- do.call(rbind, res_list)
  
  
  results
}

### test run with a subset of species
# Just the first three
d_e_sub <- d_e |>
  dplyr::filter(sci_name %in% c("Cossypha caffra", 
                         "Passer diffusus", 
                         "Vanellus coronatus",
                         "Streptopelia semitorquata",
                         "Bostrychia hagedash",
                         "Psittacula krameri"))
sub_test <- run_model_by_species_perf(d_e_sub) # Okay, test worked!

### run model for each response activity
setwd(here::here("Data/ActivityData/EveningSppModels"))
eveningData <- run_model_by_species_perf(d_e) # Takes ~ 2hrs
# Multiple 'rank deficient conditional' models 

setwd(here::here("Data/ActivityData/MorningSppModels"))
morningData <- run_model_by_species_perf(d_onset)
# Multiple 'rank deficient conditional' models 

# Save model files
saveRDS(morningData, here("Data/ActivityData/MorningSppModels/MorningActivModels.rds"))
saveRDS(eveningData, here("Data/ActivityData/EveningSppModels/EveningActivModels.rds"))
