### Scomberomorus cavalla pcrit estimation and plotting ####
# inputs - weight, temperature, oxygen
# intermediate calculations - log(V), n , Eo, est SE
# outputs- pcrit, prob po2 > pcrit, phi

rm(list = ls())
library(dplyr)
library(MASS)
library(TMB)
library(ggplot2)
library(gridExtra)
library(readxl)
library(Matrix)
library(worrms)
library(purrr)
library(gsw)
library(cmocean)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# Setup Code ####
## load functions ####
# source("code/helper/fit_model_funs.R")
kelvin <- function(x) {x + 273.15}
find_index <- function(x, y) {y <- which(y == x)}

#' Convert Dissolved Oxygen from mg/L to umol/kg
#' @param DO_mg_L   Dissolved oxygen concentration in mg/L
#' @param temp_c    In situ water temperature in degrees Celsius
#' @param sal_psu   Practical Salinity (PSU)
#' @param depth_m   Sampling depth in meters
#' @param lat       Latitude in decimal degrees (e.g., 27.0 for Gulf of Mexico)
#' @return Dissolved oxygen in umol/kg
# convert_DO_mgL_to_umol_kg <- function(DO_mg_L, temp_c, sal_psu, depth_m, lat = 27.0) {
#   # 1. Convert depth to pressure in dbar
#   press_dbar <- gsw_p_from_z(-depth_m, latitude = lat)
  
#   # 2. Convert Practical Salinity to Absolute Salinity (g/kg)
#   SA <- gsw_SA_from_SP(sal_psu, press_dbar, longitude = -90.0, latitude = lat)
  
#   # 3. Convert in situ temperature to Conservative Temperature (°C)
#   CT <- gsw_CT_from_t(SA, temp_c, press_dbar)
  
#   # 4. Calculate in situ density (kg/m^3)
#   rho_kg_m3 <- gsw_rho(SA, CT, press_dbar)
  
#   # 5. Density in kg/L
#   rho_kg_L <- rho_kg_m3 / 1000
  
#   # 6. Molecular weight of O2 (g/mol = mg/mmol)
#   mw_O2 <- 31.9988
  
#   # Convert: (mg/L / (31.9988 mg/mmol)) * 1000 umol/mmol / rho_kg_L
#   DO_umol_kg <- (DO_mg_L / mw_O2) * (1000 / rho_kg_L)
  
#   return(DO_umol_kg)
# }
convert_DO_mgL_to_umol_kg <- function(DO_mg_L, temp_c, sal_psu, depth_m, lat = 27.0, lon = -90.0) {
  press_dbar <- gsw_p_from_z(-depth_m, latitude = lat)
  SA <- gsw_SA_from_SP(sal_psu, press_dbar, longitude = lon, latitude = lat)
  CT <- gsw_CT_from_t(SA, temp_c, press_dbar)
  
  rho_kg_L <- gsw_rho(SA, CT, press_dbar) / 1000
  mw_O2 <- 31.9988
  
  DO_umol_kg <- (DO_mg_L / mw_O2) * (1000 / rho_kg_L)
  return(DO_umol_kg)
}

# convert_umol_kg_to_kPa <- function(salinity_psu, depth_m, longitude, latitude, temp_c, O2_umol_kg) {
#   gas_const = 8.31
#   partial_molar_vol = 0.000032
#   boltz = 0.000086173324
  
#   press_dbar <- gsw_p_from_z(-depth_m, latitude = latitude)
#   SA = gsw_SA_from_SP(salinity_psu, press_dbar, longitude, latitude) #absolute salinity for pot T calc
#   pt = gsw_pt_from_t(SA, temp_c, press_dbar) #potential temp at a particular depth
#   O2_Sat0 = gsw_O2sol_SP_pt(salinity_psu, pt)

#   press = exp(depth_m * 10000 * partial_molar_vol / gas_const / kelvin(temp_c))
#   O2_satdepth = O2_Sat0 * press

#   #solubility at p=0
#   sol0 = O2_Sat0 / 0.209
#   sol_Dep = sol0 * press
#   po2 = O2_umol_kg / sol_Dep
#   po2 <- po2 * 101.325 # convert to kPa

#   return(po2)
# }
convert_umol_kg_to_kPa <- function(salinity_psu, depth_m, longitude, latitude, temp_c, O2_umol_kg) {
  R <- 8.3144626          # Universal gas constant (J / mol*K)
  Vm_O2 <- 3.2e-5         # Partial molar volume of O2 (m^3 / mol)
  x_O2 <- 0.20946         # Mole fraction of O2 in dry air
  P_atm_kPa <- 101.325    # Standard atmospheric pressure (kPa)
  temp_K <- temp_c + 273.15
  
  # 1. Surface equilibrium O2 concentration (umol/kg) at 1 atm moist air
  press_dbar <- gsw_p_from_z(-depth_m, latitude = latitude)
  SA <- gsw_SA_from_SP(salinity_psu, press_dbar, longitude = longitude, latitude = latitude)
  pt <- gsw_pt_from_t(SA, temp_c, press_dbar)
  O2_sat0 <- gsw_O2sol_SP_pt(salinity_psu, pt)
  
  # 2. Water vapor pressure (kPa) via Weiss & Price (1980) to get exact surface pO2
  # Accounts for salinity depression of vapor pressure
  pH2O_atm <- exp(24.4543 - 67.4509 * (100 / temp_K) - 4.8489 * log(temp_K / 100) - 0.000544 * salinity_psu)
  pO2_sat0_kPa <- x_O2 * (P_atm_kPa - (pH2O_atm * P_atm_kPa))
  
  # 3. Henry's law solubility at surface (umol / kg / kPa)
  sol_0 <- O2_sat0 / pO2_sat0_kPa
  
  # 4. Adjust solubility for hydrostatic pressure (negative exponent!)
  # 1 dbar = 10,000 Pa
  sol_depth <- sol_0 * exp((-press_dbar * 1e4 * Vm_O2) / (R * temp_K))
  
  # 5. Calculate in situ pO2 in kPa
  po2_kPa <- O2_umol_kg / sol_depth
  return(po2_kPa)
}

filter_data <- function(data.df) {
  data.df$EstMethod_Metric <- tolower(data.df$EstMethod_Metric)
  data.df <- dplyr::filter(data.df, !Method == "unknown", !EstMethod_Metric == "unknown")
  ## Keep only oxygen consumption methods
  data.df <- dplyr::filter(data.df, Method == "OxygenConsumption")
  return(data.df)
} ### end filter_data

load_data <- function() {
  all.dat <- readRDS(file = "data/alldata_taxonomy.RDS")
  
  # remove data with missing body size
  all.dat <- all.dat %>%
    dplyr::filter(!is.na(W))

  # change name of incertae sedis orders to one name
  tmp.index <- which(all.dat$Order == "Eupercaria incertae sedis")
  all.dat$Order[tmp.index] <- "Eupercaria"
  tmp.index <- which(all.dat$Order == "Ovalentaria incertae sedis")
  all.dat$Order[tmp.index] <- "Ovalentaria"
  # remove midwater crustacean 'Pleuroncodes planipes', as authors believed it lived anaerobically
  all.dat <- dplyr::filter(all.dat, Species !="Pleuroncodes planipes")
  
  return(all.dat)
} ### end load_data

fit_model_augmented_taxa <- function(fitnew = F) {
    if (fitnew) {
      library(dplyr)
      library(MASS)
      library(TMB)
      library(ggplot2)
      library(gridExtra)
      library(readxl)
      library(Matrix)
      
      conflicted::conflict_prefer("select", "dplyr")
      conflicted::conflict_prefer("filter", "dplyr")
      
      # Setup Code ####
      ## load functions ####
      source("code/helper/fit_model_funs.R")
      
      ## load data ####
      all.dat <- load_data()
      all.dat <- filter_data(all.dat)
      augmented_df <- augment_taxa(all.dat)
      
      # Do all of the setup stuff on data for TMB
      method_mat <- model.matrix( ~ EstMethod_Metric , all.dat)
      n_methods <- ncol(method_mat) - 1
      kb <-  8.617333262145E-5
      tref <- 15
      wref <- 5
      all.dat$W <- all.dat$W / wref
      all.dat$inv.temp <- (1 / kb) * (1 / (all.dat$Temp + 273.15) - 1 / (tref + 273.15))
      all.dat$minuslogpo2 <- -log(all.dat$Pcrit)
      
      taxa.list <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
      
      taxa.info <- make_taxa_tree_augmented(augmented_df, taxa.list)
      
      ParentChild_gz <- taxa.info$ParentChild_gz
      PC_gz <- taxa.info$PC_gz
      g_i <- taxa.info$g_i
      g_i_i <- taxa.info$g_i_i
      n_k <- taxa.info$n_k
      n_j <- taxa.info$n_j
      n_g <- taxa.info$n_g
      n_i <- taxa.info$n_i
      # which rows of taxa ParentChild_gz are attached to "real" species in data
      spc_in_PC_gz <- taxa.info$spc_in_PC_gz
      ## Setup TMB ####
      data <- list(
        PC_gz = PC_gz,
        g_i = g_i - 1,
        invtemp = all.dat$inv.temp,
        logW = log(all.dat$W),
        taxa_id = g_i_i - 1,
        minuslogpo2 = -log(all.dat$Pcrit),
        spc_in_PCgz = spc_in_PC_gz - 1,
        method_mat = method_mat[, -1]
      )
      
      parameters = list(
        alpha_j = c(0, 0, 0),
        L_z = rep(1, 6),
        log_lambda = rep(-1, length(unique(PC_gz[, 2])) - 1),
        beta_gj = matrix(0, nrow = n_g, ncol = n_j),
        beta_method = 0,
        logsigma = 0
      )
      Random <- c("beta_gj")
      model <- "hierarchical_mi"
      compile(paste0("code/TMB/", model, ".cpp"))
      dyn.load(dynlib(paste0("code/TMB/", model)))
      
      ## Run TUMB ####
      obj <-
        MakeADFun(
          data = data,
          parameters = parameters,
          DLL = model,
          random = Random,
          silent = TRUE
        )
      opt <- nlminb(obj$par, obj$fn, obj$gr)
      rep = sdreport(obj,
                     getReportCovariance = TRUE,
                     getJointPrecision = TRUE)
      # save rep and ParentChild_gz in analysis/
      saveRDS(object = list(rep = rep, ParentChild_gz = ParentChild_gz),
              file = "analysis/model_fit_augmented.RDS")
      return(list(rep = rep, ParentChild_gz = ParentChild_gz) )
    }
    if (!fitnew) {
      # check to see if file exists in analysis subfolder
      check.file <- file.exists("analysis/model_fit_augmented.RDS")
      if (!check.file)
        stop("model_fit_augmented.RDS not found. Re-run setting fitnew = TRUE")
      if (check.file) {
        model.fit <- readRDS("analysis/model_fit_augmented.RDS")
        return(model.fit)
      }
    }
  } ### end fit_model_augmented_taxa

# lookup lowest taxonomic group represntation ####
  lookup_taxonomic_group <- function(taxa.name, all.dat, ParentChild_gz) {
    # lookup taxa in sealifebase
    # look to see if species is in database
    taxa_in_data <- taxa.name %in% all.dat$Species
    # do this if taxa is within the database 'alldat'
    if (taxa_in_data) {
      # lookup row
      row_in_PC <- grep(x = ParentChild_gz$ChildName, pattern = taxa.name)
      cat("This taxa is in the dataset")
      cat("\n")
      return(row_in_PC)
    }
    
    # do this if taxa is not in data (lookup lowest group represented)
    if (!taxa_in_data) {
      # lookup full taxonomy of taxa
      worms_lookup <- try(as.data.frame(wm_records_name(taxa.name)), silent = T)
      if (class(worms_lookup) == "try-error") {
        stop("taxa name not found in WoRMS database.  Check spelling")
        
      } else {
        # reduce df to include only Phylum, Class, etc. Remove unaccepted taxonomy
        taxa_lookup <- as.data.frame(worms_lookup) %>%
          dplyr::filter(status == "accepted") 
        
        # if there are no accepted names, return an error
        if (nrow(taxa_lookup)==0) stop(paste0("Error: Supplied taxa name ", taxa.name, " is not an accepted name in the WoRMS database"))
        
        # Keep only the columns for taxonomy
        
        taxa_loopup <-  as.data.frame(worms_lookup)  %>%
          dplyr::select(phylum, class, order, family, genus, scientificname)
        # Identify the lowest taxonomic order that is represented in the df all.dat
        taxa_levels <- c("phylum", "class", "order", "family", "genus")
        
        # make taxa_df from all.dat that only has taxonomic information
        taxa_df <- all.dat %>%
          dplyr::select(Phylum, Class, Order, Family, Genus)
        colnames(taxa_df) = tolower(colnames(taxa_df))
        
        
        # Start from most specific and work upward to find lowest shared rank
        lowest_shared_rank <- NULL
        lowest_shared_value <- NULL
        
        for (i in seq_along(rev(taxa_levels))) {
          level <- rev(taxa_levels)[i]
          value <- taxa_lookup[[level]]
          
          if (!is.null(value) && level %in% colnames(taxa_df)) {
            if (value %in% taxa_df[[level]]) {
              lowest_shared_rank <- level
              lowest_shared_value <- value
              break
            }
          }
        }
        
        # If match was found, construct the placeholder name accordingly
        if (!is.null(lowest_shared_rank)) {
          # Create the placeholder row string based on the level
          ranks_to_use <- taxa_levels[1:which(taxa_levels == lowest_shared_rank)]
          base_taxa <- sapply(ranks_to_use, function(r)
            taxa_lookup[[r]])
          
          # Add placeholder at appropriate level
          placeholder_suffix <- list(
            genus = "sp",
            family = "blankgenusinfamily_blankgenusinfamily sp",
            order = "blankfamilyinorder_blankgenusinorder_blankgenusinorder sp",
            class = "blankorderinclass_blankfamilyinclass_blankgenusinclass_blankgenusinclass sp",
            phylum = "blankclassinphylum_blankorderinphylum_blankfamilyinphylum_blankgenusinphylum_blankgenusinphylum sp"
          )
          
          placeholder <- placeholder_suffix[[lowest_shared_rank]]
          
          # Combine taxon string
          collapsed_taxon <- paste(c(base_taxa, placeholder), collapse = "_")
          
          # Lookup row in ParentChild_gz
          row_in_PC <- which(ParentChild_gz$ChildName == collapsed_taxon)
          
          if (length(row_in_PC) > 0) {
            cat( "using estimate from the following taxa: \n" )
            cat( collapsed_taxon )
            cat("\n")
            
            return(row_in_PC)
          } else {
            message("no matching taxonomic rank found in data")
            return(NULL)
          }
        }
      }
    }
  } ### end lookup_taxonomic_group

 # Estimate pcrit for taxa using full covariance ####
  estimate_taxa <- function(taxa.name, w , temperature, method = "smr", psigma = 0, fitnew = F) {
    
    # input checking
    method <- tolower(method)
    if (!method %in% c("routine", "smr")) stop("input argument method must be either 'routine' or 'smr'")

    # first load taxa
    all.dat <- load_data()
    all.dat <- filter_data(all.dat)
    
    # load model or fit new one
    model.fit <- fit_model_augmented_taxa(fitnew = fitnew)
    ParentChild_gz <- model.fit$ParentChild_gz
    rep = model.fit$rep
    # lookup rownumber for listed taxa within ParentChild matrix
    lookup_row_number <- lookup_taxonomic_group(taxa.name = taxa.name,
                                                all.dat = all.dat,
                                                ParentChild_gz = ParentChild_gz)
    
    # Extract fitted parameters from fitted object
    re <- summary(rep, "random")
    fixef <- summary(rep, "fixed")
    names <- rownames(summary(rep, "all"))
    
    # Lookup estimated logsigma
    logsigma <- fixef["logsigma", "Estimate"]
    sigma <- exp(logsigma)
    
    n_j <- 3 # number of traits
    n_g <- nrow(ParentChild_gz) # number of taxagroups
    
    beta_mle <- matrix(re[grep(rownames(re), pattern = "beta_gj"), 1],
                       nrow = n_g,
                       ncol = 3,
                       byrow = F)
    
    beta_se <- matrix(re[grep(rownames(re), pattern = "beta_gj"), 2],
                      nrow = n_g,
                      ncol = 3,
                      byrow = F)
    
    # get all rownumbers for beta_gj indices
    beta_gj_positions <- which(names == "beta_gj")
    all_beta_gj_matrix_indices <- matrix(beta_gj_positions,
                                         nrow = n_g,
                                         ncol = n_j,
                                         byrow = FALSE)
    
    # Extract taxonomic-specific values
    beta_mle_taxa <- beta_mle[lookup_row_number, ]
    beta_se_taxa <- beta_se[lookup_row_number, ]
    
    
    # Lookup beta_method
    beta_method_index <- which(names == "beta_method")
    beta_method <- fixef["beta_method", 1]
    beta_method_se <- fixef["beta_method", 2]
    
    # Do calculations to approximate the variance- covariance matrix
    # I assume no covariance between traits and beta_method
    
    # get sparse precision matrix
    J <- rep$jointPrecision
    # convert to var-covar
    
    V_joint <- solve(J)
    
    beta_gj_indices <- all_beta_gj_matrix_indices[lookup_row_number,]
    V <- V_joint[c(beta_gj_indices, beta_method_index), c(beta_gj_indices, beta_method_index)]
    
    # some placeholder settings
    wref <- 5
    tref <- 15
    kb <-  8.617333262145E-5
    
    betas <- as.vector(c(beta_mle_taxa, beta_method))
    logw <- log(w / wref)
    inv.temp <- 1 / kb * (1 / kelvin(temperature) - 1 / kelvin(tref))
    
    # Set up x_predict
    if (method == "smr") x_predict <- as.vector(c(1, -logw, -inv.temp, 1))
    if (!method == "smr") x_predict <- as.vector(c(1, -logw, -inv.temp, 0))
    
    # calculate prediction and SE of prediction
    log_pcrit_predict <-  x_predict %*% betas
    log_pcrit_se <- as.numeric(sqrt(t(x_predict) %*% V %*% x_predict + psigma * sigma^2))
    
    # format results nicely
    parameter_estimates <- cbind(c(beta_mle_taxa, beta_method),
                                 c(beta_se_taxa, beta_method_se))
    rownames(parameter_estimates) <- c("log(V)", "n", "Eo", "beta_method")
    colnames(parameter_estimates) <- c("Estimate", "SE")
    
    log_pcrit_estimate <- c(logpcrit = log_pcrit_predict, se = log_pcrit_se)
    return(list(parameters = parameter_estimates, log_pcrit = log_pcrit_estimate, var_covar = V))
  } ### end estimate_taxa

 # Calculate the probability that an input value po2 is greater than pcrit ####
  calc_p_po2 <- function(po2, temp, w, pcrit_type = "smr", betas, var_covar) {
    # Setup fixed values
    wref <- 5
    tref <- 15
    kb <-  8.617333262145E-5
    # Transform w and temperature
    inv_temperature <- (1/ kb) * (1 / kelvin(temp) - 1 / kelvin(tref) )
    logw <- log(w / wref )
    # calculate prediction x's
    if(pcrit_type == "smr") x_predict <- as.vector(c(1, -logw, -inv_temperature, 1))
    if(!pcrit_type == "smr") x_predict <- c(1, -logw, -inv_temperature, 0)
    # generate prediction and standard error
    log_pcrit_predict <- x_predict %*% betas
    log_pcrit_se <-as.numeric(sqrt( t(x_predict) %*% var_covar %*% x_predict ))
    # calculate probability that po2 exceeds pcrit
    p_po2_exceeds_pcrit <- pnorm(q = log(po2), mean = log_pcrit_predict, sd = log_pcrit_se)
    
    return(data.frame(p_po2_exceeds_pcrit = p_po2_exceeds_pcrit,
      log_pcrit_predict = log_pcrit_predict,
      log_pcrit_se = log_pcrit_se))
  } ### end calc_p_po2

update_df <- function(species_df, taxa.name, w, pcrit_type = "smr"){ #, rep, ParentChild_gz) {
  # get parameter estimates
  taxa_estimates <- estimate_taxa(taxa.name, w = w , temperature = 10, method = pcrit_type)
  # use pmap to run calc_p_po2
  betas <- taxa_estimates$parameters[,1]
  var_covar <- taxa_estimates$var_covar
  p_po2_out <- pmap(list(po2 = species_df$po2, temp = species_df$temperature_C),
    ~calc_p_po2(po2 = ..1, temp = ..2, 
      w = w,  pcrit_type = pcrit_type, betas = betas, var_covar = var_covar)) |>
    list_rbind()
  species_df$p_pcrit <- p_po2_out$p_po2_exceeds_pcrit
  species_df$phi <- species_df$po2 / exp(p_po2_out$log_pcrit_predict)
  species_df$pcrit <- exp(p_po2_out$log_pcrit_predict)
return(species_df)
} ### end update_df

# get fitted model
# model.fit <- fit_model_augmented_taxa(fitnew = F)

setwd('c:/Users/brendan.turley/Documents/R_projects/metabolic_index_cmp')
kgm_df <- data.frame(po2 = rep(5,21), temperature_C = seq(20,40,1))
update_df(kgm_df, taxa.name = 'Scomberomorus cavalla', w = 3000, pcrit_type = "smr")

estimate_taxa(taxa.name = 'Scomberomorus cavalla', w = 3000, temperature = 10, method = "smr")

load('c:/Users/brendan.turley/Documents/R_projects/KGM-environ-sensitivity/data/kmk_smdat.RData')

### surface
do_umol_kg <- convert_DO_mgL_to_umol_kg(DO_mg_L = kmk_pos$OXYSURF, temp_c = kmk_pos$TEMPSURF, 
  sal_psu = kmk_pos$SALSURF, depth_m = 2, lat = kmk_pos$lat, lon = kmk_pos$lon)
po2_kpa <- convert_umol_kg_to_kPa(salinity_psu = kmk_pos$SALSURF, depth_m = 2, longitude = kmk_pos$lon,
  latitude = kmk_pos$lat, temp_c = kmk_pos$TEMPSURF, O2_umol_kg = do_umol_kg)

setwd('c:/Users/brendan.turley/Documents/R_projects/metabolic_index_cmp')
kgm_df_surf <- data.frame(po2 = po2_kpa, temperature_C = kmk_pos$TEMPSURF,
  lat = kmk_pos$lat, lon = kmk_pos$lon, DO_mg_L = kmk_pos$OXYSURF, month = kmk_pos$month)
kgm_df_surf <- update_df(kgm_df_surf, taxa.name = 'Scomberomorus cavalla', w = 300, pcrit_type = "smr")

hist(kgm_df_surf$p_pcrit)
hist(kgm_df_surf$pcrit)
hist(kgm_df_surf$phi)
hist(po2_kpa)
plot(kgm_df_surf$p_pcrit,kgm_df_surf$phi)
plot(kgm_df_surf$temperature_C,kgm_df_surf$phi)
abline(h = 1, col = "red", lty = 2)
plot(kgm_df_surf$po2,kgm_df_surf$phi)
abline(h = 1, col = "red", lty = 2)
range(kgm_df_surf$phi, na.rm = TRUE)
range(kgm_df_surf$pcrit, na.rm = TRUE)

### bottom
do_umol_kg <- convert_DO_mgL_to_umol_kg(DO_mg_L = kmk_pos$OXYMAX, temp_c = kmk_pos$TEMP_BOT, 
  sal_psu = kmk_pos$SALMAX, depth_m = kmk_pos$depth, lat = kmk_pos$lat, lon = kmk_pos$lon)
po2_kpa <- convert_umol_kg_to_kPa(salinity_psu = kmk_pos$SALMAX, depth_m = kmk_pos$depth, longitude = kmk_pos$lon,
  latitude = kmk_pos$lat, temp_c = kmk_pos$TEMP_BOT, O2_umol_kg = do_umol_kg)

setwd('c:/Users/brendan.turley/Documents/R_projects/metabolic_index_cmp')
kgm_df_bot <- data.frame(po2 = po2_kpa, temperature_C = kmk_pos$TEMP_BOT, 
  lat = kmk_pos$lat, lon = kmk_pos$lon, DO_mg_L = kmk_pos$OXYMAX, month = kmk_pos$month)
kgm_df_bot <- update_df(kgm_df_bot, taxa.name = 'Scomberomorus cavalla', w = 300, pcrit_type = "smr")

hist(kgm_df_bot$p_pcrit)
hist(kgm_df_bot$phi)
hist(po2_kpa)
plot(kgm_df_bot$p_pcrit,kgm_df_bot$phi)
plot(kgm_df_bot$temperature_C,kgm_df_bot$phi)
abline(h = 1, col = "red", lty = 2)
plot(kgm_df_bot$po2,kgm_df_bot$phi)
abline(h = 1, col = "red", lty = 2)
range(kgm_df_bot$phi, na.rm = TRUE)
range(kgm_df_bot$pcrit, na.rm = TRUE)


### what combination of temperature and po2 leads to phi < 1?
range(kmk_pos$OXYMAX, na.rm = T)
range(kmk_pos$OXYSURF, na.rm = T)
range(kmk_pos$TEMP_BOT, na.rm = T)
range(kmk_pos$TEMPSURF, na.rm = T)
range(kmk_pos$SALMAX, na.rm = T)
range(kmk_pos$SALSURF, na.rm = T)

do_t <- expand.grid(DO_mg_L = seq(0,12,.1), temperature_C = seq(15,40,1))
do_t |>
  pmap(function(DO_mg_L, temperature_C) {
    do_umol_kg <- convert_DO_mgL_to_umol_kg(DO_mg_L = DO_mg_L, temp_c = temperature_C, 
      sal_psu = 35, depth_m = 2)
    po2_kpa <- convert_umol_kg_to_kPa(salinity_psu = 35, depth_m = 2, longitude = -89,
      latitude = 28, temp_c = temperature_C, O2_umol_kg = do_umol_kg)
    data.frame(po2 = po2_kpa, temperature_C = temperature_C)
  }) |>
  list_rbind() -> kgm_df_grid_do_t
setwd('c:/Users/brendan.turley/Documents/R_projects/metabolic_index_cmp')
kgm_df_grid_do_t <- update_df(kgm_df_grid_do_t, taxa.name = 'Scomberomorus cavalla', w = 1000, pcrit_type = "smr")

# hist(kgm_df_grid_do_t$p_pcrit)
# hist(kgm_df_grid_do_t$phi)
# plot(kgm_df_grid_do_t$p_pcrit,kgm_df_grid_do_t$phi)

phi_m_t <- matrix(kgm_df_grid_do_t$phi, nrow = length(unique(do_t$DO_mg_L)), 
  ncol = length(unique(do_t$temperature_C)))
pcrit_m_t <- matrix(kgm_df_grid_do_t$p_pcrit, nrow = length(unique(do_t$DO_mg_L)), 
  ncol = length(unique(do_t$temperature_C)))



do_s <- expand.grid(DO_mg_L = seq(0,12,.1), temperature_C = 28, salinity_psu = seq(5,38,.5))
do_s |>
  pmap(function(DO_mg_L, temperature_C, salinity_psu) {
    do_umol_kg <- convert_DO_mgL_to_umol_kg(DO_mg_L = DO_mg_L, temp_c = temperature_C, 
      sal_psu = salinity_psu, depth_m = 2)
    po2_kpa <- convert_umol_kg_to_kPa(salinity_psu = salinity_psu, depth_m = 2, longitude = -89,
      latitude = 28, temp_c = temperature_C, O2_umol_kg = do_umol_kg)
    data.frame(po2 = po2_kpa, temperature_C = temperature_C, salinity_psu = salinity_psu)
  }) |>
  list_rbind() -> kgm_df_grid_do_s
setwd('c:/Users/brendan.turley/Documents/R_projects/metabolic_index_cmp')
kgm_df_grid_do_s <- update_df(kgm_df_grid_do_s, taxa.name = 'Scomberomorus cavalla', w = 1000, pcrit_type = "smr")

# hist(kgm_df_grid_do_s$p_pcrit)
# hist(kgm_df_grid_do_s$phi)
# plot(kgm_df_grid_do_s$p_pcrit,kgm_df_grid_do_s$phi)

phi_m_s <- matrix(kgm_df_grid_do_s$phi, nrow = length(unique(do_s$DO_mg_L)), 
  ncol = length(unique(do_s$salinity_psu)))
pcrit_m_s <- matrix(kgm_df_grid_do_s$p_pcrit, nrow = length(unique(do_s$DO_mg_L)), 
  ncol = length(unique(do_s$salinity_psu)))



setwd("c:/Users/brendan.turley/Documents/R_projects/kgm-environ-sensitivity/figures")
png('kgm_mi.png', width = 8, height = 4, units = 'in', res = 300)
par(mfrow = c(1,2), mar = c(4, 4, 1, 1), oma = c(0, 0, 3, 0))

image(unique(do_t$DO_mg_L), unique(do_t$temperature_C),  z = phi_m_t,
  xlab = "DO (mg/L)", ylab = "Temperature (C)", 
  # main = expression(paste('King Mackerel Metabolic Index (',phi, " = 1)")),
  col = cmocean('solar')(100))
# contour(unique(do_t$DO_mg_L), unique(do_t$temperature_C), z = phi_m_t, add = TRUE, col = "white",
#   levels = 1)
contour(unique(do_t$DO_mg_L), unique(do_t$temperature_C), z = 1 - pcrit_m_t, add = TRUE, col = "white",
  levels = seq(0,.9,.1))
points(kmk_pos$OXYMAX, kmk_pos$TEMP_BOT, bg = rgb(t(col2rgb("orange"))/255, alpha = 0.7), pch = 21, cex = .7)
points(kmk_pos$OXYSURF, kmk_pos$TEMPSURF, bg = rgb(t(col2rgb("dodgerblue"))/255, alpha = 0.7), pch = 24, cex = .7)
legend("topright", legend = c("Bottom", "Surface"), 
  pch = c(21, 24), pt.bg = c(rgb(t(col2rgb("orange"))/255, alpha = 0.7), 
  rgb(t(col2rgb("dodgerblue"))/255, alpha = 0.7)), cex = .7)

image(unique(do_s$DO_mg_L), unique(do_s$salinity_psu),  z = phi_m_s,
  xlab = "DO (mg/L)", ylab = "Salinity (PSU)",
  col = cmocean('solar')(100))
# contour(unique(do_s$DO_mg_L), unique(do_s$salinity_psu), z = phi_m_s, add = TRUE, col = "white",
#   levels = 1)
contour(unique(do_s$DO_mg_L), unique(do_s$salinity_psu), z = 1 - pcrit_m_s, add = TRUE, col = "white",
  levels = seq(0,.9,.1))
points(kmk_pos$OXYMAX, kmk_pos$SALMAX, bg = rgb(t(col2rgb("orange"))/255, alpha = 0.7), pch = 21, cex = .7)
points(kmk_pos$OXYSURF, kmk_pos$SALSURF, bg = rgb(t(col2rgb("dodgerblue"))/255, alpha = 0.7), pch = 24, cex = .7)

title(main = expression(paste('King Mackerel Metabolic Index (',phi, ")")), outer = T)

dev.off()


plot(kgm_df_surf$lon, kgm_df_surf$lat, pch = 21, asp = 1)#, bg = kgm_df_surf$phi)
with(subset(kgm_df_bot, phi < 1), 
  points(lon, lat, cex = 3,
    pch = 21, bg =  rgb(t(col2rgb("orange"))/255, alpha = 0.7)))
with(subset(kgm_df_surf, phi < 1), 
  points(lon, lat, cex = 1.5,
    pch = 21, bg =  rgb(t(col2rgb("dodgerblue"))/255, alpha = 0.7)))

boxplot(phi ~ month, data = kgm_df_surf, main = "Surface Phi by Month", xlab = "Month", ylab = "Phi")
boxplot(phi ~ month, data = kgm_df_bot, main = "Bottom Phi by Month", xlab = "Month", ylab = "Phi")
