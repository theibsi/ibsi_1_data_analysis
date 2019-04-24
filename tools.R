parse_significant <- function(x, tag, n_digits=3){
  #Function to parse values to a significant numbers
  
  # Exceptions
  if(stringi::stri_detect_fixed(str=tag, pattern="_vox_count_")){
     return(x)
  }
  
  # Parse feature values to 3 significant numbers
  y <- signif(x, digits=n_digits)
  
  # If feature values lie between 1000 and 10000, use n_digits+1 significant numbers
  # y[(x<10000) & (x>=1000)] <- signif(x[(x<10000) & (x>=1000)], digits=n_digits+1)
  
  return(y)
}


parse_numbers <- function(x, tol, n_sign=3){

  # Limit x to n_sign significant digits
  if(x!=0){
    x_base <- floor(log10(abs(x)))
  } else {
    x_base <- 0
  }
  tol_base <- floor(log10(abs(tol)))
  
  # Create x and tolerance with correct number of significant digits
  if(x_base >= tol_base){
    x      <- signif(x, digits=n_sign)
    if(tol!=0){
      n_dig <- n_sign - (x_base - tol_base)
      if(n_dig > 0){
        tol <- 10^x_base * signif(tol / 10^x_base, digits=n_dig) 
      } else {
        tol <- 10^(x_base - n_sign + 1)
      }
    }
  } else if(tol!=0){
    tol    <- signif(tol, digits=n_sign)
    n_dig  <- n_sign - (tol_base - x_base)
    if(n_dig > 0){
      x      <- 10^tol_base * signif(x / 10^tol_base, digits=n_dig) 
    } else {
      x <- 0
    }
  }

  return(list("value"=x, "tolerance"=tol))
}


parse_numbers_vec <- function(x, tol, tag, n_sign=3){

  # Avoid changing the following tags
  excl_tags <- get_diagnostic_feature_tags()
  
  dt <- data.table::data.table("x"=x, "tol"=tol, "x_base"=0, "y_base"=0, "tag"=tag)
  
  # Determine base for non-zero entries
  dt[x!=0, "x_base":=floor(log10(abs(x)))]
  dt[tol!=0, "tol_base":=floor(log10(abs(tol)))]
  
  # Choose the bases that occur most
  dt[, ":="("x_base"=mode(x_base), "tol_base"=mode(tol_base)), by=tag]
  
  dt[!tag %in% excl_tags & tol==0, "x":=signif(x, digits=n_sign)]
  dt[tag %in% excl_tags & tol==0 & x_base >= n_sign, "x":=signif(x, digits=x_base+1)]
  dt[tag %in% excl_tags & tol==0 & x_base < n_sign, "x":=signif(x, digits=n_sign)]
  
  dt[!tag %in% excl_tags & x_base>=tol_base & tol!=0, "x":=signif(x, digits=n_sign)]
  dt[tag %in% excl_tags &  x_base>=tol_base & tol!=0 & x_base < n_sign, "x":=signif(x, digits=n_sign)]
  dt[tag %in% excl_tags &  x_base>=tol_base & tol!=0 & x_base >= n_sign, "x":=signif(x, digits=x_base+1)]
  
  dt[x_base>=tol_base & tol!=0, "n_dig":=n_sign - (x_base - tol_base)]
  dt[x_base>=tol_base & tol!=0 & n_dig> 0, "tol":=10^x_base * signif(tol / 10^x_base, digits=n_dig) + 10^(x_base - n_sign + 1)]
  dt[x_base>=tol_base & tol!=0 & n_dig<=0, "tol":=10^(x_base - n_sign + 1)]
  
  dt[x_base< tol_base & tol!=0, ":="("tol"=signif(tol, digits=n_sign), "n_dig"=n_sign - (tol_base - x_base))]
  dt[x_base< tol_base & tol!=0 & n_dig> 0, "x":=10^tol_base * signif(x / 10^tol_base, digits=n_dig) ]
  dt[x_base< tol_base & tol!=0 & n_dig<=0, "x":=0]
  
  return(list("value"=dt$x, "tolerance"=dt$tol))
}


check_meshing <- function(x, tag){

  # Remove values for morphological features that are not based on meshing, but do require it
  if(!"morph_volume" %in% tag){
    remove_tags <- tag %in% c("morph_area_mesh", "morph_av", "morph_comp_1", "morph_comp_2", "morph_sph_dispr", "morph_sphericity", "morph_asphericity", "morph_diam",
                              "morph_vol_dens_aabb", "morph_area_dens_aabb", "morph_vol_dens_ombb", "morph_area_dens_ombb", "morph_vol_dens_aee",
                              "morph_area_dens_aee", "morph_vol_dens_mvee", "morph_area_dens_mvee", "morph_vol_dens_conv_hull", "morph_area_dens_conv_hull",
                              "morph_integ_int")
    
    x[remove_tags] <- as.double(NA)
  }
  
  return(x)
}



get_mode <- function(dt){
  # Determine mode of a series of values
  mode        <- dt$value[which.max(dt$n)]
  return(mode)
}


mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


add_tolerance <- function(dt){
  # Check if tolerance column is present (this happens for the parsed old-style data)
  if("tolerance" %in% colnames(dt)){
    return(dt)
  }
  
  # Read datum
  datum <- dt$datum[1]
  
  # Get all tolerance files
  tol_files <- list.files(path="./tolerance", full.names=TRUE, recursive=FALSE, pattern="*.RDS")
  tol_dates <- as.integer(sub('\\..*$', '', basename(tol_files)))
  
  # Find the tolerance date that is smaller or equal to the current datum
  sel_date <- tol_dates[tol_dates <= datum]
  sel_date <- sel_date[which.max(sel_date)]
  sel_tol_file <- tol_files[tol_dates == sel_date]
  
  # Read tolerance file
  dt_tol <- readRDS(sel_tol_file)
  
  # Update morphology columns for the digital phantom data set
  morph_cols <- c("morph_volume", "morph_vol_approx", "morph_area_mesh", "morph_av",
                  "morph_comp_1", "morph_comp_2", "morph_sph_dispr", "morph_sphericity",
                  "morph_asphericity", "morph_com", "morph_diam", "morph_pca_maj_axis",
                  "morph_pca_min_axis", "morph_pca_least_axis", "morph_pca_elongation",
                  "morph_pca_flatness", "morph_vol_dens_aabb", "morph_area_dens_aabb",
                  "morph_vol_dens_ombb", "morph_area_dens_ombb", "morph_vol_dens_aee",
                  "morph_area_dens_aee", "morph_vol_dens_mvee", "morph_area_dens_mvee",
                  "morph_vol_dens_conv_hull", "morph_area_dens_conv_hull", "morph_integ_int",
                  "morph_moran_i", "morph_geary_c")
  
  # Find the mode
  # First remove all non-meshing results
  dt_data <- data.table::copy(dt)
  dt_data[, "value":=check_meshing(value, tag), by=team]
  dt_data <- dt_data[!is.na(value) & data_set=="digital phantom" & tag %in% morph_cols]
  
  # Get the number of different values for each feature and find the most common value
  dt_data <- dt_data[, .(n=.N), by=.(tag, value)]
  dt_data[, "mode":=get_mode(.SD), by=tag]
  dt_data <- dt_data[value==mode]
  dt_data[, "tolerance":=0.005*mode]
  dt_data[, ":="("value"=NULL, "n"=NULL, "mode"=NULL, "data_set"="digital phantom")]
  
  # Rearrange columns
  data.table::setcolorder(dt_data, c("data_set", "tolerance", "tag"))
  
  # Append to tolerance table
  dt_tol <- rbind(dt_data, dt_tol)
  
  # Merge with data
  dt <- merge(x=dt, y=dt_tol, by=c("data_set", "tag"), all.x=TRUE, all.y=FALSE)
  
  # Set na entries
  dt[is.na(tolerance), "tolerance":=0]
  
  # Rearrange columns
  data.table::setcolorder(dt, c("datum", "data_set", "team", "tag", "value", "tolerance"))
  
  return(dt)
}


get_benchmark <- function(dt, aggregate=TRUE){

  # Copy
  dt <- data.table::copy(dt)
  
  # Parse value and tolerance
  dt[, c("value","tolerance"):=parse_numbers_vec(value, tolerance, tag), by=.(datum, data_set)]
  
  # Set the mode
  dt[, "mode":=mode(value), by=.(datum, data_set, tag)]
  
  # Find matching
  dt[, "match":=FALSE]
  dt[between(value, mode-tolerance, mode+tolerance, incbounds=TRUE), "match":=TRUE]

  if(aggregate){
    # Count matching
    dt <- dt[, .(n_matches=sum(match), n_total=.N), by=.(datum, data_set, tag, mode, tolerance)]
  }
  
  return(dt)
}


set_data_sets <- function(dt, by_phase=FALSE){
  
  # Convert data set to factor
  if(!by_phase){
    dt$data_set <- factor(dt$data_set, levels=c("digital phantom", "configuration A", "configuration B",
                                                "configuration C", "configuration D", "configuration E"))
  } else {
    dt$data_set <- factor(dt$data_set, levels=c("digital phantom", "configuration A", "configuration B",
                                                "configuration C", "configuration D", "configuration E"),
                          labels=c("phase I", "phase II", "phase II", "phase II", "phase II", "phase II"))
  }
  
  return(dt)
}


get_diagnostic_feature_tags <- function(){
  return(c("img_dim_x_init_img", "img_dim_y_init_img", "img_dim_z_init_img", "vox_dim_x_init_img", "vox_dim_y_init_img", "vox_dim_z_init_img",
           "mean_int_init_img", "min_int_init_img", "max_int_init_img", "img_dim_x_interp_img", "img_dim_y_interp_img", "img_dim_z_interp_img",
           "vox_dim_x_interp_img", "vox_dim_y_interp_img", "vox_dim_z_interp_img", "mean_int_interp_img", "min_int_interp_img", "max_int_interp_img",
           "int_mask_dim_x_init_roi", "int_mask_dim_y_init_roi", "int_mask_dim_z_init_roi", "int_mask_bb_dim_x_init_roi", "int_mask_bb_dim_y_init_roi",
           "int_mask_bb_dim_z_init_roi", "morph_mask_bb_dim_x_init_roi", "morph_mask_bb_dim_y_init_roi", "morph_mask_bb_dim_z_init_roi",
           "int_mask_vox_count_init_roi", "morph_mask_vox_count_init_roi", "int_mask_mean_int_init_roi", "int_mask_min_int_init_roi",
           "int_mask_max_int_init_roi", "int_mask_dim_x_interp_roi", "int_mask_dim_y_interp_roi", "int_mask_dim_z_interp_roi",
           "int_mask_bb_dim_x_interp_roi", "int_mask_bb_dim_y_interp_roi", "int_mask_bb_dim_z_interp_roi", "morph_mask_bb_dim_x_interp_roi",
           "morph_mask_bb_dim_y_interp_roi", "morph_mask_bb_dim_z_interp_roi", "int_mask_vox_count_interp_roi", "morph_mask_vox_count_interp_roi",
           "int_mask_mean_int_interp_roi", "int_mask_min_int_interp_roi", "int_mask_max_int_interp_roi", "int_mask_dim_x_reseg_roi", "int_mask_dim_y_reseg_roi",
           "int_mask_dim_z_reseg_roi", "int_mask_bb_dim_x_reseg_roi", "int_mask_bb_dim_y_reseg_roi", "int_mask_bb_dim_z_reseg_roi", "morph_mask_bb_dim_x_reseg_roi",
           "morph_mask_bb_dim_y_reseg_roi", "morph_mask_bb_dim_z_reseg_roi", "int_mask_vox_count_reseg_roi", "morph_mask_vox_count_reseg_roi", "int_mask_mean_int_reseg_roi",
           "int_mask_min_int_reseg_roi", "int_mask_max_int_reseg_roi"))
}



get_teams <- function(){
  dt <- fread("./teams.csv")
  
  dt$main_language <- factor(dt$main_language, levels=c("IDL", "Java", "R", "Python", "C++", "MATLAB"))
  dt$code_source   <- factor(dt$code_source, levels=c("in-house", "web-based", "freeware", "open source"))
  
  return(dt)
  
}



add_teams <- function(dt){
  
  dt_team <- get_teams()
  
  dt_out <- merge(dt, dt_team, by="team")
  
  return(dt_out)
}


get_text_size <- function(ggtheme){
  geom_test_size_rel <- 1.0
  geom_text_size     <- ggtheme$text$size
  if(!is.null(ggtheme$axis.text$size)){
    if(class(ggtheme$axis.text$size)=="rel"){
      geom_text_size_rel <- as.numeric(ggtheme$axis.text$size)
    } else {
      geom_text_size <- as.numeric(ggtheme$axis.text$size)
      geom_text_size_rel <- 1.0
    }
  }
  if(!is.null(ggtheme$axis.text.y$size)){
    if(class(ggtheme$axis.text.y$size)=="rel"){
      geom_text_size_rel <- as.numeric(ggtheme$axis.text.y$size)
    } else {
      geom_text_size <- as.numeric(ggtheme$axis.text.y$size)
      geom_test_size_rel <- 1.0
    }
  }
  
  geom_text_size <- geom_text_size * geom_text_size_rel / 2.845276 * 0.7
  
  return(geom_text_size)
}