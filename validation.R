library(data.table)
# setwd("W:/WORK/_AlexZ/ibsi_1_data/")

.load_results <- function(file_data_path){
  
  # Read data
  data <- readRDS(file_data_path)
  
  # Reorder column names
  setcolorder(data, neworder=c("team", "id", "modality", "tag", "value"))
  
  data[!is.finite(value), "value":=NA_real_]
  
  return(data)
}



.yeo_johnson <- function(lambda, x){
  # After Yeo, I. K., & Johnson, R. A. (2000). A new family of power
  # transformations to improve normality or symmetry. Biometrika, 87(4),
  # 954-959.
  
  # Copy output
  y <- x
  
  # Determine positive and negative elements of the input vector
  pos_index <- x >= 0 & is.finite(x)
  neg_index <- x < 0 & is.finite(x)
  
  # From original value to transformed value
  if(any(pos_index)){
    if(lambda == 0.0){
      y[pos_index] <- log1p(x[pos_index])
    } else {
      y[pos_index] <- ((x[pos_index] + 1)^lambda - 1) / lambda
    }
  }
  
  if (any(neg_index)){
    if(lambda == 2.0){
      y[neg_index] <- -log1p(-x[neg_index])
    } else {
      y[neg_index] <- -((-x[neg_index] + 1)^(2-lambda) - 1) / (2-lambda)
    }
  }
  
  return(y)
}



.yeo_johnson_loglik <- function(x, y, lambda){
  # Compute log-likelihood of yeo-johnson transform
  
  # Remove non-finite values
  valid_values <- is.finite(y)
  x <- x[valid_values]
  y <- y[valid_values]
  
  # Find the number of samples
  n <- length(x)
  
  # Find the mean of y
  y_mean <- mean(y)
  
  # Compute the log-likelihood
  llf <- (lambda - 1.0) * sum(sign(x) * log1p(abs(x))) - n / 2.0 * log(sum((y - y_mean)^2.0) /n)
  
  return(llf)
}


.get_transform_parameters <- function(x){
  
  # Remove non-finite value from x
  x <- x[is.finite(x)]
  
  # Check if all values are the same
  if(all(x - x[1] == 0.0)){
    return(1.0)
  }
  
  ##### Identify lambda for transformation #####
  # Yeo-Johnson transformations (1.0 is linear)
  lambda <- c(-2.0, -1.0, -0.5, 0.0, 0.33333, 0.5, 1.0, 1.5, 2.0)
  y_list <- lapply(lambda, .yeo_johnson, x=x)
  scores <- sapply(seq_len(length(lambda)), function(ii, x, y_list, lambda){
    return(.yeo_johnson_loglik(x=x, y=y_list[[ii]], lambda=lambda[ii]))
  }, x=x, y_list=y_list, lambda=lambda)
  
  # Select optimal lambda that maximises log-likelihood score
  opt_lambda <- lambda[which.max(scores)]
  
  return(opt_lambda)
}


.compute_icc <- function(data, type="2", yeo_johnson=TRUE){
  
  # We start from the following equation: xij = mu + a_i + b_j + e_ij, with mu
  # the population mean, a_i the rater-dependent change, b_j the
  # subject-dependent change and eij an error with mean 0.
  #
  # * type = "1": ICC 1: We assume that each rated value is cased by a random
  # effect of the rater: i.e. the actual rotation is not associated with an
  # systematic change in value, and there is no interaction between rater and
  # subject. [Shrout 1979]. This means that a_i has mean 0, so that wij = ai +
  # eij. Put differently, raters are randomly chosen from a larger population
  # for each sample.
  #
  # * type = "2": ICC 2: There is a panel of raters, and the entire panel
  # evaluates each sample. However, the panel is assumed to be part of a larger
  # population of raters. Raters are assumed to be systematically biased.
  #
  # * type = "3": ICC 3: There is a panel of raters and the entire panel
  # evaluates each subject. The panel is the entire population and raters are
  # assumed to be systematically biased.
  #
  # Confidence intervals are computed after McGraw and Wong (1996), table 7.
  
  # Find tag
  tag <- data$tag[1]
  modality <- data$modality[1]
  
  # Empty placeholder table
  empty_table <- data.table::data.table("tag"=tag, "modality"=modality, "lambda"=1.0,
                                        "icc"=NA_real_, "icc_low"=NA_real_, "icc_up"=NA_real_)
  
  # Suppress NOTES due to non-standard evaluation in data.table
  value <- mu <- b_j <- a_i <- NULL
  
  # Make a local copy
  data <- data.table::copy(data)

  # Remove entries with non-finite values
  data <- data[is.finite(value)]
  
  # Compute the number of samples and raters
  n_samples <- data.table::uniqueN(data, by="id")
  n_raters <- data.table::uniqueN(data, by="team")
  
  # Perform preliminary checks
  if(n_samples <= 1) return(empty_table)
  if(n_raters <= 1) return(empty_table)
  
  if(yeo_johnson){
  
    # Obtain transformation parameter
    lambda <- .get_transform_parameters(x=data$value)

    # Perform Yeo-Johnson transformation
    data[, "value":=.yeo_johnson(lambda=lambda, x=value)]
    
  } else {
    lambda <- 1.0
  }

  # Calculate each parameter in the equation
  data[, "mu":=mean(value)]
  data[, "b_j":=mean(value) - mu, by=c("id")]
  data[, "a_i":=mean(value) - mu, by=c("team")]
  data[, "e_ij":=value - mu - b_j - a_i]
  
  # Calculate mean squared errors:
  # * mse_r between subjects (b_j)
  # * mse_c between raters (a_i)
  # * mse_e of error (e_ij)
  # * mse_w of residual error (a_i + e_ij)

  # Note that we divide by n_raters * n_samples because we populate b_j and a_i
  # across the entire dataset.
  mse_r <- sum(data$b_j^2) / (n_raters * n_samples)
  mse_c <- sum(data$a_i^2) / (n_raters * n_samples)
  mse_e <- sum(data$e_ij^2) / (n_raters * n_samples)
  mse_w <- sum((data$e_ij + data$a_i)^2) / (n_raters * n_samples)

  if(type=="1"){
  
    # Calculate icc for individual rater and rater panel
    if(mse_w == 0) {
      icc <- icc_ci_low <- icc_ci_up <- icc_panel <- icc_panel_ci_low <- icc_panel_ci_up <- 1
      
    } else {
      
      # Single rater
      icc <- (mse_r - mse_w) / (mse_r + (n_raters - 1) * mse_w)

      # Derivation of the confidence interval after McGraw and Wong, 1996, table 7.
      f_score <- mse_r / mse_w
      thresh_low <- f_score / stats::qf(0.975, n_samples-1, n_samples * (n_raters-1))
      thresh_up <- f_score / stats::qf(0.025, n_samples-1, n_samples * (n_raters-1))
      
      # Calcuate confidence intervals from fisher score
      icc_ci_low <- (thresh_low - 1.0) / (thresh_low + n_raters - 1)
      icc_ci_up <- (thresh_up - 1.0) / (thresh_up + n_raters - 1)
      
      # Multi-rater
      icc_panel <- (mse_r  - mse_w) / mse_r
      icc_panel_ci_low <- 1.0 - 1.0 / thresh_low
      icc_panel_ci_up  <- 1.0 - 1.0 / thresh_up
    }
    
  } else if(type=="2"){
    
    # Calculate icc for individual rater and rater panel
    if(mse_e == 0 & mse_c == 0) {
      icc <- icc_ci_low <- icc_ci_up <- icc_panel <- icc_panel_ci_low <- icc_panel_ci_up <- 1
      
    } else {
      # Single rater
      icc <- (mse_r - mse_e) / (mse_r + (n_raters - 1) * mse_e + n_raters / n_samples * (mse_c - mse_e))

      # Derivation of the confidence interval after McGraw and Wong, 1996, table 7.
      a <- n_raters * icc / (n_samples * (1.0 - icc))
      b <- 1.0 + (n_raters * icc * (n_samples - 1)) / (n_samples * (1.0 - icc))
      
      # Determine denominator degrees of freedom
      v <- (a * mse_c + b * mse_e)^2 / ((a * mse_c)^2 / (n_raters - 1) + (b * mse_e)^2 / ((n_samples - 1) * (n_raters - 1)))

      # Determine confidence intervals
      thresh_low <- stats::qf(0.975, n_samples-1, v)
      thresh_up <- stats::qf(0.025, n_samples-1, v)
      
      # Compute confidence intervals from F distribution.
      icc_ci_low <- n_samples * (mse_r - thresh_low * mse_e) / (thresh_low * (n_raters * mse_c + (n_raters * n_samples - n_raters - n_samples) * mse_e) + n_samples * mse_r)
      icc_ci_up <- n_samples * (mse_r - thresh_up * mse_e) / (thresh_up * (n_raters * mse_c + (n_raters * n_samples - n_raters - n_samples) * mse_e) + n_samples * mse_r)
      
      # Multi-rater
      icc_panel <- (mse_r - mse_e) / (mse_r + (mse_c - mse_e) / n_samples)
      icc_panel_ci_low <- icc_ci_low * n_raters / (1.0 + icc_ci_low * (n_raters-1))
      icc_panel_ci_up <- icc_ci_up * n_raters / (1.0 + icc_ci_up * (n_raters-1))
    }
    
  } else if(type=="3"){
    
    # Calculate icc for individual rater and rater panel
    if(mse_e == 0) {
      icc <- icc_ci_low <- icc_ci_up <- icc_panel <- icc_panel_ci_low <- icc_panel_ci_up <- 1
      
    } else {
      # Single rater
      icc <- (mse_r - mse_e) / (mse_r + (n_raters - 1) * mse_e)
      
      # Derivation of the confidence interval after McGraw and Wong, 1996, table 7.
      f_score <- mse_r / mse_e
      thresh_low <- f_score / stats::qf(0.975, n_samples-1, (n_samples-1) * (n_raters-1))
      thresh_up <- f_score / stats::qf(0.025, n_samples-1, (n_samples-1) * (n_raters-1))
      
      # Calcuate confidence intervals from fisher score
      icc_ci_low <- (thresh_low - 1.0) / (thresh_low + n_raters - 1)
      icc_ci_up <- (thresh_up  - 1.0) / (thresh_up  + n_raters - 1)
      
      # Multi-rater
      icc_panel <- (mse_r - mse_e) / mse_r
      icc_panel_ci_low <- 1.0 - 1.0 / thresh_low
      icc_panel_ci_up  <- 1.0 - 1.0 / thresh_up
    }
  }
  
  return(data.table::data.table("tag"=tag, "modality"=modality, "lambda"=lambda,
                                "icc"=icc, "icc_low"=icc_ci_low, "icc_up"=icc_ci_up))
}

.get_validation_data <- function(){
  # Find results files
  results_files <- list.files(path="./validation", pattern="validation_results", recursive=TRUE, full.names=TRUE)
  
  # Load results
  data <- rbindlist(lapply(results_files, .load_results))
  
  return(data)
}


# Obtain icc-data
# results_files <- c("./validation/cardiff/20190906_validation_results_cardiff.RDS",
#                    "./validation/mcgill/20190902_validation_results_mcgill.RDS",
#                    "./validation/oncoray/20190902_validation_results_oncoray.RDS",
#                    "./validation/pyradiomics/20190906_validation_results_pyradiomics.RDS",
#                    "./validation/umcg_beukinga/20190923_validation_results_umcg_beukinga.RDS",
#                    "./validation/umcg_pfaehler/20191014_validation_results_umcg_pfaehler.RDS")
# 
# # Load results
# data <- rbindlist(lapply(results_files, .load_results))
# 
# icc_data <- rbindlist(lapply(split(data, by=c("tag", "modality")), .compute_icc, type="2", yeo_johnson=TRUE))

# 
# # Strip diagnostic features

