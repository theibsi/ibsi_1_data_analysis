create_table_2_s1 <- function(dt){
  
  dt <- data.table::copy(dt)
  
  # Set order
  dt_A <- set_data_sets(dt=dt, by_phase=TRUE)
  dt_B <- set_data_sets(dt=dt, by_phase=FALSE)
  
  # Remove "digital phantom" from B
  dt_B <- droplevels(dt_B[data_set!="digital phantom"])
  
  # Join A and B
  dt <- rbind(dt_A, dt_B)
  
  # Add datum_id
  dt <- dt[order(datum, data_set, tag)]
  dt <- dt[, "datum_id":=.GRP, by=datum]
  
  # Remove diagnostic features
  excl_tags <- get_diagnostic_feature_tags()
  dt <- dt[!tag %in% excl_tags]
  
  # Create categories
  dt[between(n_matches, 0, 2), "match_cat":="weak"]
  dt[between(n_matches, 3, 5), "match_cat":="moderate"]
  dt[between(n_matches, 6, 9), "match_cat":="strong"]
  dt[between(n_matches, 10, Inf), "match_cat":="very strong"]
  
  dt$match_cat <- ordered(dt$match_cat, levels=c("weak", "moderate", "strong", "very strong"))
  
  # Add dissent
  dt[, "n_dissent":=n_total-n_matches]
  dt[between(n_dissent, 0, 2), "dissent_cat":="weak"]
  dt[between(n_dissent, 3, 5), "dissent_cat":="moderate"]
  dt[between(n_dissent, 6, 9), "dissent_cat":="strong"]
  dt[between(n_dissent, 10, Inf), "dissent_cat":="very strong"]
  
  dt$dissent_cat <- ordered(dt$dissent_cat, levels=c("weak", "moderate", "strong", "very strong"))

  # Find if consensus is better than dissent
  dt[, "not_exceeds_dissent":= n_matches <= 0.5 * n_total]

  # Aggregate counts for each category
  dt_n_1 <- dt[, .(n=.N, n_not_exceed_dissent=sum(not_exceeds_dissent)), by=.(datum_id, data_set, match_cat)]
  dt_n_1[, "n_total":=sum(n), by=.(datum_id, data_set)]
  
  # Moderate or better
  dt_n_2 <- dt_n_1[match_cat>="moderate", .(match_cat="ge_moderate", n=sum(n), n_not_exceed_dissent=sum(n_not_exceed_dissent), n_total=mean(n_total)), by=.(datum_id, data_set)]

  # Strong or better
  dt_n_3 <- dt_n_1[match_cat>="strong", .(match_cat="ge_strong", n=sum(n), n_not_exceed_dissent=sum(n_not_exceed_dissent), n_total=mean(n_total)), by=.(datum_id, data_set)]
  
  # All
  dt_n_4 <- dt_n_1[, .(match_cat="all", n=sum(n), n_not_exceed_dissent=sum(n_not_exceed_dissent), n_total=mean(n_total)), by=.(datum_id, data_set)]
  
  # Combine
  dt_n <- rbind(dt_n_1, dt_n_2, dt_n_3, dt_n_4)
  
  # Add percentages
  dt_n[, ":="("n_perc"=100*n/n_total, "n_perc_not_exceed_dissent"=100*n_not_exceed_dissent/n), by=.(datum_id, data_set)]
  
  # Add in consent and dissent strings
  dt_n[, ":="("consent"=enc2utf8(paste0(n, " (", format(round(n_perc, 1), nsmall=1, trim=TRUE), ")")),
              "dissent"=enc2utf8(paste0(n_not_exceed_dissent, ifelse(n_not_exceed_dissent==0, " (—)", paste0(" (", format(round(n_perc_not_exceed_dissent, 1), nsmall=1, trim=TRUE), ")")))))]
  
  # Cast wide
  dt_n <- dcast(dt_n, datum_id+data_set+n_total~match_cat, value.var=c("consent", "dissent"), fill=enc2utf8("0 (—)"))
  
  # Select only initial and final data points
  dt_n <- dt_n[datum_id %in% c(1L, 10L, max(datum_id))]
  
  # Drop dissent columns for ge_moderate and ge_strong
  dt_n[, ":="("dissent_ge_moderate"=NULL, "dissent_ge_strong"=NULL, "consent_all"=NULL)]
  
  # Reorder columns
  setcolorder(dt_n, c("datum_id", "data_set", "n_total", "consent_weak", "dissent_weak", "consent_moderate", "dissent_moderate", "consent_strong", "dissent_strong", "consent_very strong",
                      "dissent_very strong", "consent_ge_moderate", "consent_ge_strong", "dissent_all"))

  # Write to file
  fwrite(dt_n, file="./results/table_2.csv", sep=";")
  
  # Export weak or dissenting image biomarkers at the final step
  dt_excl <- dt[datum_id==max(datum_id) & data_set!="phase II" & (match_cat=="weak" | not_exceeds_dissent==TRUE)]
  fwrite(dt_excl, file="./results/table_S1.csv", sep=";")
  
}


create_table_s2_s3 <- function(dt){
  # Unique institutions (S2) and languages (S3)
  
  dt <- data.table::copy(dt)
  
  # Select only most recent time point
  dt <- dt[datum==max(datum) & match==TRUE]
  
  # Remove diagnostic features
  excl_tags <- get_diagnostic_feature_tags()
  dt <- dt[!tag %in% excl_tags]
  
  # Count the number of matches for each feature
  dt[, "n_matches":=sum(match), by=.(data_set, tag)]
  
  # Create categories
  dt[between(n_matches, 0, 2), "match_cat":="weak"]
  dt[between(n_matches, 3, 5), "match_cat":="moderate"]
  dt[between(n_matches, 6, 9), "match_cat":="strong"]
  dt[between(n_matches, 10, Inf), "match_cat":="very strong"]
  
  dt$match_cat <- ordered(dt$match_cat, levels=c("weak", "moderate", "strong", "very strong"))
  
  # Add team data
  dt <- add_teams(dt)
  
  # Count the number of unique main institutions and main languages
  dt_n <- dt[, .(uniq_instit=uniqueN(institution), uniq_language=uniqueN(main_language)), by=.(tag, data_set, match_cat)]
  
  # Count the number of unique institutions per matching category over all data
  dt_n[, .(n=.N), by=.(match_cat, uniq_instit)][order(match_cat, uniq_instit)]
  
  # Write to file
  fwrite(dt_n[, .(n=.N), by=.(match_cat, uniq_instit)][order(match_cat, uniq_instit)], "./results/table_S2.csv", sep=";")
  
  # Count the number of unique languages per matching category over all data
  fwrite(dt_n[, .(n=.N), by=.(match_cat, uniq_language)][order(match_cat, uniq_language)], "./results/table_S3.csv", sep=";")

}

create_table_s4 <- function(dt, dt_name){
  # Table with image biomarker and team information for the digital phantom
  
  dt <- data.table::copy(dt)
  
  # Focus on digital phantom
  dt <- dt[data_set=="digital phantom"]
  
  # Identify earliest data for each biomarker
  dt[, "earliest_date":=min(datum), by="tag"]
  
  # Identify the number of teams at each date
  dt[, "n_teams":=uniqueN(team), by="datum"]
  
  # Identify the number of implementers for each biomarker at each date
  dt[, "n_teams_implemented":=uniqueN(team), by=c("datum", "tag")]
  
  # Identify the number of matching teams
  dt[, "n_teams_match":=sum(match), by=c("datum", "tag")]
  
  # Drop team, value, tolerance, match and mode columns
  dt[, ":="("team"=NULL, "value"=NULL, "tolerance"=NULL, "mode"=NULL, "match"=NULL)]
  
  # Drop all duplicate rows
  dt <- unique(dt)
  
  # Keep only first and final date
  dt <- dt[datum==earliest_date | datum==max(datum)]
  
  # Simplify datum
  dt[datum==earliest_date, "status":="initial"]
  dt[datum==max(datum), "status":="final"]
  dt$status <- factor(dt$status, levels=c("initial", "final"))
  
  # Drop datum, earliest date
  dt[, ":="("datum"=NULL, "earliest_date"=NULL)]
  
  # Merge with name table
  dt <- merge(x=dt, y=dt_name, by="tag", all.x=TRUE, all.y=FALSE)
  
  # Cast columns wide
  dt <- dcast(dt, order+family+image_biomarker~status, value.var=c("n_teams_match", "n_teams_implemented", "n_teams"))
  
  # Order columns
  setorder(dt, "order")
  
  # Reorder table
  setcolorder(dt, c("order", "family", "image_biomarker", "n_teams_match_initial", "n_teams_implemented_initial", "n_teams_initial",
                    "n_teams_match_final", "n_teams_implemented_final", "n_teams_final"))
  
  # Save to file
  fwrite(dt, "./results/table_S4.csv", sep=";")
}


create_benchmark_tables <- function(dt, dt_name){
  
  # Select most recent date
  dt <- dt[datum==max(dt$datum)]

  # Create categories
  dt[between(n_matches, 0, 2), "match_cat":="weak"]
  dt[between(n_matches, 3, 5), "match_cat":="moderate"]
  dt[between(n_matches, 6, 9), "match_cat":="strong"]
  dt[between(n_matches, 10, Inf), "match_cat":="very strong"]
  
  dt$match_cat <- ordered(dt$match_cat, levels=c("weak", "moderate", "strong", "very strong"))
  
  # Remove details tags
  dt[grepl(pattern="*_2D_avg$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_2D_avg"), "subset"="2D_avg")]
  dt[grepl(pattern="*_2D_comb$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_2D_comb"), "subset"="2D_comb")]
  dt[grepl(pattern="*_2_5D_avg$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_2_5D_avg"), "subset"="2_5D_avg")]
  dt[grepl(pattern="*_2_5D_comb$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_2_5D_comb"), "subset"="2_5D_comb")]
  dt[grepl(pattern="*_3D_avg$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_3D_avg"), "subset"="3D_avg")]
  dt[grepl(pattern="*_3D_comb$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_3D_comb"), "subset"="3D_comb")]
  dt[grepl(pattern="*_2D$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_2D"), "subset"="2D")]
  dt[grepl(pattern="*_2_5D$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_2_5D"), "subset"="2_5D")]
  dt[grepl(pattern="*_3D$", x=tag), ":="("ibm"=stringr::str_remove(tag, pattern="_3D"), "subset"="3D")]
  
  dt[is.na(subset), "ibm":=tag]

  # Merge with name table
  dt <- merge(x=dt, y=dt_name, by="tag", all.x=TRUE, all.y=FALSE)
    
  # Split table by ibm and parse information
  lapply(split(dt, by="ibm"), .process_ibm_benchmark_to_table)
  
}

.process_ibm_benchmark_to_table <- function(dt){

  if("2D_avg" %in% dt$subset){
    dt$subset <- factor(dt$subset,
                        levels=c("2D_avg", "2D_comb", "2_5D_avg", "2_5D_comb", "3D_avg", "3D_comb"),
                        labels=c("2D, averaged", "2D, slice-merged", "2.5D, direction-merged",
                                 "2.5D, merged", "3D, averaged", "3D, merged"))
    
    has_subset <- TRUE
    
  } else if("2D" %in% dt$subset) {
    dt$subset <- factor(dt$subset,
                        levels=c("2D", "2_5D", "3D"),
                        labels=c("2D", "2.5D", "3D"))
    
    has_subset <- TRUE
    
  } else {
    
    has_subset <- FALSE
  }
  
  # Reorder
  dt <- dt[order(data_set, subset)]
  
  # Get current image biomarker
  current_ibm <- dt$ibm[1]
  current_ibm_name <- dt$image_biomarker[1]
  
  # Convert first character to lowercase
  if(!grepl(pattern="^Moran|^Geary", x=current_ibm_name)){
    current_ibm_name <- paste(tolower(substr(current_ibm_name, 1, 1)), substr(current_ibm_name, 2, nchar(current_ibm_name)), sep="")
  }

  # Check for percent signs
  if(grepl(pattern="%", x=current_ibm_name, fixed=TRUE)){
    current_ibm_name <- gsub(pattern="%", replacement="\\%", x=current_ibm_name, fixed=TRUE)
  }
  
  # Rename levels in data_set to abbreviations
  dt$data_set <- factor(dt$data_set, levels=c("digital phantom", "configuration A", "configuration B", "configuration C", "configuration D", "configuration E"),
                        labels=c("dig. phantom", "config. A", "config. B", "config. C", "config. D", "config. E"))
  
  # Update weak consensus
  dt[match_cat=="weak" | n_matches <= 0.5 * n_total, ":="("mode"=NA_real_, "tolerance"=NA_real_)]
  
  # Update tolerance
  dt[tolerance==0.0, "tolerance":=NA_real_]
  
  # Drop unnecessary columns
  dt[, ":="("datum"=NULL, "tag"=NULL, "n_matches"=NULL, "n_total"=NULL, "ibm"=NULL, "family"=NULL,
            "image_biomarker"=NULL, "order"=NULL, "digital.phantom"=NULL, "configuration.A"=NULL,
            "configuration.B"=NULL, "configuration.C"=NULL, "configuration.D"=NULL, "configuration.E"=NULL)]

  # Define the caption
  caption_text <- paste0("Benchmark table for the \\textit{", current_ibm_name, "} feature.")
  
  # Extend caption for 
  if(sum(is.na(dt$mode))==1){
    caption_text <- paste0(caption_text, " An unset value (\\textemdash) indicates the lack of a reliable benchmark value.")
  } else if(sum(is.na(dt$mode)) > 1){
    caption_text <- paste0(caption_text, " Unset values (\\textemdash) indicate the lack of reliable benchmark values.")
  }
  
  
  # Define the file name
  file_name <- file.path("./results/benchmark_tables", paste0(current_ibm, ".txt"))
  
  # Define a bolding function
  bold <- function(x) {paste('{\\textbf{',x,'}}', sep ='')}
  
  # Write to file
  con <- file(file_name, open="w")
  
  # Reorder and rename columns
  if(has_subset){
    data.table::setcolorder(dt, c("data_set", "subset", "mode", "tolerance", "match_cat"))
    data.table::setnames(dt, c("data", "aggr. method", "value", "tol.", "consensus"))
    
    # To xtable
    x <- xtable::xtable(dt, display=c("d", "s", "s", "g", "g", "s"), digits=c(3, 1, 1, 3, 3, 1),
                        align=c("l", "c", "c", "c", "c", "c"), caption=caption_text)
    
    # File header for spacing and setting font size
    writeLines("\\vspace{2mm}", con=con)
    writeLines("\\small{", con=con)
    
    # Note that the longtable environment is used, as these tend to be long tables.
    writeLines(xtable::print.xtable(x=x, floating=FALSE, tabular.environment="longtable", booktabs=TRUE,
                                    math.style.negative=TRUE, math.style.exponents=TRUE, include.rownames=FALSE,
                                    NA.string="\\textemdash", sanitize.colnames.function=bold, print.results=FALSE),
               con=con)
    
    writeLines("}", con=con)
    writeLines("\\FloatBarrier", con=con)

  } else {
    dt[, "subset":=NULL]
    
    data.table::setcolorder(dt, c("data_set", "mode", "tolerance", "match_cat"))
    data.table::setnames(dt, c("data", "value", "tol.", "consensus"))
    
    # To xtable
    x <- xtable::xtable(dt, display=c("d", "s", "g", "g", "s"), digits=c(3, 1, 3, 3, 1),
                        align=c("l", "c", "c", "c", "c"), caption=caption_text)
    
    # File header for spacing and setting font size
    writeLines("\\vspace{2mm}", con=con)
    writeLines("\\begin{table}[ht]", con=con)
    writeLines("\\centering", con=con)
    writeLines("\\small", con=con)
    
    # Note that the regular tabular enviroment is used, as these tables are smaller.
    writeLines(xtable::print.xtable(x=x, floating=FALSE, tabular.environment="tabular", booktabs=TRUE,
                                    math.style.negative=TRUE, math.style.exponents=TRUE, include.rownames=FALSE,
                                    NA.string="\\textemdash", sanitize.colnames.function=bold, print.results=FALSE),
               con=con)
    
    writeLines(paste0("\\caption{", caption_text, "}"), con=con)
    writeLines("\\end{table}", con=con)
    writeLines("\\FloatBarrier", con=con)
  }

  # Close file connection
  close(con)

}