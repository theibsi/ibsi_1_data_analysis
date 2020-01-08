get_comparison_table <-function(dt, dt_name){

  # Round values to 3 or 4 significant vakyes
  dt[,"value":=parse_significant(value, tag), by=tag]

  # Get the number of different values for each feature and find the most common value
  dt_n <- dt[, .(n=.N), by=.(tag, value)]
  dt_n[, "mode":=get_mode(.SD), by=tag]
  dt_n <- dt_n[value==mode]
  dt_n[, "value":=NULL]
  
  # Cast the values as function of the teams
  dt <- data.table::dcast(data=dt, tag ~ team, value.var="value") 
  
  # Merge the wide dt table with the number of values and the mode
  dt <- merge(x=dt_n, y=dt, x.all=FALSE, y.all=TRUE, by="tag")
  
  # Merge the wide dt table with the actual ib names
  dt <- merge(x=dt_name, y=dt, x.all=FALSE, y.all=TRUE, by="tag")
  
  # Order by order
  dt <- dt[order(order)]
  
  # Drop order and tag columns
  dt[, ":="("order"=NULL, "tag"=NULL, "digital.phantom"=NULL, "configuration.A"=NULL, "configuration.B"=NULL, "configuration.C"=NULL,
            "configuration.D"=NULL, "configuration.E"=NULL)]
  
  return(dt)
}


save_benchmarks <- function(dt, dt_name){
  
  dt <- data.table::copy(dt)
  
  # Melt naming table and only keep elements that are flagged for inclusion
  dt_name <- melt(dt_name, id.vars=c("order", "family", "image_biomarker", "tag"), value.name="include", variable.name="data_set")
  dt_name <- dt_name[include==1]
  dt_name[, "include":=NULL]
  
  # Rename factors
  dt_name$data_set <- factor(dt_name$data_set, levels=c("digital.phantom", "configuration.A", "configuration.B", "configuration.C", "configuration.D", "configuration.E"),
                             labels=c("digital phantom", "configuration A", "configuration B", "configuration C", "configuration D", "configuration E"))
  
  # Get recent date
  recent_date <- max(dt$datum)
  
  # Keep only the most recent data with n_match geq 3
  dt <- dt[n_matches>=3 & datum==recent_date]
  
  # Merge dt with dt_name
  dt <- merge(x=dt_name, y=dt, by=c("data_set", "tag"), all.x=TRUE, all.y=FALSE)
  
  # Order to data_set column a factor
  dt <- set_data_sets(dt)
  
  # Add consensus column
  dt[is.na(n_matches), "consensus":="< 3" ]
  dt[between(n_matches, 3, 5), "consensus":="3-5"]
  dt[between(n_matches, 6, 9), "consensus":="6-9"]
  dt[between(n_matches, 10, Inf), "consensus":="\u2265 10"]
  
  # Order rows
  dt <- dt[order(data_set, order)]
  
  # Drop unnecessary columns
  dt[, ":="("datum"=NULL, "n_matches"=NULL, "n_total"=NULL, "order"=NULL)]
  
  # Rename columns
  data.table::setnames(dt, old="mode", new="benchmark_val")
  
  # Reorder columns
  data.table::setcolorder(dt, c("data_set", "family", "image_biomarker", "consensus", "benchmark_val", "tolerance", "tag"))
  
  # Save to spreadsheet
  openxlsx::write.xlsx(x=split(dt, by="data_set"), file=paste0("./results/", recent_date, "benchmark_table.xlsx"))
}