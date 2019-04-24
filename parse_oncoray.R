library(data.table)

replace_tag <- function(tag, tag_pattern){
  out_tag <- tag_pattern[sapply(tag_pattern, grepl, x=tag, fixed=TRUE)]
  
  if(length(out_tag)==0){
    out_tag <- NA_character_
  }
  
  if(length(out_tag)>1){
    # The right tag should be equal to or smaller in length than tag
    out_tag_len <- sapply(out_tag, nchar)
    
    # Take the one with the longest matching length
    out_tag <- out_tag[which.max(out_tag_len)]
  }
  return(out_tag)
}


parse_csv <- function(file, dt_link){
  # Read data
  dt <- data.table::fread(file=file)
  
  # Drop unnecessary columns
  dt[, ":="("id_subject"=NULL, "id_cohort"=NULL, "img_data_modality"=NULL, "img_data_config"=NULL, "img_data_roi"=NULL)]
  
  # Wide to long
  dt <- melt(dt, id.vars=c("img_data_settings_id"), value.name="value", variable.name="tag_mirp", variable.factor=FALSE)
  
  # Set the current data set
  curr_data_set <- dt$img_data_settings_id[1]
  dt_link <- dt_link[data_set==curr_data_set]
  
  # Find linking tags and merge
  dt[, "own_tag":=replace_tag(tag_mirp, dt_link$own_tag), by=tag_mirp]
  dt <- merge(x=dt, y=dt_link, by="own_tag", all.x=FALSE, all.y=TRUE)[order(order)]
  
  # Drop unncessary columns
  dt[, ":="("own_tag"=NULL, "img_data_settings_id"=NULL, "tag_mirp"=NULL, "order"=NULL)]
  
  # Reorder columns
  data.table::setcolorder(dt, c("data_set", "family", "image_biomarker", "value"))
  
  return(dt)
}

process_oncoray <- function(){
  # Load linking table
  dt_link <- data.table::as.data.table(openxlsx::read.xlsx("./oncoray/linking_table/20181002.xlsx"))
  
  # Melt naming table and only keep elements that are flagged for inclusion
  dt_link <- melt(dt_link, id.vars=c("order", "family", "image_biomarker", "tag", "own_tag"), value.name="include", variable.name="data_set")
  dt_link <- dt_link[include==1]
  dt_link[, "include":=NULL]
  
  # Rename factors
  dt_link$data_set <- factor(dt_link$data_set, levels=c("digital.phantom", "configuration.A", "configuration.B", "configuration.C", "configuration.D", "configuration.E"),
                             labels=c("digital phantom", "configuration A", "configuration B", "configuration C", "configuration D", "configuration E"))
  
  data_list <- c("./oncoray/digital_phantom_CT__2018-10-02_digital phantom_features.csv",
                 "./oncoray/radiomics_ct_phantom_CT__2018-10-02_configuration A_features.csv",
                 "./oncoray/radiomics_ct_phantom_CT__2018-10-02_configuration B_features.csv",
                 "./oncoray/radiomics_ct_phantom_CT__2018-10-02_configuration C_features.csv",
                 "./oncoray/radiomics_ct_phantom_CT__2018-10-02_configuration D_features.csv",
                 "./oncoray/radiomics_ct_phantom_CT__2018-10-02_configuration E_features.csv")
  
  proc_list <- lapply(data_list, parse_csv, dt_link=dt_link)
  names(proc_list) <- c("digital phantom", "configuration A", "configuration B", "configuration C", "configuration D", "configuration E")
  
  openxlsx::write.xlsx(proc_list, "./oncoray/processed.xlsx")
}