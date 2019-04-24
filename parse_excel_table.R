parse_xlsx <- function(file_name){

  # Extract team name from file name
  team_name <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file_name))
  team_name <- sub(pattern = "[[:digit:]]*_", replacement="", team_name)
  
  dt_data <- data.table::rbindlist(lapply(seq_len(6), read_sheet, file_name=file_name))
  
  dt_data[, "team":=team_name]
  
  return(dt_data)
}



read_sheet <- function(sheet, file_name){
  
  # Read data from sheet
  dt_data <- data.table::as.data.table(openxlsx::read.xlsx(file_name, sheet=sheet))
  
  # Maintain only relevant columns
  dt_data <- dt_data[,c("data_set", "your_result", "tag")]
  
  # Rename the "your_result" column to "value"
  data.table::setnames(dt_data, "your_result", "value")
  
  # Parse the "value" column to numeric
  if(!is.numeric(dt_data$value)) { dt_data$value <- as.numeric(dt_data$value) }
  
  return(dt_data)
}