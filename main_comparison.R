setwd("W:/WORK/_AlexZ/ibsi_1_data/")

source("parse_excel_table.R")
source("tools.R")
source("reporting.R")

# Get all files in data directory
data_files <- list.files(path="./data", full.names=TRUE, recursive=FALSE, pattern="*.xlsx")

# Parse files
dt_data <- data.table::rbindlist(lapply(data_files, parse_xlsx))

# Drop NA values
dt_data <- dt_data[!is.na(value)]

# Find date
date_str <- format(Sys.Date(), format="%Y%m%d")

# Add date to dt_data
dt_data[, "datum":=as.integer(date_str)]

# Reorder
data.table::setcolorder(dt_data, c("datum", "data_set", "team", "tag", "value"))

# Handle morphological features that depend on meshing being performed
dt_data[,"value":=check_meshing(value, tag), by=.(team, data_set)]
dt_data <- dt_data[!is.na(value)]

# Save a copy of the raw data with the date string
saveRDS(dt_data, paste0("data/rawR/", date_str, ".RDS"))

# Read the naming table
dt_name <- data.table::as.data.table(openxlsx::read.xlsx("naming_tables/20181002.xlsx"))

# Iterate over the unique data_sets (configurations) and generate a list of data tables
data_sets <- c("digital phantom", "configuration A", "configuration B", "configuration C", "configuration D", "configuration E")

# Auto-parse the tables
# If this gives aggregation warnings, please check whether each team only has one data set.
report_list <- lapply(data_sets, function(sel_data_set, dt, dt_name) (get_comparison_table(dt=dt[data_set==sel_data_set], dt_name=dt_name)),
                      dt=dt_data, dt_name=dt_name)
names(report_list) <- data_sets

# Write to xlsx file
openxlsx::write.xlsx(x=report_list, file=paste0("compare_tables/", format(Sys.Date(), format="%Y%m%d"), ".xlsx"), asTable=TRUE, col.names=TRUE, row.names=FALSE,
                     firstActiveRow=2, firstActiveCol=5, borders="surrounding")
