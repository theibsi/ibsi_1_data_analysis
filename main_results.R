setwd("W:/WORK/_AlexZ/ibsi_1_data/")

library(data.table)
library(ggplot2)
library(extrafont)
library(gridExtra)

source("tools.R")
source("plots.R")
source("reporting.R")
source("tables.R", encoding="utf-8")

# Get the naming table
dt_name <- data.table::as.data.table(openxlsx::read.xlsx("naming_tables/20181002.xlsx"))

# Load raw data into a list
# Get all files in data/rawR directory
data_files <- list.files(path="./data/rawR", full.names=TRUE, recursive=FALSE, pattern="*.RDS")
data_list <- lapply(data_files, readRDS)

# Add tolerances if required
data_list <- lapply(data_list, add_tolerance)

# Combine lists
# dt <- data_list[[2]] # TEMPORARY
dt <- data.table::rbindlist(data_list)

# Read name table
dt_name <-  data.table::as.data.table(openxlsx::read.xlsx("naming_tables/20181002.xlsx"))

# Determine benchmark, but do not aggregate
dt_bench <- get_benchmark(dt, aggregate=FALSE)

create_table_s2_s3(dt=dt_bench)

create_table_s4(dt=dt_bench, dt_name=dt_name)

# Figure S1
plot_figure_S1(dt=dt_bench, file_type=c("png", "svg"), font_family="Arial", font_size=9, save_size=c(16,12))

# Figure 4
g4 <- plot_figure_4(dt=dt_bench, file_type=c("png", "svg"), font_family="Arial", font_size=9, save_size=c(16,8))

# Figure 3
plot_figure_3(dt=dt, file_type=c("png", "svg"), font_family="Arial", font_size=9, save_size=c(16,16), extra_g=g4)

# Determine benchmark
dt_bench <- get_benchmark(dt)

# Benchmark export
create_benchmark_tables(dt=dt_bench, dt_name=dt_name)

# Figure 2
plot_figure_2(dt=dt_bench, file_type=c("png", "svg"), font_family="Arial", font_size=9, save_size=c(16, 14))

# Results table
create_table_2_s1(dt=dt_bench)
