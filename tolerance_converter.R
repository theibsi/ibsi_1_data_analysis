source("parse_oncoray.R")

tol_files <- c("./oncoray/tolerance/radiomics_ct_phantom_CT__2018-10-02_ptb configuration A_features.csv",
               "./oncoray/tolerance/radiomics_ct_phantom_CT__2018-10-02_ptb configuration B_features.csv",
               "./oncoray/tolerance/radiomics_ct_phantom_CT__2018-10-02_ptb configuration C_features.csv",
               "./oncoray/tolerance/radiomics_ct_phantom_CT__2018-10-02_ptb configuration D_features.csv",
               "./oncoray/tolerance/radiomics_ct_phantom_CT__2018-10-02_ptb configuration E_features.csv")

# Read files
dt        <- data.table::rbindlist(lapply(tol_files, data.table::fread, dec=".", sep=";"), fill=TRUE, use.names=TRUE)

# Add data_set column
dt[, "data_set":=factor(img_data_settings_id, levels=c("ptb configuration A", "ptb configuration B", "ptb configuration C", "ptb configuration D", "ptb configuration E"),
                        labels=c("configuration A", "configuration B", "configuration C", "configuration D", "configuration E"))]

# Add settings_id column
dt[, "settings_id":=.I, by=data_set]

# Drop unnecessary columns
dt[, ":="("id_subject"=NULL, "id_cohort"=NULL, "img_data_settings_id"=NULL, "img_data_modality"=NULL, "img_data_config"=NULL,
          "img_data_noise_level"=NULL, "img_data_noise_iter"=NULL, "img_data_rotation_angle"=NULL, "img_data_roi_randomise_iter"=NULL, "img_data_roi_adapt_size"=NULL,  
          "img_data_translate_x"=NULL, "img_data_translate_y"=NULL, "img_data_translate_z"=NULL, "img_data_voxel_size"=NULL, "img_data_roi"=NULL)]

# Melt table
dt <- melt(data=dt, id.vars=c("data_set", "settings_id"), na.rm=TRUE, variable.name="tag_mirp", variable.factor=FALSE, value.name="value")

# Find tolerance value
dt <- dt[, .(tolerance=0.05*diff(quantile(value, c(0.25, 0.75)))), by=.(data_set, tag_mirp)]

# Load linking table
dt_link <- data.table::as.data.table(openxlsx::read.xlsx("./oncoray/linking_table/20181002.xlsx"))

# Melt naming table and only keep elements that are flagged for inclusion
dt_link <- melt(dt_link, id.vars=c("order", "family", "image_biomarker", "tag", "own_tag"), value.name="include", variable.name="data_set")
dt_link <- dt_link[include==1]
dt_link[, "include":=NULL]

# Rename factors
dt_link$data_set <- factor(dt_link$data_set, levels=c("digital.phantom", "configuration.A", "configuration.B", "configuration.C", "configuration.D", "configuration.E"),
                           labels=c("digital phantom", "configuration A", "configuration B", "configuration C", "configuration D", "configuration E"))
# Remove digital phantom
dt_link <- dt_link[data_set!="digital phantom"]

# Add own tags
dt[data_set=="configuration A", "own_tag":=replace_tag(tag=tag_mirp, tag_pattern=dt_link[data_set=="configuration A"]$own_tag), by=tag_mirp]
dt[data_set=="configuration B", "own_tag":=replace_tag(tag=tag_mirp, tag_pattern=dt_link[data_set=="configuration B"]$own_tag), by=tag_mirp]
dt[data_set=="configuration C", "own_tag":=replace_tag(tag=tag_mirp, tag_pattern=dt_link[data_set=="configuration C"]$own_tag), by=tag_mirp]
dt[data_set=="configuration D", "own_tag":=replace_tag(tag=tag_mirp, tag_pattern=dt_link[data_set=="configuration D"]$own_tag), by=tag_mirp]
dt[data_set=="configuration E", "own_tag":=replace_tag(tag=tag_mirp, tag_pattern=dt_link[data_set=="configuration E"]$own_tag), by=tag_mirp]

# Merge by own_tag to link with actual tag from dt_link
dt <- merge(dt, dt_link, by=c("own_tag", "data_set"), all.x=FALSE, all.y=TRUE)

# Drop all unused info
dt <- dt[, c("data_set", "tolerance", "tag")]

saveRDS(dt, file="tolerance/20181101.RDS")

