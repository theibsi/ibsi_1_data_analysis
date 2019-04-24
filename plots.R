plot_benchmarks <- function(dt, by_phase=FALSE, file_type="png", font_size=14, colour_scheme=c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5")){
  
  dt <- data.table::copy(dt)
  
  # Set order of data sets
  dt <- set_data_sets(dt=dt, by_phase=by_phase)
  
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
  
  dt$match_cat <- factor(dt$match_cat, levels=c("weak", "moderate", "strong", "very strong"))
  
  dt_n <- dt[, .(n=.N), by=.(datum_id, data_set, match_cat)]
  dt_n[, "n_perc":=100*n/sum(n), by=.(datum_id, data_set)]

  p <- ggplot(data=dt_n, mapping=aes(x=datum_id, y=n_perc, fill=match_cat))
  p <- p + geom_bar(stat="identity")
  p <- p + cowplot::theme_cowplot(font_family="Roboto", font_size=font_size)
  p <- p + scale_fill_manual(name="consensus", values=colour_scheme)
  p <- p + facet_wrap(~data_set, ncol=3)
  p <- p + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1) )
  p <- p + labs(x="analysis time point", y="features (%)")
  p <- p + theme(panel.background= element_blank(),
                 plot.background = element_blank(),
                 strip.background = element_blank(),
                 strip.placement = "outside")


  
  for(curr_file_type in file_type){
    cowplot::save_plot(filename=paste0("./results/", max(dt$datum),   ifelse(by_phase, "benchmark_phase_plot.", "benchmark_config_plot."), curr_file_type),
                       plot=p, ncol=3, nrow=2, dpi=200)
  }
}


plot_figure_2 <- function(dt, file_type, font_family, font_size, save_size){
  
  colour_scheme <- c("#F21A00", "#EBCC2A", "#78B7C5", "#3B9AB2")
  
  ###### Subfigure A
 
  dt_A <- data.table::copy(dt)
   
  # Set order
  dt_A <- set_data_sets(dt=dt_A, by_phase=TRUE)
  
  # Add datum_id
  dt_A <- dt_A[order(datum, data_set, tag)]
  dt_A <- dt_A[, "datum_id":=.GRP, by=datum]
  
  # Remove diagnostic features
  excl_tags <- get_diagnostic_feature_tags()
  dt_A <- dt_A[!tag %in% excl_tags]
  
  # Create categories
  dt_A[between(n_matches, 0, 2), "match_cat":="weak"]
  dt_A[between(n_matches, 3, 5), "match_cat":="moderate"]
  dt_A[between(n_matches, 6, 9), "match_cat":="strong"]
  dt_A[between(n_matches, 10, Inf), "match_cat":="very strong"]
  
  dt_A$match_cat <- factor(dt_A$match_cat, levels=c("weak", "moderate", "strong", "very strong"))
  
  dt_A_n <- dt_A[, .(n=.N), by=.(datum_id, data_set, match_cat)]
  dt_A_n[, "n_perc":=100*n/sum(n), by=.(datum_id, data_set)]
  
  p_A <- ggplot(data=dt_A_n, mapping=aes(x=datum_id, y=n_perc, fill=match_cat))
  p_A <- p_A + geom_bar(stat="identity")
  p_A <- p_A + cowplot::theme_cowplot(font_family=font_family, font_size=font_size)
  p_A <- p_A + scale_fill_manual(name="consensus", values=colour_scheme)
  p_A <- p_A + facet_wrap(~data_set, ncol=3)
  p_A <- p_A + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1) )
  p_A <- p_A + labs(x="analysis time point", y="image biomarkers (%)")
  p_A <- p_A + theme(panel.background= element_blank(),
                 plot.background = element_blank(),
                 strip.background = element_blank(),
                 strip.placement = "outside",
                 strip.text.x=element_text(margin=margin(b=5)),
                 legend.position="none")
  
  ##### Subfigure B
  dt_B <- data.table::copy(dt)
  
  # Set order
  dt_B <- set_data_sets(dt=dt_B, by_phase=FALSE)
  dt_B <- droplevels(dt_B[data_set!="digital phantom"])
  
  # Add datum_id
  dt_B <- dt_B[order(datum, data_set, tag)]
  dt_B <- dt_B[, "datum_id":=.GRP + 9L, by=datum]
  
  # Remove diagnostic features
  dt_B <- dt_B[!tag %in% excl_tags]
  
  # Create categories
  dt_B[between(n_matches, 0, 2), "match_cat":="weak"]
  dt_B[between(n_matches, 3, 5), "match_cat":="moderate"]
  dt_B[between(n_matches, 6, 9), "match_cat":="strong"]
  dt_B[between(n_matches, 10, Inf), "match_cat":="very strong"]
  
  dt_B$match_cat <- factor(dt_B$match_cat, levels=c("weak", "moderate", "strong", "very strong"))
  
  dt_B_n <- dt_B[, .(n=.N), by=.(datum_id, data_set, match_cat)]
  dt_B_n[, "n_perc":=100*n/sum(n), by=.(datum_id, data_set)]
  
  p_B <- ggplot(data=dt_B_n, mapping=aes(x=datum_id, y=n_perc, fill=match_cat))
  p_B <- p_B + geom_bar(stat="identity")
  p_B <- p_B + cowplot::theme_cowplot(font_family=font_family, font_size=font_size)
  p_B <- p_B + scale_fill_manual(name="consensus", values=colour_scheme)
  p_B <- p_B + facet_wrap(~data_set, ncol=5)
  p_B <- p_B + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1) )
  p_B <- p_B + labs(x="analysis time point", y="image biomarkers (%)")
  p_B <- p_B + theme(panel.background= element_blank(),
                     plot.background = element_blank(),
                     strip.background = element_blank(),
                     strip.placement = "outside",
                     strip.text.x=element_text(margin=margin(b=5)),
                     legend.position="bottom",
                     legend.justification="center")
  
  g <- grid.arrange(p_A, p_B, heights=c(2,1), nrow=2)
  
  for(curr_file_type in file_type){
    cowplot::save_plot(filename=paste0("./results/figure_2.", curr_file_type),
                       plot=g, base_width=save_size[1], base_height=save_size[2], units="cm", dpi=200)
  }
}


plot_figure_3 <- function(dt, file_type, font_family, font_size, save_size, extra_g){

  ggtheme <- cowplot::theme_cowplot(font_family=font_family, font_size=font_size)
  
  ###### Subplot 3A ----------------------
  dt_A <- data.table::copy(dt)
  
  # Set by phase
  dt_A <- set_data_sets(dt=dt_A, by_phase=TRUE)
  
  # Find participation at the different time points
  dt_A <- unique(dt_A[, c("datum", "data_set", "team"), with=FALSE])
  
  # Add participation marker
  dt_A[, "included":=TRUE]
  dt_A <- dcast(dt_A, datum+team ~ data_set, value.var="included", fill=FALSE)
  dt_A[get("phase II")==FALSE, "participation":="phase I"]
  dt_A[get("phase II")==TRUE, "participation":= "phase I & II"]
  dt_A$participation <- factor(dt_A$participation, levels=c("phase I", "phase I & II"))
  
  # Add datum_id
  dt_A <- dt_A[order(datum, participation)]
  dt_A <- dt_A[, "datum_id":=.GRP, by=datum]
  
  # Sum number
  dt_A_n <- dt_A[, .(n=.N), by=.(datum_id, participation)]
  
  # Create plot
  p_A <- ggplot(data=dt_A_n, mapping=aes(x=datum_id, y=n, fill=participation))
  p_A <- p_A + geom_bar(stat="identity")
  p_A <- p_A + ggtheme
  p_A <- p_A + scale_fill_manual(name="participation", values=c("#EBCC2A", "#78B7C5"))
  p_A <- p_A + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1) )
  p_A <- p_A + labs(x="analysis time point", y="number of teams")
  p_A <- p_A + theme(panel.background= element_blank(),
                     plot.background = element_blank())

  
  ###### Subplot 3C ----------------------
  dt_B <- dt_A_n[datum_id==max(datum_id)]
  
  # Add in "retired"
  dt_B <- rbind(dt_B, data.table("datum_id"=25L, "participation"="retired", "n"=1L))
  dt_B$participation <- factor(dt_B$participation, levels=c("retired", "phase I", "phase I & II"))
  
  # Create plot
  p_B <- ggplot(data=dt_B, mapping=aes(x="", y=n, fill=participation))
  p_B <- p_B + geom_bar(stat="identity", width=1)
  p_B <- p_B + coord_polar("y", start=0) + geom_text(aes(label=n), position=position_stack(vjust = 0.5), size=get_text_size(ggtheme=ggtheme))
  p_B <- p_B + scale_fill_manual(values=c("#F21A00", "#EBCC2A", "#78B7C5"))
  p_B <- p_B + labs(x=NULL, y=NULL, fill=NULL, title="participation")
  p_B <- p_B + ggtheme
  p_B <- p_B + theme(axis.line = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.background= element_blank(),
                     plot.background = element_blank())
  
  ###### Subplot 3D ----------------------
  dt_C <- get_teams()
  
  # Parse availability
  dt_C <- dt_C[, .(n=.N), by="code_source"]
  
  # Create plot
  p_C <- ggplot(data=dt_C, mapping=aes(x="", y=n, fill=code_source))
  p_C <- p_C + geom_bar(stat="identity", width=1)
  p_C <- p_C + coord_polar("y", start=0) + geom_text(aes(label=n), position=position_stack(vjust = 0.5), size=get_text_size(ggtheme=ggtheme))
  p_C <- p_C + scale_fill_manual(values=c("#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5"))
  p_C <- p_C + labs(x=NULL, y=NULL, fill=NULL, title="availability")
  p_C <- p_C + ggtheme
  p_C <- p_C + theme(axis.line = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.background= element_blank(),
                     plot.background = element_blank())

  
  
  ##### Subplot 3E ----------------------
  dt_D <- get_teams()
  
  # Parse availability
  dt_D <- dt_D[, .(n=.N), by="main_language"]
  
  # Create plot
  p_D <- ggplot(data=dt_D, mapping=aes(x="", y=n, fill=main_language))
  p_D <- p_D + geom_bar(stat="identity", width=1)
  p_D <- p_D + coord_polar("y", start=0) + geom_text(aes(label=n), position=position_stack(vjust = 0.5), size=get_text_size(ggtheme=ggtheme))
  p_D <- p_D + scale_fill_manual(values=c("#B21300", "#F21A00", "#E1AF00", "#EBCC2A", "#78B7C5", "#3B9AB2"))
  p_D <- p_D + labs(x=NULL, y=NULL, fill=NULL, title="language")
  p_D <- p_D + ggtheme
  p_D <- p_D + theme(axis.line = element_blank(),
                     axis.text = element_blank(),
                     axis.ticks = element_blank(),
                     panel.background= element_blank(),
                     plot.background = element_blank())

  g1 <- egg::ggarrange(p_A, extra_g, ncol=1, draw=FALSE)
  g2 <- egg::ggarrange(p_B, p_C, p_D, ncol=3, draw=FALSE)
  g <- arrangeGrob(grobs=list(g1, g2), nrow=2, draw=FALSE, heights=c(2,1))

  for(curr_file_type in file_type){
    cowplot::save_plot(filename=paste0("./results/figure_3.", curr_file_type),
                       plot=g, base_width=save_size[1], base_height=save_size[2], units="cm", dpi=200)
  }

}


plot_figure_4 <- function(dt, file_type, font_family, font_size, save_size){

  dt <- data.table::copy(dt)
  
  # Limit to phase 1 and the last time point
  dt <- dt[datum==max(datum) & data_set=="digital phantom"]
  
  # Determine which features have less than three matches
  dt[, "n_matches":=sum(match), by="tag"]
  
  # Determine the number of unique image biomarkers
  n_ib <- uniqueN(dt$tag)
  
  # Determine the number of tags in each category for each team
  dt <- dt[, list("no_benchmark"=sum(n_matches < 3),
                  "deviant"=sum(match==FALSE & n_matches >= 3),
                  "match"=sum(match==TRUE & n_matches >= 3)),
           by=c("data_set","team")]
  
  dt[, "not_implemented":=n_ib - no_benchmark - deviant - match]
  
  dt <- melt(dt, id.vars=c("data_set", "team"), variable.name="category", value.name="coverage")
  
  # Convert to percentage
  dt[, "coverage":=100 * coverage/n_ib]
  
  # Category as a factor
  dt$category <- factor(dt$category, levels=c("not_implemented", "no_benchmark", "deviant", "match"),
                        labels=c("not implemented", "no benchmark", "deviating", "matching"))
  
  # Add team information
  dt <- add_teams(dt=dt)
  
  # Order by total number implemented
  team_order <- dt[category=="matching"][order(-coverage)]$team_full
  dt$team_full <- factor(dt$team_full, levels=team_order)
  
  # Create plot
  p <- ggplot(data=dt, mapping=aes(x=team_full, y=coverage, fill=category))
  p <- p + geom_bar(stat="identity")
  p <- p + cowplot::theme_cowplot(font_size=font_size, font_family=font_family)
  p <- p + scale_fill_manual(name="coverage", values=c("#EBEBEB", "#F21A00", "#EBCC2A", "#78B7C5"))
  p <- p + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1) )
  p <- p + labs(x=NULL, y="image biomarkers (%)")
  p <- p + theme(panel.background= element_blank(),
                 plot.background = element_blank())
  
  # for(curr_file_type in file_type){
  #   cowplot::save_plot(filename=paste0("./results/figure_4.", curr_file_type),
  #                      plot=p, base_width=save_size[1], base_height=save_size[2], units="cm", dpi=200)
  # }
  return(p)
}


plot_figure_S1 <- function(dt, file_type, font_family, font_size, save_size){
  # Shows coverage, but also for the configurations.
  
  dt <- data.table::copy(dt)
  
  # Limit to the last time point
  dt <- dt[datum==max(datum)]
  
  # Determine which features have less than three matches
  dt[, "n_matches":=sum(match), by=c("tag", "data_set")]
  dt[, "n_tags":=uniqueN(tag), by=c("data_set")]
  
  # Determine the number of unique image biomarkers
  n_ib <- uniqueN(dt$tag)
  
  # Determine the number of tags in each category for each team
  dt <- dt[, list("no_benchmark" = sum(n_matches < 3),
                  "deviant"      = sum(match==FALSE & n_matches >= 3),
                  "match"        = sum(match==TRUE & n_matches >= 3),
                  "n_tags"       = mean(n_tags)),
           by=c("data_set", "team")]
  
  dt[, "not_implemented":=n_tags - no_benchmark - deviant - match]

  # Add missing data sets
  dt <- rbindlist(lapply(split(dt, by="team"), function(dt_team){
    present_data_sets <- unique(dt_team$data_set)
    missing_data_sets <- levels(dt_team$data_set)[!levels(dt_team$data_set) %in% present_data_sets]
    
    # Append missing data sets
    if(length(missing_data_sets) > 0){
      
      dt_team_missing <- data.table("data_set" = missing_data_sets,
                                    "team" = dt_team$team[1],
                                    "no_benchmark" = 0L,
                                    "deviant" = 0L,
                                    "match" = 0L,
                                    "n_tags" = dt_team$n_tags[1],
                                    "not_implemented" = dt_team$n_tags[1])
      
      dt_team <- rbind(dt_team, dt_team_missing)
      
    } 
      
    return(dt_team)
    
  }))
  
  # Prevent a warning
  dt$not_implemented <- as.integer(dt$not_implemented)
  
  dt <- melt(dt, id.vars=c("data_set", "team", "n_tags"), variable.name="category", value.name="coverage")
  
  # Convert to percentage
  dt[, "coverage":=100 * coverage/n_tags]
  
  # Category as a factor
  dt$category <- factor(dt$category, levels=c("not_implemented", "no_benchmark", "deviant", "match"),
                        labels=c("not implemented", "no benchmark", "deviating", "matching"))
  
  # Add team information
  dt <- add_teams(dt=dt)
  
  # Order by total number implemented
  team_order <- dt[category=="matching" & data_set=="digital phantom"][order(-coverage)]$team_full
  dt$team_full <- factor(dt$team_full, levels=rev(team_order))

  # Create plot
  p <- ggplot(data=dt, mapping=aes(x=team_full, y=coverage, fill=category))
  p <- p + geom_bar(stat="identity")
  p <- p + cowplot::theme_cowplot(font_size=font_size, font_family=font_family)
  p <- p + scale_fill_manual(name="coverage", values=c("#ECEBEB", "#F21A00", "#E1AF00", "#78B7C5"))
  p <- p + theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1) )
  p <- p + labs(x=NULL, y="image biomarkers (%)")
  p <- p + facet_wrap(~data_set, ncol=nlevels(dt$data_set))
  p <- p + coord_flip()
  p <- p + theme(panel.background= element_blank(),
                 plot.background = element_blank(),
                 strip.background = element_blank(),
                 strip.placement = "outside",
                 strip.text.x=element_text(margin=margin(b=5)),
                 legend.position="bottom",
                 legend.justification = "center")

  for(curr_file_type in file_type){
    cowplot::save_plot(filename=paste0("./results/figure_S1.", curr_file_type),
                       plot=p, base_width=save_size[1], base_height=save_size[2], units="cm", dpi=300)
  }
}
