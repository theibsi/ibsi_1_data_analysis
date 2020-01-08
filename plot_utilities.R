.arrange_plots <- function(plots_or_grobs, ggtheme, n_cols=NULL, n_rows=NULL, widths=NULL, heights=NULL, byrow=TRUE, 
                                   common_legend=FALSE, common_axis_title_x=FALSE, common_axis_title_y=FALSE,
                                   common_strip_x=FALSE, common_strip_y=FALSE, grob_elements=NULL){
  
  # Partly inspired by egg::ggarange by Auguie Baptiste.
  
  if(is.null(n_cols) & is.null(n_rows)){
    stop("The number of columns and the number of rows cannot both be NULL.")
    
  } else if(is.null(n_rows)){
    n_rows <- ceiling(length(plots_or_grobs) / n_cols)
  } else if(is.null(n_cols)){
    n_cols <- ceiling(length(plots_or_grobs) / n_rows)
  }
  
  # Convert to grobs
  grobs <- .convert_to_grob(plots_or_grobs=plots_or_grobs)
  
  # Determine the number of grobs
  n_grobs <- length(grobs)
  
  # Determine whether the grobs would completely span the n_col * n_row matrix.
  if((n_rows * n_cols) > n_grobs){
    
    # Add empty dummy grobs
    grobs <- c(grobs, rep(list(egg::.dummy_gtable), nrow * ncol - n_grobs))
    
    # Update the number of grobs
    n_grobs <- length(grobs)
  }
  
  # Parse widths.
  if (is.numeric(widths) && !inherits(widths, "unit")) {
    # Check whether widths is numeric and not of the grid::unit class.
    widths <- lapply(widths, grid::unit, "null")
    
  } else if(is.null(widths)){
    # Generate widths with width 1 for each grob.
    widths <- lapply(rep(1, n_grobs), grid::unit, "null")
  }
  
  # Convert to unit.list.
  # TODO: this code may be superfluous.
  if (grid::is.unit(widths)){
    widths <- egg::as.unit.list(widths)
  }
  if (grid::is.unit(widths) && length(widths) == 1) {
    widths <- list(widths)
  }
  
  # Parse heights.
  if (is.numeric(heights) && !inherits(heights, "unit")) {
    # Check whether heights is numeric and not of the grid::unit class. 
    heights <- lapply(heights, grid::unit, "null")
    
  } else if(is.null(heights)) {
    # Generate heights with height 1 for each grob
    heights <- lapply(rep(1, n_grobs), grid::unit, "null")
  }
  
  # Convert to unit.list.
  # TODO: this code may be superfluous.
  if (grid::is.unit(heights)){
    heights <- egg::as.unit.list(heights)
  }
  if (grid::is.unit(heights) && length(heights) == 1) {
    heights <- list(heights)
  }
  
  # Determine breaks.
  n_break_points <- if(byrow) n_rows else n_cols
  
  # Generate splits
  if(n_break_points == 1){
    splits <- rep(1, n_grobs)
    
  } else {
    # Grob sequence
    grob_sequence <- seq_along(grobs) 
    
    # Create splits
    splits <- cut(x=grob_sequence, breaks=n_break_points, labels=seq_len(n_break_points))
    
    # Determine indices for ordering widths and heights
    width_ind <- rep_len(seq_along(widths), length.out=n_grobs)
    height_ind <- rep_len(seq_along(heights), length.out=n_grobs)
    
    # Reorder widths and heights
    widths <- c(matrix(widths[width_ind], ncol=n_break_points, byrow=!byrow))
    heights <- c(matrix(heights[height_ind], ncol=n_break_points, byrow=byrow))
  }
  
  # Create gtables for alignment.
  gtable_frames = mapply(egg::gtable_frame, g=grobs, width=widths, height=heights, SIMPLIFY=FALSE)
  
  # Split the frames according to splits.
  split_gtables <- split(gtable_frames, splits)
  
  # Perform alignment
  if(byrow){
    # Combine by row
    rows <- lapply(split_gtables, function(gtable_row_list) (do.call(gridExtra::gtable_cbind, gtable_row_list)))
    g <- do.call(gridExtra::gtable_rbind, rows)
    
  } else {
    # Combine by column
    columns <- lapply(split_gtables, function(gtable_col_list) (do.call(gridExtra::gtable_rbind, gtable_col_list)))
    g <- do.call(gridExtra::gtable_cbind, columns)
  }
  
  return(g)
}


.convert_to_grob <- function(plots_or_grobs){
  
  # Convert to list if the input is a single grob or 
  if(grid::is.grob(plots_or_grobs) | ggplot2::is.ggplot(plots_or_grobs)){
    plots_or_grobs <- list(plots_or_grobs)
  }
  
  # Initialise list of grobs
  grobs <- list()
  
  for(p in plots_or_grobs){
    if(any(class(p) == "familiar_ggplot")){
      
      # Convert to grob
      g <- ggplot2::ggplotGrob(p)
      
      # Make changes to g according to p$custom_grob.
      if(!is.null(p$custom_grob)){
        if(!is.null(p$custom_grob$heights)){
          
          # Iterate over the elements that need to be updated.
          for(ii in seq_len(length(p$custom_grob$heights$name))){
            # Extract name and height
            name <- p$custom_grob$heights$name[ii]
            height <- p$custom_grob$heights$height[ii]
            
            # Find the row containing the grob cell indicated by the name and
            # update the heights.
            g$heights[g$layout[g$layout$name == name, "t"]] <- height
          }
        }
        
        if(!is.null(p$custom_grob$widths)){
          
          # Iterate over the elements that need to be updated.
          for(ii in seq_len(length(p$custom_grob$widths$name))){
            name <- p$custom_grob$widths$name[ii]
            width <- p$custom_grob$widths$width[ii]
            
            # Find the column containing the grob cell indicated by the name and
            # update the widths.
            g$widths[g$layout[g$layout$name == name, "l"]] <- width
          }
        }
      }
    } else if(any(class(p) == "ggplot")){
      
      # Convert to grob
      g <- ggplot2::ggplotGrob(p)
      
    } else if(any(class(p) == "grob")){
      
      # Assign grob
      g <- p
      
    } else {
      warning(paste0("Could not convert an object of class ", class(p), " to a grob."))
      g <- NULL
    }
    
    grobs <- append(grobs, list(g))
  }
  
  return(grobs)
}