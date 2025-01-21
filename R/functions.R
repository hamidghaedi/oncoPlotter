library(tidyverse)
library(ComplexHeatmap)
library(patchwork)
library(readxl)


# Gene Summary 
getGeneSummary <- function(data, gene_col = "GENE", id_col = "ID", mutation_type_col = "Mutation_type") {
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  
  # Check if specified columns exist in the data frame
  if (!all(c(gene_col, id_col, mutation_type_col) %in% colnames(data))) {
    stop("One or more specified columns do not exist in the data frame.")
  }
  
  # Check for NA values in the specified columns
  if (any(is.na(data[[gene_col]]))) {
    stop(paste("The", gene_col, "column contains NA values."))
  }
  if (any(is.na(data[[id_col]]))) {
    stop(paste("The", id_col, "column contains NA values."))
  }
  if (any(is.na(data[[mutation_type_col]]))) {
    stop(paste("The", mutation_type_col, "column contains NA values."))
  }
  
  # Mutation summary: Count the occurrences of each mutation type for each gene
  mutation_summary <- data %>%
    group_by(across(all_of(c(gene_col, mutation_type_col)))) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(mutation_type_col), values_from = Count, values_fill = list(Count = 0))
  
  # Add the 'Total' column: Count total mutations (sum of all mutation types for each gene)
  mutation_summary <- mutation_summary %>%
    mutate(Total = rowSums(select(., -all_of(gene_col)), na.rm = TRUE))
  
  # Add the 'Mutated_Samples' column: Count unique samples for each gene
  mutated_samples <- data %>%
    group_by(across(all_of(gene_col))) %>%
    summarise(Mutated_Samples = n_distinct(.data[[id_col]]), .groups = "drop")
  
  # Join the mutation summary with the mutated samples count
  result <- mutation_summary %>%
    left_join(mutated_samples, by = gene_col) %>%
    arrange(desc(Mutated_Samples))  # Arrange by 'Mutated_Samples' in descending order
  
  return(result)
}
#
getSampleSummary <- function(data, id_col = "ID", gene_col = "GENE", mutation_type_col = "Mutation_type") {
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyr))
  
  # Check if specified columns exist in the data frame
  if (!all(c(id_col, gene_col, mutation_type_col) %in% colnames(data))) {
    stop("One or more specified columns do not exist in the data frame.")
  }
  
  # Check for NA values in the specified columns
  if (any(is.na(data[[id_col]]))) {
    stop(paste("The", id_col, "column contains NA values."))
  }
  if (any(is.na(data[[gene_col]]))) {
    stop(paste("The", gene_col, "column contains NA values."))
  }
  if (any(is.na(data[[mutation_type_col]]))) {
    stop(paste("The", mutation_type_col, "column contains NA values."))
  }
  
  # Mutation summary: Count the occurrences of each mutation type for each sample (ID)
  mutation_summary <- data %>%
    group_by(across(all_of(c(id_col, mutation_type_col)))) %>%
    summarise(Count = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(mutation_type_col), values_from = Count, values_fill = list(Count = 0))
  
  # Add the 'Total' column: Count total mutations (sum of all mutation types for each sample)
  mutation_summary <- mutation_summary %>%
    mutate(Total = rowSums(select(., -all_of(id_col)), na.rm = TRUE))
  
  # Add the 'Mutated_Genes' column: Count unique genes with mutations for each sample
  mutated_genes <- data %>%
    group_by(across(all_of(id_col))) %>%
    summarise(Mutated_Genes = n_distinct(.data[[gene_col]]), .groups = "drop")
  
  # Join the mutation summary with the mutated genes count
  result <- mutation_summary %>%
    left_join(mutated_genes, by = id_col) %>%
    arrange(desc(Total))  # Arrange by 'Total' in descending order
  
  return(result)
}
#
createMutMat <- function(data, 
                         gene_col = "GENE", 
                         id_col = "ID", 
                         mutation_type_col = "Mutation_type", 
                         top_n = NULL, 
                         sample_order_by = NULL) {
  # 1. Validate input arguments
  if (!all(c(gene_col, id_col, mutation_type_col) %in% colnames(data))) {
    stop("The specified column names must exist in the input data.")
  }
  
  # 2. Generate gene summary
  gene_summary <- getGeneSummary(data,
                                 gene_col = gene_col,
                                 id_col = id_col,
                                 mutation_type_col = mutation_type_col)
  
  # Identify genes with a single mutation and genes mutated in only one patient
  single_mut_gene <- gene_summary$GENE[gene_summary$Total == 1]
  gene_mut_a_pts <- gene_summary$GENE[gene_summary$Mutated_Samples == 1]
  message(sprintf("There are %d genes with one mutation.", length(single_mut_gene)))
  message(sprintf("There are %d genes mutated in one patient.", length(gene_mut_a_pts)))
  
  # 3. Select top genes based on `top_n`
  n_gene <- if (is.null(top_n)) nrow(gene_summary) else min(top_n, nrow(gene_summary))
  top_genes <- gene_summary$GENE[seq_len(n_gene)]
  
  # 4. Create mutation profile matrix
  mut_mat <- data %>%
    rename(GENE = {{gene_col}}, ID = {{id_col}}, Mutation_type = {{mutation_type_col}}) %>%
    group_by(ID, GENE) %>%
    filter(GENE %in% top_genes) %>%
    summarise(Type = paste(sort(unique(Mutation_type)), collapse = ";"), .groups = 'drop') %>%
    mutate(Type = if_else(grepl(";", Type), "MultiHit", Type)) %>%
    pivot_wider(names_from = GENE, values_from = Type, values_fill = "") %>%
    column_to_rownames(var = "ID") %>%
    as.matrix()
  
  # 5. Generate patient summary
  pts_summary <- getSampleSummary(data,
                                  gene_col = gene_col,
                                  id_col = id_col,
                                  mutation_type_col = mutation_type_col)
  pts_no_mut <- setdiff(data[[id_col]], pts_summary$ID)
  if (length(pts_no_mut) > 0) {
    message(sprintf("There are %d samples with no mutation.", length(pts_no_mut)))
  }
  
  # 6. Reorder rows in the matrix based on mutation counts
  temp <- ifelse(mut_mat == "", 0, 1)
  row_sums <- rowSums(temp))
ordered_rows <- names(sort(row_sums, decreasing = TRUE))

# Reorder mutation matrix
t_mut_mat <- mut_mat %>%
  t() %>%
  as.data.frame() %>%
  select(all_of(ordered_rows))

# # 7. Handle sample ordering if sample_order_by is provided
# if (!is.null(sample_order_by)) {
#   # Ensure the columns in sample_order_by exist in the data
#   if (!all(sample_order_by %in% colnames(data))) {
#     stop("Some of the columns specified in sample_order_by do not exist in the data.")
#   }
#   
#   # Create ordered dataframe
#   ord_df <- data %>%
#     select(ID, all_of(sample_order_by)) %>%
#     distinct(ID, .keep_all = TRUE) %>%
#     filter(ID %in% colnames(t_mut_mat))  # Keep only IDs present in t_mut_mat column names
#   
#   # Sort based on sample_order_by columns
#   ex_ordered_ids <- ord_df %>%
#     arrange(across(all_of(sample_order_by))) %>%
#     pull(ID)
#   
#   # Re-order mutation matrix columns based on the ordered IDs
#   t_mut_mat <- t_mut_mat %>%
#     select(all_of(ex_ordered_ids))
# }

return(t_mut_mat)
}


# variant classification color setting 
vcColorSetter <- function(mut_mat, col_pal = NULL) {
  # Predefined colors from maftools and custom additions
  predefined_colors <- c(
    RColorBrewer::brewer.pal(11, name = "Paired"),
    RColorBrewer::brewer.pal(11, name = "Spectral")[1:3],
    'black', 'violet', 'royalblue', '#7b7060', '#535c68', 
    '#1F78B4', '#FF7F00', '#9E0142', '#660066', 'skyblue', '#33A02C', 'black', "#ABDDA4" ,"#66C2A5"
  )
  
  # Define names for predefined colors
  predefined_labels <- c(
    'Nonstop_Mutation', 'Frame_Shift_Del', 'IGR', 'Missense_Mutation', 'Silent', 'Nonsense_Mutation',
    'RNA', 'Splice_Site', 'Intron', 'Frame_Shift_Ins', 'In_Frame_Del', 'ITD', 'In_Frame_Ins',
    'Translation_Start_Site', "Multi_Hit", 'Amp', 'Del', 'Complex_Event', 'pathway',
    'Frameshift', 'Splice site', 'Stopgain', 'Gain', 'Loss', 'Missense', "MultiHit", 
    "Non-frameshift","Non-frameshift deletion"
  )
  
  # Create a named color vector
  predefined_color_map <- setNames(predefined_colors, predefined_labels)
  
  # Extract unique variant classifications from the mutation matrix
  unique_variants <- unique(mut_mat[mut_mat != ""])
  
  # Identify unmatched categories
  unmatched_variants <- setdiff(unique_variants, names(predefined_color_map))
  
  # If there are unmatched categories and no color palette provided, stop with a message
  if (length(unmatched_variants) > 0 && is.null(col_pal)) {
    message("The following variant classifications are missing in predefined categories:")
    print(unmatched_variants)
    message("Available predefined categories:")
    print(sort(names(predefined_color_map)))
    stop("Please provide a named vector of colors for the missing categories using the 'col_pal' argument.")
  }
  
  # Validate provided `col_pal` (if not NULL)
  if (!is.null(col_pal)) {
    if (!is.named(col_pal)) {
      stop("'col_pal' must be a named vector where names correspond to the  variant classifications.")
    }
    # Ensure `col_pal` only includes missing variants
    col_pal <- col_pal[names(col_pal) %in% unmatched_variants]
    # Merge predefined colors with the new palette
    predefined_color_map <- c(predefined_color_map, col_pal)
  }
  
  # Return the complete color mapping
  return(predefined_color_map)
}

#alter_fun generator for oncoprint function
alterFun <- function(mut_mat, background_color = "#CCCCCC", dim = 0.9) {
  # Get the color palette for variant classifications
  color_pallets <- vcColorSetter(mut_mat)
  
  # Define the function
  function(x, y, w, h, mutation_type) {
    # If mutation_type is empty, use the background color
    fill_color <- ifelse(mutation_type == "", background_color, color_pallets[mutation_type])
    
    # Draw the rectangle for each mutation type
    grid.rect(x, y, w * dim, h * dim, 
              gp = gpar(fill = fill_color, col = NA))
  }
}

# color for annotations 
heatmapAnnotColorGenerator <- function(var_data, pal) {
  suppressPackageStartupMessages(library(ggsci))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(RColorBrewer))
  
  # Function to process the large list
  extract_character_vectors <- function(large_list) {
    result <- list()
    for (name in names(large_list)) {
      sublist <- large_list[[name]]
      if (is.list(sublist)) {
        first_nested <- sublist[[1]]
        if (is.character(first_nested)) {
          result[[name]] <- as.character(first_nested)
        }
      } else if (is.character(sublist)) {
        result[[name]] <- as.character(sublist)
      }
    }
    return(result)
  }
  
  # Unique levels and count
  levels <- unique(var_data)
  num_levels <- length(levels)
  
  # Load palettes
  source("./R/palettes.R")
  ggsci_db <- extract_character_vectors(ggsci_db)
  
  # Add RColorBrewer palettes to ggsci_db
  R_palettes <- c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
  get_palette_colors <- function(palette_name) {
    brewer.pal(brewer.pal.info[palette_name, "maxcolors"], palette_name)
  }
  R_palette_colors <- lapply(R_palettes, get_palette_colors)
  names(R_palette_colors) <- R_palettes
  ggsci_db <- c(ggsci_db, R_palette_colors)
  
  # Palette metadata
  paletts_names <- data.frame(
    Name = names(ggsci_db),
    SOURCE = c(
      rep("ggsci", length(ggsci_db) - length(R_palettes)),
      rep("RColorBrewer", length(R_palettes))
    ),
    COLOR_COUNT = sapply(ggsci_db, length)
  )
  
  # Ensure `pal` is a single valid name
  if (!pal %in% paletts_names$Name) {
    cat("ERROR: pal is not a valid palette name. See available palette names and color count: \n")
    print(paletts_names)
    stop("ERROR: Invalid palette name.")
  }
  
  # Check if the palette has enough colors
  palette_row <- paletts_names[paletts_names$Name == pal, ]
  if (num_levels > palette_row$COLOR_COUNT) {
    cat("ERROR: Number of needed colors exceeds available colors in the selected palette. \n")
    print(paletts_names)
    stop("ERROR: Insufficient colors in the selected palette.")
  }
  
  # Select and assign colors
  picked_colors <- ggsci_db[[pal]][1:num_levels]
  names(picked_colors) <- levels
  return(picked_colors)
}

# heatmap column annotation creator:
HeatmapColAnnCreator <- function(df, annot_colors = NULL, palettes = NULL, ...) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  
  # Ensure all columns are character (or factor) as annotations typically map to categories
  if (!all(sapply(df, is.character) | sapply(df, is.factor))) {
    stop("All columns in the dataframe must be of character or factor type.")
  }
  
  # Default palettes if none are provided
  default_palettes <- c(
    "Accent","Paired", "Dark2", "npg", "aaas", "nejm", "lancet", "jama", "bmj", "jco", "ucscgb", "observable", "igv", 
    "locuszoom", "gsea", "material", "simpsons", "futurama", "rickandmorty", "startrek", "tron",
    "Pastel1", "Pastel2", "Set1", "Set2", "Set3", 
    "uchicago", "cosmic", "frontiers", "flatui"
  )
  
  if (is.null(palettes)) {
    palettes <- default_palettes
  }
  
  # Initialize color list
  color_list <- list()
  an_vars <- colnames(df)
  
  # Assign colors if annot_colors is provided
  if (!is.null(annot_colors) && is.list(annot_colors)) {
    for (i in names(annot_colors)) {
      color_list[[i]] <- annot_colors[[i]]
    }
    # Assign colors for the remaining variables not in annot_colors
    rest_an_vars <- an_vars[!an_vars %in% names(color_list)]
    for (i in rest_an_vars) {
      color_list[[i]] <- heatmapAnnotColorGenerator(df[[i]], pal = palettes[which(an_vars == i)])
    }
  } else {
    # Generate colors for all variables if annot_colors is NULL
    for (i in an_vars) {
      color_list[[i]] <- heatmapAnnotColorGenerator(df[[i]], pal = palettes[which(an_vars == i)])
    }
  }
  
  # Debugging: Print the color_list content
  #cat("Debug: color_list content\n")
  #print(color_list)
  
  # Validation: Ensure color_list values are named vectors
  for (key in names(color_list)) {
    if (!is.vector(color_list[[key]]) || is.null(names(color_list[[key]]))) {
      stop(paste("ERROR: color_list entry for", key, "is not a named vector. Ensure it maps levels to colors."))
    }
  }
  
  # Dynamically create HeatmapAnnotation arguments
  annotation_data <- lapply(names(df), function(col_name) df[[col_name]])
  names(annotation_data) <- names(df)
  
  # Create HeatmapAnnotation object, passing additional arguments via ...
  annotation <- do.call(columnAnnotation, c(annotation_data, list(col = color_list), list(...)))
  
  return(annotation)
}

