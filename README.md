# oncoPlotter
This repository contains a collection of personal functions used to visualize the mutational landscape in oncology samples. These scripts are intended for personal record-keeping and may not be optimized for production use, then use at your own risk!:) 

## Steps Toward Generating a Plot

### 1. Reading the Input File:

The input file should be a dataframe with a column for genes, a column for sample IDs, and a column specifying the mutation/event type.

### 2. Creating a Mutation Matrix Using `createMutMat()`:

The mutation matrix will have genes as rows and samples as columns.

### 3. Reordering the Mutation Matrix Based on Variables of Interest:

Currently, this is a manual step where you can reorder the mutation matrix based on the variables of interest.

### 4. Creating the Patient Annotation Dataframe:

The row order of this dataframe must match the column names of the mutation matrix.

### 5. Creating Oncoplot Annotations:

This step is similar to preparing annotations for `ComplexHeatmap` plots. There are functions to help with the creation of sample/patient annotation objects. Examples:
```R
# Function from the ComplexHeatmap package
annotation_legend_param = list(
    var1 = list(direction = "vertical"), 
    var2 = list(direction = "vertical"),
    var3 = list(direction = "horizontal", ncol = length(unique(pts_annot_df[["var3"]]))),
    var4 = list(direction = "horizontal", ncol = length(unique(pts_annot_df[["var4"]]))),
    var5 = list(direction = "horizontal", ncol = length(unique(pts_annot_df[["var5"]])))
)

# Color setting: semi-manual
v1_color = setNames(c("#A6CEE3","#1F78B4","#B2DF8A"), c("Yes","No","etc"))
v2_color = heatmapAnnotColorGenerator(pts_annot_df[["v2"]], pal = "Dark2")

pts_ann <- HeatmapColAnnCreator(pts_annot_df,
                                  annot_colors = list(v1 = v1_color, v2 = v2_color),
                                  annotation_name_gp = gpar(fontsize = 8),
                                  annotation_legend_param = annotation_legend_param) %>%
    re_size(annotation_height = unit(0.3, "pt")) # Re-size the annotation bars

# Automated color setting: colors are assigned to each level of variables in pts_annot_df
pts_ann <- HeatmapColAnnCreator(pts_annot_df,
                                annotation_name_gp = gpar(fontsize = 8),
                                annotation_legend_param = annotation_legend_param) %>%
                                re_size(annotation_height = unit(0.3, "pt"))

```

### 6 - Plotting Using oncoPrint() from the `ComplexHeatmap` Package:

```
onc <- oncoPrint(
  mut_mat, # Created by `createMutMat()`
  alter_fun = alterFun(mut_mat), # Provided in functions.R
  left_annotation = NULL,  # Apply row annotations
  right_annotation = NULL,
  bottom_annotation = pts_ann, # Created by HeatmapColAnnCreator()
  top_annotation = top_annot, # A ComplexHeatmap annotation object
  col = vcColorSetter(mut_mat), # To color variants in the oncoplot
  heatmap_legend_param = list(
    title = "Alteration",
    legend_direction = "horizontal"  # Set legend direction to horizontal
  ),
  show_pct = FALSE,
  show_column_names = TRUE,
  column_title = column_title,
  column_title_side = "top",
  column_names_gp = gpar(fontsize = 8),
  # row_names_gp = gpar(fontsize = 8, col=row_text_colors),
  row_names_gp = gpar(fontsize = 8),
  # row_names_gp = row_name_styles,  # Apply different styles to row names
  row_names_side = "right",
  row_order = colnames(mut_mat),
  column_order = pts_order # Desired patient order
)

# Plot and save
draw(
  onc,
  heatmap_legend_side = "bottom",           # Place heatmap legend at bottom
  annotation_legend_side = "bottom",        # Place annotation legend at bottom
  merge_legends = TRUE                      # Merge heatmap and annotation legends
)

```
