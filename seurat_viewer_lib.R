library("dplyr")
library("ggplot2")


# get_markers_1 <- function(markers_file) {
#   markers_yaml <- yaml::read_yaml(markers_file)
#   markers_df <- list()
#   markers_df[['cell_type']] <- c()
#   markers_df[['cell_type_label']] <- c()
#   markers_df[['ensembl']] <- c()
#   markers_df[['symbol']] <- c()
#   
#   for (cell_type in names(markers_yaml)) {
#     cell_type_label <- markers_yaml[[cell_type]][['label']]
#     ensembl_vec <- c()
#     symbol_vec <- c()
#     for (ensembl in names(markers_yaml[[cell_type]][['markers']])) {
#       ensembl_vec <- c(ensembl_vec, ensembl)
#       symbol_vec <- c(symbol_vec, markers_yaml[[cell_type]][['markers']][[ensembl]])
#     }
#     
#     markers_df[['cell_type']] <- c(markers_df[['cell_type']],
#                                    rep(cell_type, length(ensembl_vec)))
#     
#     markers_df[['cell_type_label']] <- c(markers_df[['cell_type_label']],
#                                          rep(cell_type_label, length(ensembl_vec)))
#     
#     markers_df[['ensembl']] <- c(markers_df[['ensembl']], 
#                                  ensembl_vec)
#     
#     markers_df[['symbol']] <- c(markers_df[['symbol']], 
#                                 symbol_vec)
#   }
#   
#   markers_df <- data.frame(markers_df)
#   return(markers_df)
# }


ui.cell_types <- function(input_id, markers_file) {
  df <- get_cell_type_markers(markers_file) %>%
    dplyr::distinct(cell_type, cell_type_label) %>%
    dplyr::arrange(cell_type_label)
  
  choices <- list()
  
  for (i in 1:nrow(df)) {
    val <- df[i, 'cell_type']
    label <- df[i, 'cell_type_label']
    
    choices[[label]] <- val
  }
  
  obj <- selectInput(inputId = input_id, 
                     label = h3("Cell Type Markers"), 
                     choices = choices, 
                     selected = NULL, 
                     multiple = TRUE, 
                     selectize = FALSE,
                     size = 8)
  
  return(obj)
}


ui.marker_plot_types <- function(input_id) {
  obj <- selectInput(inputId = input_id, 
                     label = h3("Plot type"), 
                     choices = list("Violin" = "violin", "Ridge" = "ridge", "Bar" = "bar", "t-SNE" = "tsne"), 
                     selected = "violin")
  
  return(obj)
}


ui.clusters <- function(input_id, sce, size=4) {
  cluster_names <- seurat.get_cluster_names(sce)
  
  obj <- selectInput(inputId = input_id, 
                     label = h3("Clusters"), 
                     choices = cluster_names, 
                     selected = NULL, 
                     multiple = TRUE, 
                     selectize = FALSE,
                     size = size)
  
  return(obj)
}


get_cell_type_markers <- function(markers_file) {
  markers_yaml <- yaml.load_file(markers_file)
  
  markers_df <- list()
  markers_df[['cell_type']] <- c()
  markers_df[['cell_type_label']] <- c()
  markers_df[['symbol']] <- c()
  
  for (cell_type in names(markers_yaml)) {
    marker_symbols <- markers_yaml[[cell_type]][['markers']]
    cell_type_label <- markers_yaml[[cell_type]][['name']]
    
    vals <- rep(cell_type, length(marker_symbols))
    markers_df[['cell_type']] <- c(markers_df[['cell_type']], vals)
    
    vals <- rep(cell_type_label, length(marker_symbols))
    markers_df[['cell_type_label']] <- c(markers_df[['cell_type_label']], vals)
    
    markers_df[['symbol']] <- c(markers_df[['symbol']], marker_symbols)
  }
  
  markers_df <- data.frame(markers_df, stringsAsFactors = FALSE)
  
  return(markers_df)
}


get_cell_type_markers_by_type <- function(markers_file, selected_cell_types) {
  markers <- get_cell_type_markers(markers_file) %>%
    dplyr::filter(cell_type %in% selected_cell_types) %>% .$symbol
  
  markers <- as.character(markers)
  return(markers)
}


