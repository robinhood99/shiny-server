library("assertthat")
library("dplyr")
library("ggplot2")
library("ggrepel")
library("ggridges")
library("ggthemes")
library("Seurat")
library("tidyr")

source(file.path(dirname(parent.frame(2)$ofile), "marker_viz.R"))


#' Get the the cell metadata from a Seurat object.
#'
#' @param sce A Seurat single cell experiment object
seurat.get_cell_metadata <- function(sce) {
  cell_metadata <- sce@meta.data
  cell_metadata$cell_id <- rownames(cell_metadata)
  cell_metadata <- cell_metadata %>% dplyr::select(cell_id, everything())
  
  return(cell_metadata)
}


#' Get the cluster assignments from the Seurat object. This function
#' assumes that there is a 'cluster' column in the Seurat metadata.
#'
#' @param sce A Seurat single cell experiment object
seurat.get_cluster_membership <- function(sce) {
  cell_metadata <- seurat.get_cell_metadata(sce)
  cluster_membership <- cell_metadata %>% dplyr::select(cell_id, cluster)
  
  return(cluster_membership)
}


#' Get the cluster names from the Seurat object. This function
#' assumes that there is a 'cluster' column in the Seurat metadata.
#'
#' @param sce A Seurat single cell experiment object
seurat.get_cluster_names <- function(sce) {
  cell_metadata <- seurat.get_cell_metadata(sce)
  cluster_names <- sort(unique(cell_metadata$cluster))
  
  return(cluster_names)
}


seurat.get_counts <- function(sce, format="dataframe_wide", add_cell_metadata=FALSE, 
                                         cell_metadata_columns=NULL, gene_filter=NULL,
                                         data_type = "normalized") {
  
  if (is.null(gene_filter)) {
    if (data_type == "raw") {
      counts_mat <- as.matrix(sce@raw.data)
    } else if (data_type == "normalized") {
      counts_mat <- as.matrix(sce@data)
    } else if (data_type == "scaled") {
      counts_mat <- as.matrix(sce@scale.data)
    }
  } else {
    idx <- tolower(rownames(sce@raw.data)) %in% tolower(gene_filter)
    # gene_filter <- gene_filter[gene_filter %in% tolower(rownames(sce@raw.data))]
    
    if (data_type == "raw") {
      counts_mat <- as.matrix(sce@raw.data[idx, ])
    } else if (data_type == "normalized") {
      counts_mat <- as.matrix(sce@data[idx, ])
    } else if (data_type == "scaled") {
      counts_mat <- as.matrix(sce@scale.data[idx, ])
    }
  }
  
  if (format == "matrix") {
    the_counts <- counts_mat
  }
  
  if (format == "dataframe_wide") {
    the_counts <- data.frame(counts_mat, stringsAsFactors = FALSE)
  }
  
  if (format == "dataframe_long") {
    the_counts <- data.frame(counts_mat, stringsAsFactors = FALSE)
    the_counts$gene_id <- rownames(the_counts)
    the_counts <- the_counts %>% tidyr::gather(cell_id, expr_val, -gene_id)
    
    if (add_cell_metadata) {
      cell_metadata <- seurat.get_cell_metadata(sce)
      
      if (!is.null(cell_metadata_columns)) {
        cell_metadata <- cell_metadata[, c("cell_id", cell_metadata_columns)]
      }
      
      the_counts <- inner_join(the_counts,
                                cell_metadata,
                                by = "cell_id")
    }
  }
  
  return(the_counts)
}


seurat.get_raw_counts <- function(sce, format="dataframe_wide", add_cell_metadata=FALSE, 
                                         cell_metadata_columns=NULL, gene_filter=NULL) {
  
  raw_counts <- seurat.get_counts(sce = sce, 
                                   format = format, 
                                   add_cell_metadata = add_cell_metadata, 
                                   cell_metadata_columns = cell_metadata_columns, 
                                   gene_filter = gene_filter, 
                                   data_type = "raw")
  
  return(raw_counts)
}


seurat.get_normalized_counts <- function(sce, format="dataframe_wide", add_cell_metadata=FALSE, 
                                         cell_metadata_columns=NULL, gene_filter=NULL) {
  
  norm_counts <- seurat.get_counts(sce = sce, 
                                   format = format, 
                                   add_cell_metadata = add_cell_metadata, 
                                   cell_metadata_columns = cell_metadata_columns, 
                                   gene_filter = gene_filter, 
                                   data_type = "normalized")
  
  return(norm_counts)
}


seurat.get_scaled_counts <- function(sce, format="dataframe_wide", add_cell_metadata=FALSE, 
                                  cell_metadata_columns=NULL, gene_filter=NULL) {
  
  scaled_counts <- seurat.get_counts(sce = sce, 
                                  format = format, 
                                  add_cell_metadata = add_cell_metadata, 
                                  cell_metadata_columns = cell_metadata_columns, 
                                  gene_filter = gene_filter, 
                                  data_type = "scaled")
  
  return(scaled_counts)
}


#' Get the t-SNE coordinates from a Seurat object
#'
#' @param sce A Seurat single cell experiment object
seurat.get_tsne_coords <- function(sce, add_cell_metadata=FALSE, 
                                   cell_metadata_columns=NULL) {
  tsne_coords <- data.frame(sce@dr$tsne@cell.embeddings, 
                            stringsAsFactors = FALSE)
  
  tsne_coords$cell_id <- rownames(tsne_coords)
  tsne_coords <- tsne_coords %>% dplyr::select(cell_id, 
                                               tsne_1 = tSNE_1, 
                                               tsne_2 = tSNE_2)
  
  if (add_cell_metadata) {
    cell_metadata <- seurat.get_cell_metadata(sce)

    if (!is.null(cell_metadata_columns)) {
      cell_metadata <- cell_metadata[, c("cell_id", cell_metadata_columns)]
    }

    tsne_coords <- inner_join(tsne_coords,
                              cell_metadata,
                              by = "cell_id")
  }
  
  return(tsne_coords)
}


#' Create a t-SNE plot from a Seurat object
#'
#' @param sce A Seurat single cell experiment object
seurat.plot_tsne <- function(sce, cell_attribute=NULL, cell_attribute_label=NULL, 
                             show_facets=FALSE, facets_per_row=2, label_clusters=FALSE, alpha=0.25) {  
  
  # Create the plot data for the t-SNE scatterplot
  tsne_coords <- seurat.get_tsne_coords(sce)
  cell_metadata <- seurat.get_cell_metadata(sce)

  tsne_plot_data <- inner_join(tsne_coords, cell_metadata, by = "cell_id")
  
  p <- plot_tsne(tsne_plot_data = tsne_plot_data,
                 cell_attribute = cell_attribute, 
                 cell_attribute_label = cell_attribute_label, 
                 show_facets = show_facets, 
                 facets_per_row = facets_per_row, 
                 label_clusters = label_clusters, 
                 alpha = alpha)
  
  return(p)
}


seurat.plot_cell_attribute_across_clusters <- function(sce, cell_attribute,
                                                       cell_attribute_label=NULL,
                                                       show_facets=FALSE) {
  cell_metadata <- seurat.get_cell_metadata(sce)
  
  cell_metadata <- cell_metadata %>% 
    dplyr::select_(cell_attribute, "cluster")
  
  colnames(cell_metadata)[1] <- "cell_attribute"
  
  plot_data <- cell_metadata %>% dplyr::group_by(cluster, cell_attribute) %>%
    dplyr::summarize(num_cells = n())
  
  #     print(head(plot_data))
  
  cell_counts_per_meta_col <- cell_metadata %>% 
    dplyr::select(cell_attribute) %>% 
    dplyr::group_by(cell_attribute) %>%
    dplyr::summarize(total_cell_count = n())
  
  #     print(head(cell_counts_per_meta_col))
  
  plot_data <- inner_join(plot_data, cell_counts_per_meta_col, by = "cell_attribute") %>%
    dplyr::mutate(perc_cells = num_cells/total_cell_count)
  
  print(head(plot_data))
  
  plot_data <- ungroup(plot_data)
  
  p <- ggplot(plot_data, aes(x = cell_attribute, 
                             y = perc_cells))
  
  if (!show_facets) {
    p <- p + geom_col(aes(fill = cluster))
  } else {
    p <- p + geom_col()
  }
  
  
  n <- length(unique(plot_data$cell_attribute))
  
  if (n < 10) {
    p <- p + scale_fill_brewer(palette = "Set1")
  } else {
    p <- p + scale_fill_viridis_d()
  }
  
  p <- p + scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    coord_flip()
  
  if (show_facets) {
    p <- p + facet_wrap(~ cluster, ncol = 3)
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  if (is.null(cell_attribute_label)) {
    p <- p + xlab(cell_attribute)
  } else {
    p <- p + xlab(cell_attribute_label)
  }
  
  p <- p + ylab("# of cells")
  
  return(p)
}


seurat.plot_markers_violin <- function(sce, markers, data_type = "normalized",
                               clusters=NULL, y_intercept=NULL, ncol=2) {
  
  the_counts <- seurat.get_counts(sce = sce, 
                                   format = "dataframe_long", 
                                   add_cell_metadata = TRUE, 
                                   cell_metadata_columns = NULL, 
                                   gene_filter = markers, 
                                   data_type = data_type)
  
  p <- plot_markers_violin(norm_counts = the_counts, 
                           markers = markers, 
                           clusters = clusters,
                           cluster_column = "cluster",
                           y_intercept = y_intercept, 
                           ncol = ncol)
  
  return(p)
}


seurat.plot_markers_ridge <- function(sce, markers, data_type = "normalized",
                              clusters=NULL, x_intercept=NULL, ncol=2) {
  
  the_counts <- seurat.get_counts(sce = sce, 
                                  format = "dataframe_long", 
                                  add_cell_metadata = TRUE, 
                                  cell_metadata_columns = NULL, 
                                  gene_filter = markers, 
                                  data_type = data_type)
  
  p <- plot_markers_ridge(norm_counts = the_counts, 
                          markers = markers, 
                          clusters = clusters,
                          cluster_column = "cluster",
                          x_intercept = x_intercept, 
                          ncol = ncol)
  
  return(p)
}


seurat.plot_markers_bar <- function(sce, markers, data_type = "normalized",
                            clusters=NULL, y_intercept=NULL, ncol=2) {
  
  the_counts <- seurat.get_counts(sce = sce, 
                                  format = "dataframe_long", 
                                  add_cell_metadata = TRUE, 
                                  cell_metadata_columns = NULL, 
                                  gene_filter = markers, 
                                  data_type = data_type)
  
  p <- plot_markers_bar(norm_counts = the_counts, 
                        markers = markers, 
                        clusters = clusters,
                        cluster_column = "cluster",
                        y_intercept = y_intercept, 
                        ncol = ncol)
  
  return(p)
}


seurat.plot_markers_tsne <- function(sce, markers, data_type = "normalized",
                                    clusters=NULL, y_intercept=NULL, ncol=2) {
  FeaturePlot(object = UpdateSeuratObject(sce),
              features = markers)
}


seurat.find_markers <- function(sce, ident_1, ident_2=NULL, logfc_threshold=log2(1.5), adj_pval_cutoff=0.05) {
  df <- FindMarkers(sce, ident.1 = ident_1, ident.2 = ident_1)
  df$gene_id <- rownames(df)
  df <- df %>% dplyr::select(gene_id, pval, avg_logfc, percent_1, percent_2, pval_adj_bonferroni)
}


#' Get the number of cells per cluster
#'
#' @param sce A Seurat single cell experiment object
seurat.num_cells_per_cluster <- function(sce) {
  df <- seurat.get_cell_metadata(sce) %>% 
    dplyr::group_by(cluster) %>%
    dplyr::summarise(num_cells = n())
  
  return(df)
}


#' Plot the number of cells per cluster
#'
#' @param sce A Seurat single cell experiment object
seurat.plot_cells_per_cluster <- function(sce) {
  df <- seurat.num_cells_per_cluster(sce)
    
  p <- ggplot(data = df, aes(x = reorder(cluster, num_cells), y = num_cells)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    xlab("Cluster") +
    ylab("# of cells")
  
  return(p)
}