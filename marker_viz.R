library("dplyr")
library("ggplot2")
library("ggridges")
library("viridis")


#' Calcuate the centroids for each cluster (based on the t-SNE coordinates) 
#' for a Seurat object
#'
#' @param tsne_coords A dataframe with t-SNE coordinates for each cell. At
#' a minimum, it must contain the columns 'cell_id', 'tsne_1', and 'tsn_2'. 
#' It can also contain other cell metadata columns, which will be ignored.
#' @param cluster_assignments A dataframe with the cluster assignments for
#' each cell. It must include the columns 'cell_id' and 'cluster'.
get_cluster_coords_for_tsne_plot <- function(tsne_coords, cluster_assignments) {
  df <- inner_join(tsne_coords, cluster_assignments, by = "cell_id")
  
  cluster_centroids <- df %>% dplyr::select(cluster, tsne_1, tsne_2) %>%  
    dplyr::group_by(cluster) %>%
    dplyr::summarize(avg_tsne_1 = mean(tsne_1), 
                     avg_tsne_2 = mean(tsne_2))
  
  cluster_centroids <- ungroup(cluster_centroids)
  
  
  levels <- sort(unique(as.numeric(cluster_centroids$cluster)))
  levels <- as.character(levels)
  labels <- levels
  cluster_centroids$cluster <- factor(cluster_centroids$cluster,
                                      levels = levels, 
                                      labels = labels)
  
  return(cluster_centroids)
}


#' Create a t-SNE plot
#'
#' @param tsne_plot_data A dataframe with with t-SNE coordinates
#' and cell metadata
plot_tsne <- function(tsne_plot_data, cell_attribute=NULL, cell_attribute_label=NULL, 
                      show_facets=FALSE, facets_per_row=2, label_clusters=FALSE, alpha=0.25) {  
  
  # !!!!! Should split this into two functions: one faceted, one not !!!!!
  
  if (show_facets & !is.null(cell_attribute)) {
    p <- ggplot() + geom_point(data = tsne_plot_data, 
                               aes(x = tsne_1, 
                                   y = tsne_2, 
                                   color = cluster), 
                               alpha = alpha)
    
    p <- p + facet_wrap(cell_attribute, ncol = facets_per_row)
  } else {
    # t-SNE scatterplot layer (color by specified column)
    if (!is.null(cell_attribute)) {
      # https://stackoverflow.com/questions/15323269/
      #addressing-x-and-y-in-aes-by-variable-number
      p <- ggplot() + geom_point(data = tsne_plot_data, 
                                 aes(x = tsne_1, 
                                     y = tsne_2, 
                                     color = !!ensym(cell_attribute)), 
                                 alpha = alpha)
    } else {
      p <- ggplot() + geom_point(data = tsne_plot_data, 
                                 aes(x = tsne_1, 
                                     y = tsne_2), 
                                 alpha = alpha)
    }
  }
  
  n <- length(unique(tsne_plot_data[, cell_attribute]))
  
  if (n < 10) {
    p <- p + scale_color_brewer(palette = "Set1")
  } else {
    p <- p + scale_color_viridis_d()
  }
  
  # Cluster labels layer
  if (label_clusters) {
    tsne_coords <- tsne_plot_data %>% dplyr::select(cell_id, tsne_1, tsne_2)
    cluster_assignments <- tsne_plot_data %>% dplyr::select(cell_id, cluster)
    
    # Create the plot data for the cluster labels
    cluster_centroids <- get_cluster_coords_for_tsne_plot(tsne_coords, cluster_assignments)

    # Layer for cluster labels
    p <- p + geom_text_repel(data = cluster_centroids,
                             aes(x = avg_tsne_1,
                                 y = avg_tsne_2,
                                 label = cluster),
                             size = 5)

  }
  
  p <- p + theme_minimal()
  p <- p + xlab("t-SNE 1")
  p <- p + ylab("t-SNE 2")
  
  return(p)
}


plot_markers_violin <- function(norm_counts, markers, clusters=NULL, 
                                cluster_column="cluster", y_intercept=NULL, ncol=2) {
  
  plot_data <- norm_counts %>%
    dplyr::rename(cluster_id = !!cluster_column) %>% 
    dplyr::filter(tolower(gene_id) %in% tolower(markers))
  
  if (!is.null(clusters)) {
    plot_data <- plot_data %>% dplyr::filter(cluster_id %in% as.character(clusters))
  }
  
  p <- ggplot(plot_data, aes(x = gene_id, y = expr_val, fill = gene_id)) +
    geom_violin(scale = "width",
                adjust = 1,
                trim = TRUE,
                alpha = 0.75)
  
  if (!is.null(y_intercept)) {
    p <- p + geom_hline(aes(yintercept = y_intercept),
                        colour = "red",
                        alpha = 0.75)
  }
  
  p <- p + scale_fill_viridis_d() +
    facet_wrap(~ cluster_id, ncol = ncol) +
    theme_minimal() +
    coord_flip()
  
  return(p)
}


plot_markers_ridge <- function(norm_counts, markers, clusters=NULL,
                               
                       cluster_column="cluster", x_intercept=NULL, ncol=2) {
  
  plot_data <- norm_counts %>%
    dplyr::rename(cluster_id = !!cluster_column) %>% 
    dplyr::filter(tolower(gene_id) %in% tolower(markers))
  
  if (!is.null(clusters)) {
    plot_data <- plot_data %>% dplyr::filter(cluster_id %in% as.character(clusters))
  }
  
  p <- ggplot(plot_data, aes(x = expr_val, y = gene_id, fill = gene_id)) +
    geom_density_ridges(alpha = 0.75) +
    #               scale_fill_brewer(palette = "Set1") +
    theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
    scale_fill_viridis_d() +
    facet_wrap(~ cluster_id, ncol = ncol)
  
  if (!is.null(x_intercept)) {
    p <- p + geom_vline(aes(xintercept = x_intercept),
                        colour = "red",
                        alpha = 0.75)
  }
  
  return(p)
}


plot_markers_bar <- function(norm_counts, markers, clusters=NULL, 
                     cluster_column="cluster", y_intercept=NULL, ncol=2) {
  
  plot_data <- norm_counts %>%
    dplyr::rename(cluster_id = !!cluster_column) %>% 
    dplyr::filter(tolower(gene_id) %in% tolower(markers))
  
  if (!is.null(clusters)) {
    plot_data <- plot_data %>% dplyr::filter(cluster_id %in% as.character(clusters))
  }
  
  plot_data <- plot_data %>% group_by(cluster_id, gene_id) %>% 
    summarise(avg_expr_val = mean(expr_val))
  
  p <- ggplot(plot_data, aes(x = gene_id, y = avg_expr_val, fill = gene_id)) + 
    geom_col()
  
  if (!is.null(y_intercept)) {
    p <- p + geom_hline(aes(yintercept = y_intercept),
                        colour = "red",
                        alpha = 0.75)
  }
  
  p <- p + scale_fill_viridis_d() +
    facet_wrap(~ cluster_id, ncol = ncol) +
    theme_minimal() +
    coord_flip()
  
  return(p)
}