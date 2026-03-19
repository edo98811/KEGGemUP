
#' Functions to modify KEGG pathway graphs by adding nodes, edges, or moving nodes.
#' @param g An igraph object representing the KEGG pathway graph.
#' @param new_vertices Data frame of new vertices to add,
#' with necessary columns: name, label, type, KEGG, 
#' graphics_type, x, y, width, height.
#' @param new_vertices A data frame containing the new 
#' vertices to be added. It must have the following columns: 
#' name, label, type, KEGG, graphics_type, x, y, width, height.
#' @param node_positions Data frame of node positions to move nodes, with necessary columns: name, x, y.
#' @return Modified igraph graph object with added nodes, edges, or moved nodes.
#' @details These functions allow you to modify an existing KEGG pathway graph by adding new nodes or edges, 
#' or by moving existing nodes to new positions. The new vertices and edges must have the necessary 
#' attributes to be properly integrated into the graph. 
#' The functions also include checks to ensure that the modifications
#' are valid and do not introduce inconsistencies into the graph structure.
#' @importFrom igraph add_vertices add_edges delete_edges as_data_frame V vertex_attr
#' @importFrom stats setNames
#' @export
add_vertices_kegg <- function(g, new_vertices) {

  # Define necessary columns for new vertices
  necessary_cols <- c("label", "type", "KEGG", "graphics_type", "x", "y", "width", "height"
  )

  # Check if the graph is an igraph object
  if (!inherits(g, "igraph")) {
    warning("Graph is not an igraph object. Cannot add vertices.")
    return(g)
  }

  # Check if new_vertices is a data frame
  if (!inherits(new_vertices, "data.frame")) {
    warning("new_vertices should be a data frame.")
    return(g)
  }

  # Check if new_vertices is empty
  if (nrow(new_vertices) == 0) {
    warning("No new vertices to add.")
    return(g)
  }

  # Check for missing columns
  missing_cols <- setdiff(necessary_cols, names(new_vertices))
  if (length(missing_cols) > 0) {
    warning(
      paste(
        "The following necessary columns are missing in new_vertices:",
        paste(missing_cols, collapse = ", ")
      )
    )
    return(g)
  }

  # If no new vertices to add, return the graph
  if (nrow(new_vertices) == 0) {
    warning("No new vertices to add after removing duplicates and existing nodes.")
    return(g)
  }

  # Get the number of new vertices to add and create a data frame with default attributes
  n_rows <- nrow(new_vertices)
  vertices_to_add <- data.frame(
    lapply(kegg_vertex_defaults(), rep, each = n_rows)
  )

  # Map new_vertices columns to vertices_to_add, ignoring unrecognized columns
  for (col in names(new_vertices)) {
    if (col %in% names(vertices_to_add)) { 
      vertices_to_add[[col]] <- new_vertices[[col]]
    } else {
      warning(
        paste("Column", col, "is not a recognized vertex attribute. It will be ignored.")
      )
    }
  }

  # Assign unique names to new vertices if not already provided
  max_name <- max(as.numeric(igraph::V(g)$name), na.rm = TRUE)
  vertices_to_add$name <- as.character(seq(max_name + 1, max_name + nrow(vertices_to_add)))

  # Style new vertices based on KEGG graphics_type
  vertices_to_add <- style_vertices_igraph(vertices_to_add)

  # Add new vertices to the graph
  g <- igraph::add_vertices(g, nrow(vertices_to_add))

  # Set vertex attributes for the new vertices
  for (col in names(vertices_to_add)) {
    values <- vertices_to_add[[col]]

    indexes_init <- igraph::vcount(g) - nrow(vertices_to_add) + 1
    indexes_end <- igraph::vcount(g)

    igraph::vertex_attr(g, col, indexes_init:indexes_end) <- values
  }

  return(g)
}

#' Add edges to an igraph graph object.
#' @param g An igraph graph object to which edges will be added.
#' @param new_edges A data frame containing the edges to be added. It must have the following columns: from, to, type.
#' @return An igraph graph object with the new edges added.
#' @details This function adds new edges to an existing igraph graph object. It performs several checks to ensure that the new edges are valid and do not introduce inconsistencies into the graph.
#' The function checks if the graph is an igraph object, if the new_edges data frame has the necessary columns, and if the 'from' and 'to' vertices exist in the graph. It also removes duplicate edges based on 'from' and 'to' pairs before adding them to the graph.
#' @importFrom igraph add_edges delete_edges as_data_frame V vertex_attr
#' @importFrom stats setNames
#' @export
add_edges_kegg <- function(g, new_edges) {
  necessary_cols <- c("from", "to", "type", "type")

  if (!inherits(g, "igraph")) {
    warning("Graph is not an igraph object. Cannot add edges.")
    return(g)
  }

  if (!inherits(new_edges, "data.frame")) {
    warning("new_edges should be a data frame.")
    return(g)
  }

  if (nrow(new_edges) == 0) {
    warning("No new edges to add.")
    return(g)
  }

  missing_cols <- setdiff(necessary_cols, names(new_edges))
  if (length(missing_cols) > 0) {
    warning(
      paste(
        "The following necessary columns are missing in new_edges:",
        paste(missing_cols, collapse = ", ")
      )
    )
    return(g)
  }

  if (any(is.na(new_edges$from) | new_edges$from == "" |
          is.na(new_edges$to) | new_edges$to == "")) {
    warning("Columns 'from' or 'to' contain missing or empty values.")
    return(g)
  }

  nrows <- nrow(new_edges)
  edges_to_add <- data.frame(
    lapply(kegg_edge_defaults(), rep, each = nrows),
    stringsAsFactors = FALSE
  )

  # Map new_edges columns to edges_to_add, ignoring unrecognized columns
  for (col in names(new_edges)) {
    if (col %in% names(edges_to_add)) {
      edges_to_add[[col]] <- new_edges[[col]]
    } else {
      warning(
        paste("Column", col, "is not a recognized edge attribute. It will be ignored.")
      )
    }
  }

  # Remove duplicate edges based on from-to pairs
  edges_to_add <- edges_to_add[
    !duplicated(edges_to_add[, c("from", "to")]), ,
    drop = FALSE
  ]

  # Check that referenced nodes exist
  all_vertex_names <- igraph::V(g)$name
  missing_from <- !edges_to_add$from %in% all_vertex_names
  missing_to <- !edges_to_add$to %in% all_vertex_names

  if (any(missing_from) || any(missing_to)) {
    warning(
      paste(
        "Some edges reference non-existent vertices.",
        "These edges will be skipped."
      )
    )
    edges_to_add <- edges_to_add[!missing_from & !missing_to, , drop = FALSE]
  }

  if (nrow(edges_to_add) == 0) {
    warning("No valid edges to add after validation.")
    return(g)
  }

  # Remove existing edges with same from-to pairs to avoid duplicates
  existing_edges <- igraph::as_data_frame(g, what = "edges")
  if (nrow(existing_edges) > 0) {
    duplicate_mask <- paste(existing_edges$from, existing_edges$to) %in%
                      paste(edges_to_add$from, edges_to_add$to)
    if (any(duplicate_mask)) {
      edges_to_delete <- which(duplicate_mask)
      g <- igraph::delete_edges(g, edges_to_delete)
    }
  }

  # Style new edges based on KEGG type
  edges_to_add <- style_edges_igraph(edges_to_add)
  
  # Add new edges with attributes
  edge_list <- as.vector(t(as.matrix(edges_to_add[, c("from", "to")])))
  attr_list <- as.list(edges_to_add[, !names(edges_to_add) %in% c("from", "to"), drop = FALSE])
  
  # Add edges with attributes
  g <- igraph::add_edges(g, edge_list, attr = attr_list)

  return(g)
}

#' Move nodes in an igraph graph object to new positions.
#' @param g An igraph graph object whose nodes will be moved.
#' @param node_positions A data frame containing the new positions for nodes. 
#' It must have the following columns: name, x, y.
#' @return An igraph graph object with the specified nodes moved to new positions.
#' @details This function updates the x and y coordinates of specified 
#' nodes in an igraph graph object. It performs checks to ensure that 
#' the graph is an igraph object, that the node_positions data frame has the necessary columns, 
#' and that the specified nodes exist in the graph. The function will only update the positions of nodes that 
#' are specified in the node_positions data frame and will ignore any nodes that do not exist in the graph.
#' @importFrom igraph as_data_frame V vertex_attr
#' @importFrom stats setNames
#' @export
move_vertices <- function(g, node_postions) {
  expected_cols <- c("name", "x", "y")

  if (!inherits(g, "igraph")) {
    warning("Graph is not an igraph object. Cannot move nodes.")
    return(g)
  }

  if (!inherits(node_postions, "data.frame")) {
    warning("node_positions should be a data frame.")
    return(g)
  }

  missing_cols <- setdiff(expected_cols, names(node_postions))
  if (length(missing_cols) > 0) {
    warning(
      paste(
        "The following necessary columns are missing in node_positions:",
        paste(missing_cols, collapse = ", ")      )
    )
    return(g)
  }

  if (nrow(node_postions) == 0) {
    warning("No node positions provided to move nodes.")
    return(g)
  }

  if (!is.numeric(node_postions$x) || !is.numeric(node_postions$y)) {
    warning("Columns 'x' and 'y' in node_positions must be numeric.")
    return(g)
  }

  if (any(is.na(node_postions$name) | node_postions$name == "")) {
    warning("Column 'name' in node_positions contains missing or empty values.")
    return(g)
  }

  if (any(!node_postions$name %in% igraph::V(g)$name)) {
    warning("Some node names in node_positions do not exist in the graph. These will be ignored.")
    node_postions <- node_postions[node_postions$name %in% igraph::V(g)$name, , drop = FALSE]
  }

  vertices_df <- as_data_frame(g, what = "vertices")
  vertices_df$x[vertices_df$name %in% node_postions$name] <- node_postions$x
  vertices_df$y[vertices_df$name %in% node_postions$name] <- node_postions$y

  # Update the graph with the new vertex positions
  for (col in names(vertices_df)) {
    igraph::vertex_attr(g, col) <- vertices_df[[col]]
  }

  return(g)
}

