
x <- model
args <- capture_args(ignore = c("x", "to_file", "show_trainable"),
                     force = c("show_layer_names"))
args$model <- x

nodes_df <- model_nodes(model)
if (is.keras_model_sequential(model))
  edges_df <- model_edges_sequential(nodes_df)
else
  edges_df <- model_edges_network(model, nodes_df)

graph <- DiagrammeR::create_graph(nodes_df, edges_df)
graph <- DiagrammeR::set_node_attrs(graph, "fixedsize", FALSE)
graph <- DiagrammeR::set_node_attrs(graph, "nodesep", 2)

coords <- local({
  (igraph::layout_with_sugiyama(DiagrammeR::to_igraph(graph)))[[2]] %>%
    dplyr::as_tibble() %>%
    dplyr::rename(
      x = V1,
      y = V2
    ) %>%
    dplyr::mutate(x = 1.5 * x)
})

graph$nodes_df <- graph$nodes_df %>%
  dplyr::bind_cols(coords)

DiagrammeR::render_graph(graph)
}

model_nodes <- function(x){
  assert_that(is.keras_model(x))
  if (is.keras_model_sequential(x)) {
    model_layers <- x$get_config()$layers
    l_name <- map_chr(model_layers, ~purrr::pluck(., "config", "name"))
  } else {
    model_layers <- x$get_config()$layers
    l_name <- model_layers %>% map_chr("name")
  }
  l_type <- model_layers %>% map_chr("class_name")
  l_unit <- map_chr(model_layers, function(layer) if ("units" %in% names(layer$config)) layer$config$units else NA_integer_)
  shape <- layer$config$batch_input_shape
  l_input_shape <- map_chr(model_layers, function(layer) {
      if (is.null(shape)) {
        return(NA_character_)
      } else {
        shape <- purrr::map_chr(shape, ~ .x %||% NA_integer_)
        return(paste(shape, collapse = "x"))
      }
  })

  l_activation <- model_layers %>%
    map_chr(
      ~(purrr::pluck(., "config", "activation") %||% "")
    )
  
  create_node_df(
    n = length(model_layers),
    name = l_name,
    type = l_type,
    units = l_unit,
    label = glue::glue("{l_name}\n{l_type}\n{l_activation}"),
    shape = "rectangle",
    activation = l_activation
  )
}

is.keras_model <- function(x){
  inherits(x, "keras.engine.training.Model")
}

is.keras_model_sequential <- function(x){
  is.keras_model(x) && inherits(x, "keras.engine.sequential.Sequential")
}

is.keras_model_network <- function(x){
  is.keras_model(x) && !is.keras_model_sequential(x)
}
