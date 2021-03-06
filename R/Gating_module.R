#' Edit, visualize and show statistics from a gating hierarchy
#' @param id shiny id
#' @import shiny
#' @importFrom shinydashboard box tabBox
#' @importFrom DT DTOutput
#' @export
#' @examples
#' \dontrun{
#' library(shiny)
#' library(shinydashboard)
#' library(flowWorkspace)
#' library(flowCore)
#' library(viridis)
#' library(scales)
#' library(ggplot2)
#' library(ggrepel)
#' library(plotly)
#' library(ggridges)
#' 
#' if (interactive()){
#' 
#'   ui <- dashboardPage(
#'     dashboardHeader(title = "Gating"),
#'     sidebar = dashboardSidebar(disable = TRUE),
#'     body = dashboardBody(
#'       GatingUI("module")
#'     )
#'   )
#' 
#'   server <- function(input, output, session) {
#' 
#'     rval <- reactiveValues()
#'     observe({
#'       utils::data("GvHD", package = "flowCore")
#'       rval$gating_set <- GatingSet(GvHD)
#'     })
#'     res <- callModule(Gating, "module", rval = rval)
#'   }
#' 
#'   shinyApp(ui, server)
#' 
#' }
#' }
GatingUI <- function(id) {
  
  ns <- NS(id)
  
  fluidRow(
    column(width = 4,
           box(title = "Gates", width = NULL, height = NULL, 
               collapsible = TRUE, collapsed = FALSE,
               actionButton(ns("show_gate"), label = "Show defining gate"),
               br(),
               br(),
               box(title = "Select",width = NULL, height = NULL, 
                   collapsible = TRUE, collapsed = TRUE,
                   selectInput(ns("gate_selection"), 
                               label = "Select gate", 
                               choices = NULL, 
                               selected = NULL),
               ),
               box(title = "Create",width = NULL, height = NULL, 
                   collapsible = TRUE, collapsed = TRUE,
                   actionButton(ns("get_brush_coords"), "Get path from brush"),
                   #checkboxInput(ns("convex_hull"), "Get convex hull from points", value = FALSE),
                   br(),
                   br(),
                   textInput(ns("gate_name"), label = "Enter gate name", value = ""),
                   actionButton(ns("create_gate"), "create gate"),
                   actionButton(ns("reset_gate"), "reset gate")
               ),
               box(title = "Delete",width = NULL, height = NULL, 
                   collapsible = TRUE, collapsed = TRUE,
                   selectInput(ns("gate_to_delete"), 
                               label = "Select gate", 
                               choices = NULL, 
                               selected = NULL),
                   actionButton(ns("delete_gate"), "Delete")
               ),
               box(title = "Rename", width = NULL, height = NULL, 
                   collapsible = TRUE, collapsed = TRUE,
                   selectInput(ns("gate_to_rename"), 
                               label = "Select gate", 
                               choices = NULL, 
                               selected = NULL),
                   textInput(ns("new_name"), label = "Enter new name", value = ""),
                   actionButton(ns("rename_gate"), "rename gate")
               ),
               box(title = "Copy", width = NULL, height = NULL, 
                   collapsible = TRUE, collapsed = TRUE,
                   selectInput(ns("gate_to_copy"), 
                               label = "Select gate to copy", 
                               choices = NULL,
                               selected = NULL),
                   selectInput(ns("new_parent"), 
                               label = "Select parent gate", 
                               choices = NULL,
                               selected = NULL),
                   checkboxInput(ns("copy_children"), "Copy children gates", TRUE),
                   actionButton(ns("copy_gate"), "copy gate")
               )
                  
                   
           ),
           box(width = NULL, height = NULL, title = "Parameters", 
               collapsible = TRUE, collapsed = FALSE,
               plotCytoUI(id = ns("plot_module"))
           ),
           
           box(title = "Message_gate",
               width = NULL, height = NULL,
               verbatimTextOutput(ns("message_gate"))
           )
    ),
    column(width = 8,
           tabBox(title = "",
                  width = NULL, height = NULL,
                  tabPanel(title = "Plot",
                           simpleDisplayUI(ns("simple_display_module"))
                  )
           ),
           tabBox(title = "Hierarchy",
                  width = NULL, height = NULL,
                  tabPanel("Tree",
                           simpleDisplayUI(ns("simple_display_module_tree")),
                           box(title = "Plot Options", width = NULL, collapsible = TRUE, collapsed = TRUE,
                               uiOutput(ns("tree_ui_options"))
                           )
                  ),
                  tabPanel(title = "Gates",
                           simpleDisplayUI(ns("simple_display_module_2"))
                  ),
                  tabPanel(title = "Stats",
                           br(),
                           downloadButton(ns("download_data"), label = "Download Stats"),
                           br(),
                           br(),
                           div(style = 'overflow-x: scroll', DT::DTOutput(ns("pop_stats"))),
                           br() 
                  ),
                  tabPanel(title = "Export",
                           br(),
                           downloadButton(ns("export_gh"), label = "Export Gating Hierarchy"),
                           br(),
                           br()
                  ),
                  tabPanel(title = "Import", 
                           fileInput(inputId = ns("import_gh"),
                                     label = "Import Gating Hierarchy",
                                     multiple = FALSE
                           ),
                           uiOutput(ns("import_options")),
                           box(title = "Transform", width = NULL, collapsible = TRUE, collapsed = TRUE,
                               textInput(ns("pattern"), label = "Pattern", value  = "Comp-"),
                               textInput(ns("replacement"), label = "Replacement", value  = ""),
                               numericInput(ns("time_step"), label = "Time step", value = 1),
                               actionButton(ns("transform_gates"), label = "Apply transformation")
                           ),
                           box(title = "Preview", width = NULL, collapsible = TRUE, collapsed = TRUE,
                               verbatimTextOutput(ns("import_gh_summary"))
                           ),
                           actionButton(ns("apply_gh"), label = "Import gates")
                  )
                  
           )
    )
    
  )
  
}


#' Gating module server function
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param rval A reactive values object
#' @return The updated reactiveValues object \code{rval}
#' @import shiny
#' @importFrom flowWorkspace gs_pop_add gs_pop_remove setNode recompute gs_get_pop_paths
#' @importFrom flowCore polygonGate
#' @importFrom graph addEdge
#' @importFrom Rgraphviz renderGraph layoutGraph
#' @importFrom methods new
#' @importFrom grDevices chull
#' @importFrom DT datatable renderDT
#' @importFrom dplyr rename
#' @importFrom utils write.table
#' @export
#' @rdname GatingUI
Gating <- function(input, output, session, rval) {
  
  
  plot_params <- reactiveValues() # parameters controlling the main plot
  plot_params_gh <- reactiveValues() # parameters controlling the gating hierarchy plot
  gate <- reactiveValues(x = NULL, y = NULL) # polygon gate represented on plot
  brush_coords <- reactiveValues(x = NULL, y = NULL)
  display_params <- reactiveValues(nrow = 2, width = 200, height = 200)
  rval_mod <- reactiveValues()
  
  observe({
    rval$update_gs <- 0
  })
  
  ### Call modules ###############################################################################
  res <- callModule(plotCyto, "plot_module", 
                    rval=rval,
                    plot_params = plot_params,
                    simple_plot = TRUE, 
                    show_gates = TRUE,
                    use_ggcyto = TRUE,
                    polygon_gate = gate)
  
  res_display <- callModule(simpleDisplay, "simple_display_module", 
                            plot_list = res$plot, 
                            params = reactiveValues(width = 350, height = 350))
  
  plot_all_gates <- callModule(plotGatingHierarchy, "plot_hierarchy_module",
                               rval = rval,
                               plot_params = plot_params_gh)

  callModule(simpleDisplay, "simple_display_module_2",
             plot_list = plot_all_gates,
             params = display_params)
  
  res_tree <- callModule(simpleDisplay, "simple_display_module_tree", 
             plot_list = graph,
             params = reactiveValues(width = 500, height = 500))
  
  ### Get parameters from GatingSet ################################################################
  
  choices <- reactive({
    rval$update_gs
    validate(need(class(rval$gating_set) == "GatingSet", "input is not a GatingSet"))
    get_parameters_gs(rval$gating_set)
  })

  ### Passing plot parameters to the gating hierarchy plot #######################################
  observe({
      for(var in intersect(names(res$params), c("sample",
                                              "plot_type",
                                              "auto_focus",
                                              "legend.position",
                                              "theme",
                                              "bins",
                                              "size", 
                                              "alpha", 
                                              "option",
                                              "color_var",
                                              "facet_var",
                                              "show_outliers")) ){
      plot_params_gh[[var]] <- res$params[[var]]
    }
    display_params$title <- paste(res$params$sample, collapse = " + ")
  })
  
  observe({
    gates <- choices()$subset
    updateSelectInput(session, "gate_to_delete", choices = setdiff(gates, "root"))
    updateSelectInput(session, "gate_to_rename", choices = setdiff(gates, "root"))
    updateSelectInput(session, "gate_to_copy", choices = setdiff(gates, "root"))
    updateSelectInput(session, "new_parent", choices = gates)
    updateSelectInput(session, "gate_selection", choices = gates)
  })
  
  ### Observe environments for gating ############################################################

  observeEvent(input$gate_selection, {
    #update plot_params
    for(var in names(plot_params)){
      plot_params[[var]] <- res$params[[var]]
    }
    plot_params$subset <- input$gate_selection
  })
  
  observe({
    updateSelectInput(session, "gate_selection", selected = res$params$subset)
  })
  
  ### update polygon coordinates upon mouse click on main plot ####################################
  observeEvent(res_display$params$plot_click, {
    if(is.null(res_display$params$plot_brush)) {
      xvar <- res$params$xvar
      yvar <- res$params$yvar
      
      x_coord <- res_display$params$plot_click$x
      if(xvar %in% names(rval$gating_set@transformation)){
        x_coord <- rval$gating_set@transformation[[xvar]]$inverse(x_coord)
      }
      gate$x <- c(gate$x, x_coord)
      
      y_coord <- res_display$params$plot_click$y
      if(yvar %in% names(rval$gating_set@transformation)){
        y_coord <- rval$gating_set@transformation[[yvar]]$inverse(y_coord)
      }
      gate$y <- c(gate$y, y_coord)
  
      #idx <- grDevices::chull(gate$x, gate$y)
      #gate$x <- gate$x[idx]
      #gate$y <- gate$y[idx]
    }
  })
  
  ### update polygon coordinates upon mouse brush on main plot ####################################
  observeEvent(res_display$params$plot_brush, {
    brush <- res_display$params$plot_brush
    
    if (!is.null(brush)) {
      
      xvar <- res$params$xvar
      yvar <- res$params$yvar
      
      x_coord <- c(brush$xmin, brush$xmax, brush$xmax, brush$xmin)
      if(xvar %in% names(rval$gating_set@transformation)){
        x_coord <- rval$gating_set@transformation[[xvar]]$inverse(x_coord)
      }
      brush_coords$x <- x_coord
      
      y_coord <- c(brush$ymin, brush$ymin, brush$ymax, brush$ymax)
      if(yvar %in% names(rval$gating_set@transformation)){
        y_coord <- rval$gating_set@transformation[[yvar]]$inverse(y_coord)
      }
      brush_coords$y <- y_coord
      
    }
    
  })
  
  observeEvent(input$get_brush_coords, {
    gate$x <- brush_coords$x
    gate$y <- brush_coords$y
  })
  
  ### reset polygon upon mouse double click on main plot ##########################################
  observeEvent(res_display$params$plot_dblclick, {
    gate$x <- NULL
    gate$y <- NULL
  })
  
  ### reset polygon using action button ###########################################################
  observeEvent(input$reset_gate, {
    gate$x <- NULL
    gate$y <- NULL
    rval$gate <- NULL
  })
  
  ### reset polygon when plot x axis is modified ##################################################
  observeEvent(res$params$xvar, {
    gate$x <- NULL
    gate$y <- NULL
  })
  
  ### reset polygon when plot y axis is modified ##################################################
  observeEvent(res$params$yvar, {
    gate$x <- NULL
    gate$y <- NULL
  })
  
  ### Create gate from polygon, update rval$gating_set ############################################

  observeEvent(input$create_gate, {
    
    if(input$gate_name %in% basename(flowWorkspace::gs_get_pop_paths(rval$gating_set))){
      showModal(modalDialog(
        title = "Error",
        "Gate name already exists! Please choose another name.",
        easyClose = TRUE,
        footer = NULL
      ))
    }else{
      if(!is.null(gate$x)){
        polygon <- data.frame(x = gate$x, y = gate$y)
        
        #hpts <- grDevices::chull(polygon)
        #polygon <- polygon[hpts, ]
        polygon <- as.matrix(polygon)

        var_names <- c(res$params$xvar, res$params$yvar)
        names(var_names) <- NULL
        colnames(polygon) <- var_names
        print(polygon)
        poly_gate <- flowCore::polygonGate(.gate = polygon, filterId=input$gate_name)
        rval$gate <- poly_gate
          
        flowWorkspace::gs_pop_add(rval$gating_set, poly_gate, parent = res$params$subset)
        flowWorkspace::recompute(rval$gating_set)
        rval$update_gs <- rval$update_gs + 1
        
        gate$x <- NULL
        gate$y <- NULL

        if(res$params$subset != "root"){
          gate_name <- paste(res$params$subset, "/", input$gate_name, sep = "")
        }else{
          gate_name <- paste("/", input$gate_name, sep = "")
        }
        #update plot_params
        for(var in names(plot_params)){
          plot_params[[var]] <- res$params[[var]]
        }
        plot_params$subset <- gate_name
      }
    }
    
  })
  
  ### Delete selected gate, update rval$gating_set ################################################
  observeEvent(input$delete_gate, {
    if(input$gate_to_delete != "root"){
      
      idx_gh <- which( flowWorkspace::gs_get_pop_paths(rval$gating_set) == input$gate_to_delete )
      target_gate <- flowWorkspace::gs_get_pop_paths(rval$gating_set)[idx_gh]

      flowWorkspace::gs_pop_remove(gs = rval$gating_set, node = target_gate)
      flowWorkspace::recompute(rval$gating_set)
      rval$update_gs <- rval$update_gs + 1
      #update plot_params
      for(var in names(plot_params)){
        plot_params[[var]] <- res$params[[var]]
      }
      plot_params$subset <- "root"
      
    }
    
  })
  
  
  ### Rename selected gate, update rval$gating_set ################################################
  observeEvent(input$rename_gate, {
    
    if(input$gate_to_rename != "root"){
      
      idx_gh <- which( flowWorkspace::gs_get_pop_paths(rval$gating_set) == input$gate_to_rename )
      target_gate <- flowWorkspace::gs_get_pop_paths(rval$gating_set)[idx_gh]

      flowWorkspace::setNode(rval$gating_set, target_gate, input$new_name)
      rval$update_gs <- rval$update_gs + 1
      
    }
    
  })

  ### Copy selected gate, update rval$gating_set ################################################
  observeEvent(input$copy_gate, {
    
    
    gs <- rval$gating_set
    gs <- try(copy_gate(gs = gs, 
                    name = input$gate_to_copy, 
                    parent = input$new_parent, 
                    copy_children_gates = input$copy_children), silent = TRUE)
                                 
                    
    if(class(gs) == "try-error"){
      showModal(modalDialog(
        title = "Error",
        print(gs),
        easyClose = TRUE,
        footer = NULL
      ))
    }else{
      rval$gating_set <- gs
      rval$update_gs <- rval$update_gs + 1
    }
    

  })
  
  ### Update plot parameters to show defining gate ################################################
  observeEvent(input$show_gate, {
    
    if(res$params$subset != "root"){
      
      gates <- get_gates_from_gs(gs = rval$gating_set)
      rval$gate <- gates[[res$params$subset]]$gate
      
      #update plot_params
      for(var in names(plot_params)){
        plot_params[[var]] <- res$params[[var]]
      }
      plot_params$subset <- gates[[res$params$subset]]$parent
      
      if(class(rval$gate[[1]]) != "booleanFilter"){
        gate_params <- names(rval$gate[[1]]@parameters)
        
        if(length(gate_params) > 0){
          plot_params$xvar <- gate_params[1]
        }
        if(length(gate_params) > 1){
          plot_params$yvar <- gate_params[2]
        }
      }
      
      
      gate$x <- NULL
      gate$y <- NULL

    }

  })
  
  
  ### Output messages #############################################################################
  
  output$message_gate <- renderPrint({
    print(gate$x)
  })
  
  ### Plot gating tree ############################################################################
  
  gate_list <- reactive({
    rval$update_gs
    validate(need(class(rval$gating_set)=="GatingSet", "No GatingSet available"))
    validate(need(setdiff(choices()$subset, "root"), "No gates to display"))
    gates <- get_gates_from_gs(rval$gating_set)
    
    # if(!input$show_all_subsets){
    #   idx_cluster <- grep("^cluster[0-9]+", basename(names(gates)))
    #   idx_bins <- grep("^bin[0-9]+", basename(names(gates)))
    #   idx_hide <- union(idx_cluster, idx_bins)
    #   if(length(idx_hide)>0){
    #     gates <- gates[-idx_hide]
    #   }
    # }
    
    gates
  })
  
  # observe({
  #   plot_params_gh$selected_subsets <- names(gate_list())
  # })
  

  graph <- reactive({

    gates <- gate_list()
    rankdir <- "LR"
    if(!is.null(input$horizontal_tree)){
      if(!input$horizontal_tree){
        rankdir <- NULL
      }
    }

    fontsize <- ifelse(is.null(input$fontsize), 10, input$fontsize)
    
    p <- plot_tree(gates,
                   fontsize = fontsize,
                   rankdir = rankdir,
                   shape = ifelse(is.null(input$shape), "ellipse", input$shape),
                   fixedsize = ifelse(is.null(input$fixedsize), FALSE, input$fixedsize))

    return(p)
  })
  
  output$tree_ui_options <- renderUI({
  
      ns <- session$ns
      x <- list()
      x[["fontsize"]] <- numericInput(ns("fontsize"), "fontsize", value = 10)
      x[["horizontal_tree"]] <- checkboxInput(ns("horizontal_tree"), label = "horizontal layout", TRUE)
      x[["shape"]] <- selectInput(ns("shape"), "node shape", 
                                  choices = c("ellipse", "circle", "rectangle"), selected = "ellipse")
      x[["fixedsize"]] <- checkboxInput(ns("fixedsize"), label = "fixed node size", FALSE)

      tagList(x)

  })
  
  ### Population statistics ######################################################################
  
  pop_stats <- reactive({
    validate(need(class(rval$gating_set)=="GatingSet", "No GatingSet available"))
    
    spill <- rval$gating_set@compensation
    if(!is.null(rval$apply_comp)){
      if(!rval$apply_comp){
        spill <- NULL
      }
    }
    
    df <- getPopStatsPlus(rval$gating_set, spill = spill)
    df <- df[df$name %in% res$params$sample, ]
    df[['% parent']] <- sprintf("%.1f", df$Count / df$ParentCount * 100)
    df <- df[, c("name", "Population", "Parent", "% parent", "Count", "ParentCount")] 
    df <- dplyr::rename(df, subset = "Population")
    #df <- df[df$subset %in% plot_params_gh$selected_subsets, ]
    df
  })
  
  output$pop_stats <- DT::renderDT({
    DT::datatable(pop_stats(), rownames = FALSE)
  })
  
  output$download_data <- downloadHandler(
    filename = "pop_stats.txt",
    content = function(file) {
      utils::write.table(pop_stats(), file = file, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  )
  

  ### Import/Export Gating Hierarchy ##############################################################
  
  output$import_options <- renderUI({
    ns <- session$ns
    x <- list()
    file_path <- input$import_gh$datapath
    if(!is.null(file_path)){
      if(file_ext(file_path) %in% c("xml")){
        choices <- get_templates_from_ws_diva(file_path)
        x[[1]] <- selectInput(ns("template"), "Template", choices = choices, selected = choices[1])
      }else if(file_ext(file_path) %in% c("wsp")){
        choices <- get_groups_from_ws(file_path)
        x[[1]] <- selectInput(ns("group"), "Group", choices = choices, selected = choices[1])
      }
    }
    tagList(x)
  })
  
  output$export_gh <- downloadHandler(
    filename = "gates.rda",
    content = function(file) {
      gates <- get_gates_from_gs(rval$gating_set)
      save(gates, file=file)
    }
  )
  
  observeEvent(c(input$import_gh, input$template, input$group), {
    validate(
      need(input$import_gh$datapath, "Please select a file")
    )
    file_path <- input$import_gh$datapath
    if(file_ext(file_path) %in% c("wsp") ){
      rval_mod$gating_hierarchy <- get_gates_from_ws(
        ws_path = file_path,
        group = input$group)
    }else if(file_ext(file_path) %in% c("xml") ){
      rval_mod$gating_hierarchy <- get_gates_from_ws_diva(
        ws_path = file_path,
        template = input$template)
    }else if(file_ext(file_path) %in% c("rda", "Rda") ){
      res_name <- load(file_path)
      res <- get(res_name)
      if(class(res) == "list"){
        if("gate" %in% names(res[[1]])){
          rval_mod$gating_hierarchy <- res
        }
      }
    }
  })
  
  output$import_gh_summary <- renderPrint({
    print(paste("Number of gates imported :",length(rval_mod$gating_hierarchy)))
    print(rval_mod$gating_hierarchy)
  })
  
  observe({
    validate(need(class(rval$gating_set)=="GatingSet", "No GatingSet available"))
    ff <- rval$gating_set@data[[1]]
    time_step <- as.numeric(description(ff)[["$TIMESTEP"]])
    updateNumericInput(session, "time_step", value = time_step)
  })
  
  observeEvent(input$transform_gates, {
    validate(need(length(rval_mod$gating_hierarchy)>0, "No gating hierarchy to apply"))
    new_gates <- transform_gates(gates = rval_mod$gating_hierarchy, 
                                 transformation = NULL,
                                 pattern = input$pattern, 
                                 replacement = input$replacement,
                                 time_step = input$time_step)
    rval_mod$gating_hierarchy <- new_gates
  })
  
  observeEvent(input$apply_gh, {
    validate(need(class(rval$gating_set)=="GatingSet", "No GatingSet available"))
    validate(need(length(rval_mod$gating_hierarchy)>0, "No gates to import"))
    
    # old_gates <- gs_pop_get_children(obj = rval$gating_set[[1]], y = "root")
    # for(gate in old_gates){
    #     flowWorkspace::gs_pop_remove(gs = rval$gating_set, node = gate)
    # }
    
    new_gates <- rval_mod$gating_hierarchy
    res <- try(add_gates_flowCore(gs = rval$gating_set, gates = new_gates), silent = TRUE)
    if(class(res) == "try-error"){
      showModal(modalDialog(
        title = "Error",
        print(res),
        easyClose = TRUE,
        footer = NULL
      ))
    }
    rval$update_gs <- rval$update_gs + 1

  })
    
  return(rval)
  
  
}


### Tests #########################################################################################
# 
# library(shiny)
# library(shinydashboard)
# library(flowWorkspace)
# library(flowCore)
# library(viridis)
# library(scales)
# library(ggplot2)
# library(ggrepel)
# library(plotly)
# library(ggridges)
# 
# if (interactive()){
# 
#   ui <- dashboardPage(
#     dashboardHeader(title = "Gating"),
#     sidebar = dashboardSidebar(disable = TRUE),
#     body = dashboardBody(
#       GatingUI("module")
#     )
#   )
# 
#   server <- function(input, output, session) {
# 
#     rval <- reactiveValues()
# 
#     observe({
#       #dataDir <- system.file("extdata",package="flowWorkspaceData")
#       #gs <- load_gs(list.files(dataDir, pattern = "gs_bcell_auto",full = TRUE))
#        fs <- read.ncdfFlowSet(files = "../flowR_utils/demo-data/JL04BMVLG-Valentin/Tumor_T_001_012.fcs")
#        #gates <- gates <- get_gates_from_ws_diva("../flowR_utils/demo-data/JL04BMVLG-Valentin/JL04BMVLG.xml", template = "Gating")
#        gates <- get_gates_from_ws("../flowR_utils/demo-data/JL04BMVLG-Valentin/Tumor-testFlowR.wsp", group = "Tumor")
#        gates <- transform_gates(gates, pattern = "Comp-", replacement = "")
#        gs <- flowWorkspace::GatingSet(fs)
#        
#        add_gates_flowCore(gs, gates = gates)
#       transfo <- lapply(colnames(gs), logicle_trans)
#       names(transfo) <- colnames(gs)
#       gs@transformation <- transfo
#       spill <- gs@data[[1]]@description[["SPILL"]]
#       comp <- lapply(pData(gs)$name, function(x){spill})
#       names(comp) <- pData(gs)$name
#       gs@compensation <- comp
#       rval$gating_set <- gs
#       #load("../flowR_utils/demo-data/Rafa2Gui/analysis/cluster.rda")
#       #fs <- build_flowset_from_df(df = res$cluster$data, origin = res$cluster$flow_set)
#       #gs <- flowWorkspace::GatingSet(fs)
#       #gs@transformation <-  res$cluster$transformation
#       #add_gates_flowCore(gs, res$cluster$gates)
#       #rval$gating_set <- gs
#       #plot_params$sample <- pData(gs)$name[1]
#       #utils::data("GvHD", package = "flowCore")
#       #rval$gating_set <- GatingSet(GvHD)
#       #gs <- load_gs("./inst/ext/gs")
#       #rval$gating_set <- gs
#     })
# 
#     res <- callModule(Gating, "module", rval = rval)
# 
#   }
# 
#   shinyApp(ui, server)
# 
# }

