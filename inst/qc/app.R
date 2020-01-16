library(shiny)
library(flowWorkspace)
library(ggcyto)

wd <- "/fh/fast/gottardo_r/jkim2345_working"
studies <- grep("^SDY\\d+$", list.files(wd), value = TRUE)

ui <- fluidPage(
  titlePanel("QC gating set"),
  sidebarLayout(
    sidebarPanel(
      width = 2,
      selectInput("study", "Select a study:", studies),
      actionButton("load_study", "Load Study"),
      selectInput("gs", "Select a gating set:", NULL),
      actionButton("load_gs", "Load Gating Set"),
      selectInput("node", "Select a node:", NULL),
      actionButton("visualize", "Visualize")
    ),
    mainPanel(
      width = 10,
      tabsetPanel(
        # tabPanel("Study Summary"),
        tabPanel("Gating Set Summary",
          verbatimTextOutput("gs_summary"),
          plotOutput("gs_plot")
        ),
        tabPanel("Plot", uiOutput("gatePlot", height = "800px")),
        tabPanel("Pop",
          plotOutput("distPlot"),
          plotOutput("distPlot_pct"),
          dataTableOutput("popTab")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  gs <- reactiveVal()
  p <- reactiveVal()
  pop <- reactiveVal()

  observeEvent(input$load_study, {
    choices <- list.files(file.path(jkim, input$study, "GatingSet"))
    updateSelectInput(session, "gs", choices = choices)
  })

  observeEvent(input$load_gs, {
    if (input$gs != "") {
    path <- file.path(jkim, input$study, "GatingSet", input$gs)

      if (file.exists(path)) {
        gs(load_gs(path))
        choices <- gs_get_pop_paths(gs(), path = 1)[-1]
        updateSelectInput(session, "node", choices = choices)
      }
    }
  })

  output$gs_summary <- renderPrint({
    if (!is.null(gs())) {
      print(length(gs()))
      print(markernames(gs()))
    }
  })

  output$gs_plot <- renderPlot({
    if (!is.null(gs())) {
      plot(gs())
    }
  })

  output$gatePlot <- renderUI({
    p()
  })

  observeEvent(input$visualize, {
    GS <- isolate(gs())
    node <- isolate(input$node)
    by <- 20
    n <- length(GS)
    len <- seq_len(ceiling(n / by))
    plot_output_list <- lapply(len, function(i) {
      from <- 1 + (i - 1) * by
      to <- ifelse(i * by > n, n, i * by)
      plotname <- paste("plot", i, sep = "")
      output[[plotname]] <- renderPlot({
        autoplot(GS[seq(from, to)], node)
      })
      plotOutput(plotname)
    })

    viz <- do.call(tagList, plot_output_list)

    count <- gs_pop_get_stats(GS, node)
    percent <- gs_pop_get_stats(GS, node, type = "percent")
    stat <- merge(count, percent[, pop := NULL], by = "sample")
    pop(stat)
    p(viz)
  })

  output$distPlot <- renderPlot({
    if (!is.null(pop())) {
      hist(pop()$count)
    }
  })

  output$distPlot_pct <- renderPlot({
    if (!is.null(pop())) {
      hist(pop()$percent)
    }
  })

  output$popTab <- renderDataTable({
    if (!is.null(pop())) {
      pop()
    }
  })
}

shinyApp(ui = ui, server = server)

