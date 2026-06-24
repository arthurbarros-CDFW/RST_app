# app.R
# CDFW RST passage estimator
# for rotary screw trap data analysis

library(shiny)
library(shinyWidgets)
library(shinythemes)
library(tidyverse)
library(DT)
library(ggplot2)
library(mvtnorm)
library(splines)

#source function files
sapply(list.files("scripts/functions", pattern = "\\.R$", full.names = TRUE), source)

#global variables
if (!exists("gap_threshold_days")) gap_threshold_days <- 7
if (!exists("time_zone")) time_zone <- "America/Los_Angeles"
if (!exists("unassd.sig.digit")) unassd.sig.digit <- 1
if (!exists("knotMesh")) knotMesh <- 15
if (!exists("max.ok.gap")) max.ok.gap <- 2
if (!exists("bootstrap.CI.fx")) bootstrap.CI.fx <- "f.ci"

ui <- fluidPage(
  theme = shinytheme("cosmo"),
  titlePanel("CDFW RST passage estimator"),
  sidebarLayout(
    sidebarPanel(
      #file input for multiple files
      fileInput("files", "Upload Data Files", accept = c(".csv"), multiple = TRUE),
      h5(strong("Data that must be uploaded:")),
      h5("catch, recapture, release, visit"),
      
      actionButton("run_estimate", "Estimate Passage"),
      
      dateInput("survey_start","Survey Start Date:", value = "2022-01-19"),
      dateInput("survey_end","Survey End Date:", value = "2022-06-22"),
      selectInput("sum.by","Sum by:",c("day",
                                       "week",
                                       "month",
                                       "year"),selected = "week"),
      
      selectInput("target_species", "Target Species",
                  choices = c("Chinook salmon","NA"),selected="Chinook salmon"),
      selectInput("target_run", "Target Run (optional)",
                  choices = c("All runs" = "", "Fall", "Spring", "Winter", "Summer"),
                  selected = ""),
      
      checkboxInput("impute_all", "Impute All Efficiency Values", value = FALSE),
      
      uiOutput("fileList"),  # Output for the list of uploaded files
      h4("Data List:"),
      verbatimTextOutput("dataListNames"),  #output to display the names in dataList
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Passage plot", plotOutput("p_passage",height = "400px"),
                 downloadButton("download_passage_plot", "Download Plot")),
        tabPanel("Passage table", DTOutput("passage_result"),
                 downloadButton("download_passage_table", "Download Estimates")),
        tabPanel("Catch plot", plotOutput("p_catch", height = "400px"),
                 downloadButton("download_catch_plot", "Download Plot")),
        tabPanel("Efficiency plot", plotOutput("p_eff", height = "400px"),
                 downloadButton("download_eff_plot", "Download Plot")),
        tabPanel("User Guide", 
                 tags$iframe(src = "RST_app_documentation.html", 
                             width = "100%", 
                             height = "800px", 
                             frameborder = "0",
                             style = "border: none;"))
      )
    )
  )
)
server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100*1024^2) # sets max file size 10 100MB
  
  #reactive values to store uploaded files, assigned datasets, and estimates
  uploadedFiles <- reactiveVal(list())
  datasets <- reactiveValues(catch = NULL, 
                             recapture = NULL, 
                             release = NULL,
                             visit = NULL)
  
  passage_result <- reactiveVal(NULL)
  plot_catch<-reactiveVal(NULL)
  plot_eff<-reactiveVal(NULL)
  plot_passage<-reactiveVal(NULL)

  status_message <- reactiveVal("Waiting for data upload...")
  
  #update list of uploaded files
  observe({
    files <- input$files
    if (is.null(files)) {
      status_message("Waiting for data upload...")
      return()
    }
    
    status_message("Reading uploaded files...")
    
    uploadedFiles(files)
    
    #read uploaded files into a list
    dataList <- list()
    for (i in 1:nrow(files)) {
      file_path <- files$datapath[i]
      file_name <- tools::file_path_sans_ext(files$name[i])
      dataList[[file_name]] <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    }
    
    #display names in data list
    output$dataListNames <- renderPrint({
      cat("Uploaded files:\n")
      for (name in names(dataList)) {
        cat("  - ", name, "\n")
      }
      cat("\nAssigned datasets:\n")
      cat("  catch: ", ifelse(!is.null(datasets$catch), "✓ Loaded", "✗ Missing"), "\n")
      cat("  recapture: ", ifelse(!is.null(datasets$recapture), "✓ Loaded", "✗ Missing"), "\n")
      cat("  release: ", ifelse(!is.null(datasets$release), "✓ Loaded", "✗ Missing"), "\n")
      cat("  visit: ", ifelse(!is.null(datasets$visit), "✓ Loaded", "✗ Missing"), "\n")
    })
    
    #safely assign data to reactive variables based on file names
    file_names_lower <- tolower(names(dataList))
    
    datasets$catch <- if(any(grepl("catch", file_names_lower))) 
      dataList[[which(grepl("catch", file_names_lower))[1]]] else NULL
    
    datasets$recapture <- if(any(grepl("recapture", file_names_lower))) 
      dataList[[which(grepl("recapture", file_names_lower))[1]]] else NULL
    
    datasets$release <- if(any(grepl("release", file_names_lower))) 
      dataList[[which(grepl("release", file_names_lower))[1]]] else NULL
    
    datasets$visit <- if(any(grepl("visit", file_names_lower)) | any(grepl("trap", file_names_lower))) {
      idx <- which(grepl("visit", file_names_lower) | grepl("trap", file_names_lower))
      dataList[[idx[1]]]
    } else NULL
    
    if (!is.null(datasets$catch) && !is.null(datasets$recapture) && 
        !is.null(datasets$release) && !is.null(datasets$visit)) {
      status_message("All datasets loaded. Ready to estimate passage.")
    } else {
      status_message("Missing some datasets. Please upload all four files.")
    }
    
    })
  
  #trigger passage estimation on button click
  observeEvent(input$run_estimate, {
    
    #ensure all datasets are available
    missing_datasets <- c()
    if (is.null(datasets$catch)) missing_datasets <- c(missing_datasets, "catch")
    if (is.null(datasets$recapture)) missing_datasets <- c(missing_datasets, "recapture")
    if (is.null(datasets$release)) missing_datasets <- c(missing_datasets, "release")
    if (is.null(datasets$visit)) missing_datasets <- c(missing_datasets, "visit")
    
    if (length(missing_datasets) > 0) {
      status_message(paste("ERROR: Missing datasets:", paste(missing_datasets, collapse = ", ")))
      showNotification(paste("Missing datasets:", paste(missing_datasets, collapse = ", ")), 
                       type = "error", duration = 5)
      return()
    }
    
    #set parameters
    survey_start <- input$survey_start
    survey_end <- input$survey_end
    target_species <- input$target_species
    target_run <- if(input$target_run == "") NA else input$target_run
    sum.by <- input$sum.by
    impute_all <- input$impute_all
    bootstrap <- input$bootstrap
    
    target_species<-"Chinook salmon"
    #target_run<-"Fall"
    
    sum.by=input$sum.by
    
    #display a progress bar while running the sourced scripts
    withProgress(message = "Running passage estimation...", value = 0, {
      #increment progress for each script
      
      incProgress(0.5, detail = "Running est_passage.R")
      results<-est_passage(catch=datasets$catch,
                  visits=datasets$visit,
                  release=datasets$release,
                  recapture=datasets$recapture,
                  summarize.by=sum.by, 
                  impute_all=impute_all,
                  bootstrap=T,
                  survey_start,survey_end,
                  target_species,target_run=target_run,
                  file.name="test")
      
    })
    rounded_results <-results$passage_output
    numeric_cols <- names(rounded_results)[sapply(rounded_results, is.numeric)]
    rounded_results[numeric_cols] <- lapply(rounded_results[numeric_cols], round, digits = 2)
    
    passage_result(rounded_results)
    
    
    #store plots
    if (!is.null(results$p_catch)) {
      plot_catch(results$p_catch)
    }
    if (!is.null(results$p_eff)) {
      plot_eff(results$p_eff)
    }
    if (!is.null(results$p_passage)) {
      plot_passage(results$p_passage)
    }
  })
    
    #display passage results
    output$passage_result <- renderDataTable({
      req(passage_result())
      passage_result()
    })
    
    #renderp_catch
    output$p_catch <- renderPlot({
      req(plot_catch())
      print(plot_catch())  # Make sure to print the ggplot object
    })
    
    #render p_eff
    output$p_eff <- renderPlot({
      req(plot_eff())
      print(plot_eff())
    })
    
    #render p_passage
    output$p_passage <- renderPlot({
      req(plot_passage())
      print(plot_passage())
    })
    
  
  #download handlers
  output$download_passage_table <- downloadHandler(
    filename = function() {
      paste("passage_estimates_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      #save the results to a temporary file
      write.csv(passage_result(), file,row.names = FALSE)
    }
  )
  
  output$download_passage_plot <- downloadHandler(
    filename = function() {
      paste("passage_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_passage(), device = "png", width = 10, height = 6, dpi = 300)
    }
  )
  
  output$download_catch_plot <- downloadHandler(
    filename = function() {
      paste("catch_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_catch(), device = "png", width = 10, height = 6, dpi = 300)
    }
  )
  
  output$download_eff_plot <- downloadHandler(
    filename = function() {
      paste("efficiency_plot_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = plot_eff(), device = "png", width = 10, height = 6, dpi = 300)
    }
  )
  
}

shinyApp(ui, server)