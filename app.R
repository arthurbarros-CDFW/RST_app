# app.R
# CDFW RST passage estimator
# for rotary screw trap data

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
  tags$head(
    tags$style(HTML("
    .shiny-notification {
      opacity: 1; /* 0 is fully transparent, 1 is fully opaque */
    }
  "))
  ),
  sidebarLayout(
    sidebarPanel(
      #file input for multiple files
      fileInput("files",label=NULL,
                buttonLabel="Step 1: Upload Data",
                accept = c(".csv"), multiple = TRUE),
      uiOutput("fileList"),  #output for the list of uploaded files

      verbatimTextOutput("dataListNames"),  #output to display the names in dataList
      
      actionButton("run_models", "Step 2: Run Model Comparison"),
      
      dateInput("survey_start","Survey Start Date:", value = "2022-01-19"),
      dateInput("survey_end","Survey End Date:", value = "2022-06-22"),
      selectInput("sum.by","Sum by:",c("day",
                                       "week",
                                       "month",
                                       "year"),selected = "week"),
      
      selectInput("target_species", "Target Species",
                  choices = c("Chinook salmon","NA"),selected="Chinook salmon"),
      selectInput("target_run", "Target Run (optional)",
                  choices = c("All runs" = "", "Fall", "Late Fall", "Spring", "Winter"),
                  selected = ""),
      numericInput("min_sample_size", "Minimum sample size for spline modeling",
                   value=10, min=5,max=20),
      
      checkboxInput("impute_all", "Impute All Efficiency Values", value = FALSE),
      checkboxInput("use_discharge","Use discharge as a covariate of efficiency?", value=FALSE),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Model reports", 
                 uiOutput("model_results"),
                 actionButton("run_selected_models", "Step 3: Estimate Passage with Selected Models
                              ")
                 ),
        tabPanel("Passage plot", plotOutput("p_passage",height = "400px"),
                 downloadButton("download_passage_plot", "Download Plot")),
        tabPanel("Passage table", DTOutput("passage_result"),
                 downloadButton("download_passage_table", "Download Estimates")),
        tabPanel("Catch plot", plotOutput("p_catch", height = "400px"),
                 downloadButton("download_catch_plot", "Download Plot")),
        tabPanel("Efficiency plot", plotOutput("p_eff", height = "400px"),
                 downloadButton("download_eff_plot", "Download Plot")),
        tabPanel("User Guide", 
                 tags$iframe(src = "./RST_app_documentation.html", 
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
  
  ###########################
  #reactive values to store uploaded files,
  #assigned datasets, and estimates
  ###########################
  uploadedFiles <- reactiveVal(list())
  datasets <- reactiveValues(catch = NULL, 
                             recapture = NULL, 
                             release = NULL,
                             visit = NULL)
  
  model_results<-reactiveVal(NULL)
  eff_data_unimputed <- reactiveVal(NULL)
  selected_models <- reactiveVal(list())
  eff_warnings <- reactiveVal(list())
  model_objects<-reactiveVal(NULL)
  model_types<-reactiveVal(NULL)
  
  catch_results<-reactiveVal(NULL)
  passage_result <- reactiveVal(NULL)
  
  plot_catch<-reactiveVal(NULL)
  plot_eff<-reactiveVal(NULL)
  plot_passage<-reactiveVal(NULL)

  status_message <- reactiveVal("Waiting for data upload...")
  
  #initialize dataListNames before upload
  output$dataListNames <- renderPrint({
    cat("Uploaded files:\n")
    cat("  (No files uploaded yet)\n")
    cat("\nAssigned datasets:\n")
    cat("  catch: ", ifelse(!is.null(datasets$catch), "✓ Loaded", "✗ Missing"), "\n")
    cat("  recapture: ", ifelse(!is.null(datasets$recapture), "✓ Loaded", "✗ Missing"), "\n")
    cat("  release: ", ifelse(!is.null(datasets$release), "✓ Loaded", "✗ Missing"), "\n")
    cat("  visit: ", ifelse(!is.null(datasets$visit), "✓ Loaded", "✗ Missing"), "\n")
  })
  
  ###########################
  #update list of uploaded files
  ###########################
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
  
  ###########################
  #trigger eff model comparison on button click
  ###########################
  observeEvent(input$run_models, {
    
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
    
    if(input$use_discharge == TRUE){
      #check if discharge exists
      discharge_cols <- grep("discharge", names(datasets$visit), value = TRUE, ignore.case = TRUE)
      
      if(length(discharge_cols) == 0){
        showNotification(
          "Error: 'Include discharge' option is selected but no discharge data found in the visit data. Please either:\n
        1. Upload a file with a 'discharge' column\n
        2. Uncheck the 'Use discharge as a covariate' option",
          type = "error", 
          duration = 10
        )
        status_message("ERROR: Discharge data required but not found")
        return()
      }
      
      #check if discharge is all NA or empty
      discharge_data <- datasets$visit[[discharge_cols[1]]]
      if(all(is.na(discharge_data)) || length(discharge_data) == 0){
        showNotification(
          "Error: 'Include discharge' option is selected but the discharge data is empty or all NA. Please either:\n
        1. Upload a file with valid discharge data\n
        2. Uncheck the 'Use discharge as a covariate' option",
          type = "error", 
          duration = 10
        )
        status_message("ERROR: Discharge data is empty or all NA")
        return()
      }
    }
    
    if(input$use_discharge==TRUE){
      #fix discharge function
      fix_discharge_field <- function(data) {
        col_idx <- which(tolower(names(data)) == "discharge")
        if (length(col_idx) > 0) {
          #rename the first match to "discharge"
          names(data)[col_idx[1]] <- "discharge"
        } else {
          warning("No 'discharge' field found in the dataset")
        }
        return(data)
      }
      
      datasets$visit <- fix_discharge_field(datasets$visit)
    }
    
    
    #set parameters
    survey_start <- input$survey_start
    survey_end <- input$survey_end
    target_species <- input$target_species
    target_run <- if(input$target_run == "") NA else input$target_run
    sum.by <- input$sum.by
    impute_all <- input$impute_all
    use_discharge<-input$use_discharge
    
    #display a progress bar while running the sourced scripts
    withProgress(message = "Developing models...", value = 0, {
      #increment progress for each script
      
      ###############################################
      #run catch model
      ###############################################
      incProgress(0.25, detail = "Modeling catch")
      catch.output <- tryCatch({
        est_catch(target_species = target_species,
                  target_run = target_run,
                  catch_data = datasets$catch,
                  visit_data = datasets$visit,
                  survey_start = survey_start,
                  survey_end = survey_end)
      }, error = function(e) {
        showNotification(paste("Error in catch estimation:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      catch_plot<-catch.output$p_catch
      
      if (!is.null(catch_plot)) {
        plot_catch(catch_plot)
      }
      
      if(is.null(catch.output)) return()
      
      catch.fits<-catch.output$models
      catch.results<-catch.output$results
      catch.X.miss<-catch.output$X.miss
      
      catch_results<-catch_results(catch.output)
      
      ###############################################
      #run eff model comparison
      ###############################################
      incProgress(0.5, detail = "Comparing efficiency models")
      eff.comparison <- tryCatch({
        compare_efficiency_models(release_data = datasets$release,
                       recapture_data = datasets$recapture,
                       visit_data = datasets$visit,
                       impute_all = impute_all,
                       survey_start = survey_start,
                       survey_end = survey_end,
                       min_sample_size = input$min_sample_size,
                       use_discharge = use_discharge)
      }, error = function(e) {
        showNotification(paste("Error in efficiency estimation:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      if(is.null(eff.comparison)) return()

      eff.model.comparisons<-eff.comparison$comparison.df
      eff.candidate_models <- eff.comparison$candidate_models
      eff.data.unimputed<-eff.comparison$eff
      eff.model.types<-eff.comparison$eff.type
      eff_data_unimputed(eff.data.unimputed)
      
      #get unique id for traps
      traps <- unique(catch.results$trap_ID_decimal)
      
      all_model_results <- list()
      all_model_objects <- list()
      all_model_types <- list()
      
      #build results for each trap
      for(trap in traps){
        trap_label=as.character(trap)
        if(!is.null(eff.model.comparisons[[trap_label]])){
          comp_df<-eff.model.comparisons[[trap_label]]
          
          #add trap info to df
          if(is.data.frame(comp_df) && nrow(comp_df)>0){
            comp_df$Trap<-trap_label
            comp_df<-comp_df[,c("Trap",setdiff(names(comp_df),"Trap"))]
            all_model_results[[trap_label]]<-comp_df
            
            candidate_models <- eff.candidate_models[[trap_label]]
            all_model_objects[[trap_label]] <- candidate_models
            
            candidate_types<-eff.model.types[[trap_label]]
            all_model_types[[trap_label]]<-candidate_types
          }
        } else {
          #create placeholder for traps without eff models
          placeholder<-data.frame(
            Trap=trap_label,
            Model="No efficiency model",
            AIC=NA,
            AICc=NA,
            AICc_diff=NA
          )
          all_model_results[[trap_label]]<-placeholder
          all_model_objects[[trap_label]] <- NULL
          all_model_types[[trap_label]] <- NULL
        }
      }
      
      #store model results
      if (length(all_model_results)>0) {
        model_results(all_model_results)
        model_objects(all_model_objects)
        model_types(all_model_types)
        default_selected<-list()
        for(trap_name in names(all_model_results)){
          model_data<-all_model_results[[trap_name]]
          if(!is.null(model_data) && nrow(model_data)>0){
            #select first model as default
            default_selected[[trap_name]]<-model_data[1, ,drop=FALSE]
          }
        }
        
        selected_models(default_selected) #
      }
      
    })
  })
  
  ###########################  
  #display model comparisons
  ###########################
  output$model_results <- renderUI({
    req(model_results())
    
    model_list<-model_results()
    
    #if no models message
    if(is.null(model_list) || length(model_list) == 0) {
      return(h4("No model results available"))
    }
    
    #create list of ui elements
    output_list<-list()
    
    #add header
    output_list[[length(output_list) + 1]] <- h3("Efficiency Model Results")
    output_list[[length(output_list) + 1]] <- tags$hr()
    
    #for each trap display its table
    for(trap_name in names(model_list)){
      output_list[[length(output_list) + 1]] <- h4(paste("Trap:", trap_name))
      
      #create uid for traps table
      table_id <- paste0("model_table_", gsub("[^A-Za-z0-9]", "_", trap_name))
      
      #use renderUI with dtoutput
      output_list[[length(output_list) + 1]] <- DTOutput(table_id)
      
      #add some spacing
      output_list[[length(output_list) + 1]] <- tags$br()
      
      #render the table for this trap
      local({
        trap<-trap_name
        table_id_local<-table_id
        model_data<-model_list[[trap]]
        
        #store model data in a separate reactive value to avoid re-rendering
        #which caused screen flickering/flashing
        output[[table_id_local]]<-renderDT({
          if(is.null(model_data)||nrow(model_data)==0){
            #if no data show message
            datatable(
              data.frame(Message="No model comparison data for this trap"),
              options=list(dom='t'),
              rownames = FALSE
            )
          } else {
            #format columns
            numeric_cols<-which(sapply(model_data,is.numeric))
            
            #determine which row is currently selected
            selected_trap <- isolate(selected_models()[[trap]])
            selected_row <- 1  #default to first row
            
            #try to find current selection
            if(is.null(selected_trap) || nrow(selected_trap)==0){
              selected_row<-1 #default to 1st row
              } else {
              for(i in 1:nrow(model_data)){
                if(identical(model_data[i,],selected_trap[1,])){
                  selected_row<-i
                    break
                }
              }
                if(is.null(selected_row)){
                  selected_row<-1 #if selected model isn't found default to 1
                }
            }
            
            datatable(
              model_data,
              options=list(
                pageLength=10,
                dom="Bfrtip",
                scrollX=TRUE,
                columnDefs=list(
                  list(className="dt-center",targets='_all')
                )
              ),
              rownames = FALSE,
              selection = list(mode = 'single', selected = selected_row),
              class = 'display compact stripe hover'
            ) %>%
              formatRound(columns = numeric_cols, digits = 3)
          }
        })
      })
    }
    return(output_list)
  })
  ###########################
  #observe row selections from model tables
  ###########################
  observe({
    #get all model results
    model_list<-model_results()
    if (is.null(model_list) || length(model_list)==0) return()
      
    #only update if there's a selection change
    new_selected<-list()
    selection_changed <- FALSE
    
    #check trap tables for selections
    for(trap_name in names(model_list)){
      table_id<-paste0("model_table_",gsub("[^A-Za-z0-9]", "_", trap_name))
      selected_rows<-input[[paste0(table_id,"_rows_selected")]]
      
      if(!is.null(selected_rows) && length(selected_rows) > 0) {
        #get selected model data
        model_data<-model_list[[trap_name]]
        if(!is.null(model_data) && nrow(model_data) >= selected_rows[1]) {
          new_selected[[trap_name]] <- model_data[selected_rows[1], , drop=FALSE]
        }
      } else {
        #if no selection, keep default
        model_data<-model_list[[trap_name]]
        if(!is.null(model_data) && nrow(model_data) > 0){
          new_selected[[trap_name]] <- model_data[1, , drop=FALSE]
        }
      }
    }
    #only update if the selection actually changed
    current_selected <- selected_models()
    if(length(new_selected) != length(current_selected) || 
       !identical(new_selected, current_selected)) {
      selected_models(new_selected)
    }
  })
  
  ###########################
  #run eff imputation and passage estimation after model selection
  ###########################
  observeEvent(input$run_selected_models,{
    
    #get selected models
    selected<-selected_models()
    if (is.null(selected) || length(selected) == 0) {
      cat("No models selected.\n")
      cat("Click on a row in the model comparison tables above to select a model.")
      return()
    }
    
    #get model objects
    model_objects <- model_objects() 
    
    #get model types
    model_types <- model_types()
    
    selected_model_objects <- list()
    selected_model_types <-list()
    model_names=list()
    
    for(trap_name in names(selected)){
      #get selected model name and candidate models for trap
      selected_model_name<-selected[[trap_name]]$Model[1]
      candidate_models<-model_objects[[trap_name]]
      candidate_types<-model_types[[trap_name]]
      
      #find candidate model that matches selected
      if(!is.null(candidate_models) && selected_model_name %in% names(candidate_models)){
        selected_model_objects[[trap_name]]<-candidate_models[[selected_model_name]]
        selected_model_types[[trap_name]]<-candidate_types[[selected_model_name]]
        model_names[[trap_name]]<-selected_model_name
      } else {
        #if model not found get first available
        if(!is.null(candidate_models) && length(candidate_models)>0){
          first_name<-names(candidate_models)[1]
          selected_model_objects[[trap_name]]<-candidate_models[[first_name]]
          model_names[[trap_name]]<-first_name
          showNotification(paste("Using fallback model for",trap_name,":", first_name),type="warning")
        }else{
          showNotification(paste("No model object found for", trap_name), type = "error")
          return()
        }
      }
    }
    
    
    #debug output to verify
    print("Selected models:")
    print(model_names)
    
    #display a progress bar while running the sourced scripts
    withProgress(message = "Running selected efficiency models...", value = 0, {
      #increment progress for each script
      
      ################################################
      #Impute efficiency
      ################################################
      incProgress(0.25, detail = "Imputing efficiency")
      model_eff_output<-tryCatch({
        impute_efficiency(
          efficiency_data = eff_data_unimputed(),
          selected_models = selected_model_objects,
          model_names = model_names,
          impute_all = input$impute_all
        )
      }, error = function(e) {
        showNotification(paste("Error in efficiency imputation:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      #display discharge warnings if exist
      if (!is.null(model_eff_output) && !is.null(model_eff_output$warnings) && 
          length(model_eff_output$warnings) > 0) {
        #show as notification
        for (warn in model_eff_output$warnings) {
          if(length(warn)>0){
            showNotification(
              HTML(paste0("<b>Warning:</b><br>", gsub("\n", "<br>", warn))),
              type = "warning",
              duration = 100
            )
          }
        }
      }
      
      ################################################
      #plot efficiency
      ################################################
      incProgress(0.5, detail = "Plotting efficiency")
      eff_plot<-tryCatch({
        plot_efficiency(
          impute_eff_results=model_eff_output$results,
          use_discharge=input$use_discharge,
          impute_all = input$impute_all
        )
      }, error = function(e) {
        showNotification(paste("Error in efficiency plotting:", e$message), type = "error", duration = 10)
        return(NULL)
      })
      
      if (!is.null(eff_plot)) {
        plot_eff(eff_plot)
      }
      
      ################################################
      #Estimate and bootstrap passage
      ################################################
      incProgress(0.75, detail = "Estimating and bootstrapping passage")
      tryCatch({
        est_passage_output<-est_passage(
          catch.results=catch_results(),
          eff.results=model_eff_output,
          summarize.by=input$sum.by,
          survey_start=input$survey_start,
          survey_end=input$survey_end,
          catch.fits=catch_results()$models,
          eff.fits=model_eff_output$models,
          eff.types=selected_model_types,
          target_species=input$target_species,
          target_run=input$target_run,
          use_discharge=input$use_discharge,
          min_sample_size = input$min_sample_size
          )
        }, error = function(e) {
          showNotification(paste("Error in passage estimation:", e$message), type = "error", duration = 10)
          return(NULL)
      })
      
      passage_plot<-est_passage_output$p_passage
      if (!is.null(passage_plot)) {
        plot_passage(passage_plot)
      }
      
      passage_result(est_passage_output$passage_output)
      
    })

  })
  
  ###########################
  #output plots and tables
  ###########################
  
  #render p_catch
  output$p_catch <- renderPlot({
    req(plot_catch())
    print(plot_catch())
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
  
  #render passage table
  #display passage results
  output$passage_result <- renderDataTable({
    req(passage_result())
    passage_result()
  })
  
  ###########################
  #download handlers
  ###########################
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