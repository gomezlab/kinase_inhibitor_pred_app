library(shiny)
library(tidyverse)
library(dqshiny)
library(here)
library(tidymodels)
library(reactable)
library(shinyjs)
library(digest)
library(rhdf5)

options(shiny.maxRequestSize=30*1024^2)

dir.create(here('data/uploads'), showWarnings = F)

all_geo_archs_ids = read_rds(here('data/ARCHS_GEO_IDs.rds'))

source(here('functions.R'))

# Define UI for application that draws a histogram
ui <- fluidPage(
	
	shinyjs::useShinyjs(), 
	
	# Application title
	titlePanel("Kinase Inhibitor Human Cell Viability Prediction"),
	
	tags$p("Hello and welcome to the kinase inhibitor cell viability prediction web server. This website is a companion 
				 to a forthcoming publication concerning prediction of cell viability after exposure to kinase inhibitors. 
				 The primary model developed in this paper uses RNAseq gene expression and kinase inhibitor target profiles
				 to make these predictions. We've made this server available to allow interested biologists to submit gene expression
				 data they have gathered where there is some interest in how a set of kinase inhbitors would affect their model 
				 system."),
	
	tags$p("This model is built to work with human cell line RNAseq data. Unfortunately, the model will not work with 
				 data from any other organism."),
	
	tags$hr(),
	
	# Sidebar with a slider input for number of bins 
	sidebarLayout(
		sidebarPanel(
			
			tags$h2("Option 1: Upload RNAseq Results:"),
			fileInput("RNAseq_file", "Please select your quant.sf file from salmon",
								multiple = TRUE),
			tags$hr(),
			
			tags$h2("Option 2: Specify a GEO ID:"),
			autocomplete_input("GEO_ARCHS_ID", 
												 "All IDs start with GSM",
												 all_geo_archs_ids,
												 placeholder = "Start Typing to Find Your GEO ID",
												 max_options = 10),
			
			actionButton("submit_geo", "Submit GEO ID"),
			actionButton("submit_random_geo", "Submit Random GEO ID")
		),
		
		mainPanel(
			tags$div(id = "instructions",
							 tags$h2("Application Instructions:"),
							 
							 tags$p("There are currently two ways to input your RNAseq data to make kinase inhibitor cell 
							 				 viability predictions. The first is to upload the \"quant.sf\" file from the salmon (or compatible) 
							 				 read aligner. Alternatively, you can input a GEO database ID which has been preprocessed by the ARCHS 
							 			 project."),
							 
							 tags$p("After inputing a data set, the processing pipeline will organize your data and search for the data 
							 			 related to the genes used in the model. Then the model will be loaded and cell viability predictions 
							 			 made for your data set for all 229 compounds in the Klaeger et al set. Finally, a preview of the 
							 			 results will be displayed with an option to download the full predictions and a summary document 
							 			 highlighting some of the most interesting compounds."),
							 
							 tags$p("The processing should take less than a minute and progress indicators will appear in the bottom 
							 			 corner."),
			),
			
			
			tags$div(id = "results", 
							 tags$h1("RNAseq Data Checks")			 
			),
			
			textOutput("RNAseq_qc_text"),
			fluidRow(
				column(6,tableOutput("RNAseq_sample")),
				
				column(6,tableOutput("prediction_sample"))
			),
			
			downloadButton("model_predictions_download", label = "Download Model Predictions"),
		)
	)
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
	global_data <- reactiveValues(RNAseq = NULL)
	
	observeEvent(input$RNAseq_file, {
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(1/3, detail = "Processing RNAseq Data")
		
		hide("instructions")
		
		TPM_data = read_delim(input$RNAseq_file$datapath, delim = "\t") %>%
			convert_salmon_to_HGNC_TPM()
		
		file.copy(input$RNAseq_file$datapath, here('data/uploads',paste0(substr(digest(TPM_data), 1, 6))))
		global_data$RNAseq = TPM_data
	})
	
	observeEvent(input$submit_geo, {
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(1/3, detail = "Processing GEO Data")
		
		hide("instructions")
		
		archs_data = H5Fopen(here('data/ARCHS_subset/matt_model_matrix.h5'))
		
		GEO_col = which(archs_data$meta$samples$geo_accession == input$GEO_ARCHS_ID)
		
		global_data$RNAseq = data.frame(hgnc_symbol = archs_data$meta$genes$genes, TPM = archs_data$data$expression[GEO_col,])
	})
	
	observeEvent(input$submit_random_geo, {
		random_geo_id = sample(all_geo_archs_ids,1)
		
		update_autocomplete_input(session, "GEO_ARCHS_ID", value = random_geo_id)
		
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(1/3, detail = "Processing GEO Data")
		
		hide("instructions")
		
		archs_data = H5Fopen(here('data/ARCHS_subset/matt_model_matrix.h5'))
		
		GEO_col = which(archs_data$meta$samples$geo_accession == random_geo_id)
		
		global_data$RNAseq = data.frame(hgnc_symbol = archs_data$meta$genes$genes, TPM = archs_data$data$expression[GEO_col,])
	})
	
	output$RNAseq_qc_text <- renderText({
		if (is.null(global_data$RNAseq)) return()
		
		paste0("Your file contains ", dim(global_data$RNAseq)[1], '/110 genes.')
	})
	
	model_predictions <- reactive({
		if (is.null(global_data$RNAseq)) return()
		
		prediction_results = make_predictions(global_data$RNAseq)
		
		shinyjs::show("results")
		shinyjs::show("model_predictions_download")
		
		prediction_results
	})
	
	output$RNAseq_sample <- renderTable({
		if (is.null(global_data$RNAseq)) return()
		
		return(head(global_data$RNAseq))
	})
	
	output$prediction_sample <- renderTable({
		if (is.null(global_data$RNAseq)) return()
		
		return(head(model_predictions()))
	})
	
	output$model_predictions_download <- downloadHandler(
		filename = function() {
			paste0("kinase_inhbitor_model_predictions_",paste0(substr(digest(global_data$RNAseq), 1, 6)),".csv")
		}, 
		content = function(file) {
			write_csv(model_predictions(), file)
		})
	
	shinyjs::hide("model_predictions_download")
	shinyjs::hide("results")
}

# Run the application 
shinyApp(ui = ui, server = server)