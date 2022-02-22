library(shiny)
library(tidyverse)
library(dqshiny)
library(here)
library(tidymodels)
library(reactable)
library(shinyjs)
library(digest)
library(rmarkdown)
library(rhdf5)

options(shiny.maxRequestSize=30*1024^2)

dir.create(here('data/uploads'), showWarnings = F)
dir.create(here('www'), showWarnings = F)

all_geo_archs_ids = read_rds(here('data/ARCHS_GEO_IDs.rds'))

source(here('functions.R'))

# Define UI for application that draws a histogram
ui <- fluidPage(
	
	shinyjs::useShinyjs(), 
	
	titlePanel("Kinase Inhibitor Human Cell Viability Prediction"),
	
	tags$p("Hello and welcome to the kinase inhibitor cell viability prediction web server. This website is a companion 
				 to a forthcoming publication concerning prediction of cell viability after exposure to kinase inhibitors. 
				 The primary model developed in this paper uses RNAseq gene expression and kinase inhibitor target profiles
				 to make these predictions. We've made this server available to allow interested biologists to submit gene expression
				 data they have gathered where there is some interest in how a set of kinase inhbitors would affect their model 
				 system."),
	
	tags$p("This model is built to work with human cell line RNAseq data and makes cell viability predictions for 229 kinase inhibitors."),
	
	tags$hr(),
	
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
							 				 read aligner produced using the ensemble transcript set (ENST). Alternatively, you can input a GEO 
							 			 database ID which has been preprocessed by the ARCHS project."),
							 
							 tags$p("After inputing a data set, the processing pipeline will organize your data and search for the data 
							 			 related to the genes used in the model. Then the model will be loaded and cell viability predictions 
							 			 made for your data set for all 229 compounds in the Klaeger et al set. Finally, a preview of the 
							 			 results will be displayed with an option to download the full predictions and a summary document 
							 			 highlighting some of the most interesting compounds."),
							 
							 tags$p("The processing should take less than a minute and progress indicators will appear in the bottom 
							 			 corner."),
			),
			
			fluidRow(
				column(12,textOutput("RNAseq_qc_text"))
			),
			
			fluidRow(
				downloadButton("model_predictions_download", label = "Download Model Predictions"),
				downloadButton("predictions_summary_docx_download", label = "Download Predictions Summary - DOCX Format")
				# downloadButton("predictions_summary_pdf_download", label = "Download Predictions Summary - PDF Format")
			)
			
		)
	),
	
	hr(),
	
	fluidRow(id = "footnotes",
		column(11,
					 p("The model used in this system is based on data from ", 
					 	a(href="https://www.theprismlab.org/","PRISM", .noWS = "outside"), ", ", 
					 	a(href="https://sites.broadinstitute.org/ccle/","CCLE", .noWS = "outside")," and ", 
					 	a(href="http://dx.doi.org/10.1126/science.aan4368","Klaeger et. al", .noWS = "outside"), 
					 	". Submission through GEO ID facilitated through preprocesing of RNAseq data by ", 
					 	a(href="https://maayanlab.cloud/archs4/","ARCHS", .noWS = "outside"), ".",
					 	.noWS = c("after-begin", "before-end"))
		),
		
		column(1,
					 a(href="https://github.com/mbergins/kinase_inhibitor_pred_app", icon("github", class="fa-2x", style="float:right;")))
	)
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
	global_data <- reactiveValues(RNAseq = NULL,
																model_predictions = NULL,
																model_id = NULL)
	
	observeEvent(input$RNAseq_file, {
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(1/3, detail = "Processing RNAseq Data")
		
		hide("instructions")
		
		TPM_data = read_delim(input$RNAseq_file$datapath, delim = "\t") %>%
			convert_salmon_to_HGNC_TPM()
		
		global_data$model_id = substr(digest(TPM_data), 1, 6)
		
		file.copy(input$RNAseq_file$datapath, here('data/uploads',global_data$model_id))
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

	observeEvent(global_data$RNAseq, {
		global_data$model_id = substr(digest(global_data$RNAseq), 1, 6)
	})
		
	observeEvent(global_data$model_id, {
		run_model()
	})
	
	output$RNAseq_qc_text <- renderText({
		if (is.null(global_data$model_predictions)) return()
		
		paste0("Your file contains ", dim(global_data$RNAseq)[1], '/110 genes and ', dim(global_data$model_predictions)[1], " were made.")
	})
	
	run_model <- reactive({
		if (is.null(global_data$model_id)) return()
		
		global_data$model_predictions = make_predictions(global_data$RNAseq)
		
		render('build_inhibitor_overview.Rmd', 
					 output_file = here('www/',paste0("kinase_inhibitor_summary_",global_data$model_id,".docx")), 
					 params = list(predictions = global_data$model_predictions, RNAseq_data = global_data$RNAseq, model_id = global_data$model_id))
		
		# render('build_inhibitor_overview.Rmd', 
		# 			 output_file = here('www/',paste0("kinase_inhibitor_summary_",global_data$model_id,".pdf")), 
		# 			 params = list(predictions = prediction_results, RNAseq_data = global_data$RNAseq, model_id = global_data$model_id))

		shinyjs::show("model_predictions_download")
		shinyjs::show("predictions_summary_docx_download")
		# shinyjs::show("predictions_summary_pdf_download")
	})
	
	output$RNAseq_sample <- renderTable({
		if (is.null(global_data$model_id)) return()
		
		return(head(global_data$RNAseq))
	})
	
	output$prediction_sample <- renderTable({
		if (is.null(global_data$model_predictions)) return()
		
		return(head(global_data$model_predictions))
	})
	
	output$model_predictions_download <- downloadHandler(
		filename = function() {
			paste0("kinase_inhbitor_model_predictions_",global_data$model_id,".csv")
		}, 
		content = function(file) {
			write_csv(global_data$model_predictions %>% pivot_wider(names_from = concentration_M, values_from = predicted_viability), file)
		})
	
	output$predictions_summary_docx_download <- downloadHandler(
		filename = function() {
			paste0("kinase_inhibitor_summary_",paste0(global_data$model_id),".docx")
		}, 
		content = function(file) {
			file.copy(here('www/',paste0("kinase_inhibitor_summary_",global_data$model_id,".docx")), file)
		})
	
	# output$predictions_summary_pdf_download <- downloadHandler(
	# 	filename = function() {
	# 		paste0("kinase_inhibitor_summary_",paste0(global_data$model_id),"pdf")
	# 	}, 
	# 	content = function(file) {
	# 		file.copy(here('www/',paste0("kinase_inhibitor_summary_",global_data$model_id,"pdf")), file)
	# 	})
	
	shinyjs::hide("model_predictions_download")
	shinyjs::hide("predictions_summary_docx_download")
	# shinyjs::hide("predictions_summary_pdf_download")
	shinyjs::hide("results")
}

# Run the application 
shinyApp(ui = ui, server = server)