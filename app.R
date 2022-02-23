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
CCLE_preds = read_rds(here('data/CCLE_prediction_summary.rds'))

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
			
			div(id = 'results',
					fluidRow(
						column(12,p("The model has finished running and a summary of your results follows. You will find two buttons
												at the bottom of this section to download a CSV file with the model predictions and a report with
												more background on understanding your results."))
					),
					fluidRow(
						column(12,textOutput("RNAseq_qc_text"))
					),
					
					hr(),
					
					fluidRow(
						column(12,
									 h2("Compound Viability Predictions"),
									 
									 p("To provide context for the predictions from your data, all of the following plots also show a summary
										of the predictions from the cell lines in the CCLE. The gray shaded region shows the range (95% coverage)
										of predictions for that compound, while the black line shows the average prediction. Predictions for your
										data appear as a blue line."),
									 
									 p("To help pick out potentially interesting compounds from the model predictions, we've sorted the 
										predictions using four different methods:"),
									 tags$ul(
									 	tags$li("The compound is predicted to have minimal effects on cell viability."),
									 	tags$li("The compound is predicted to have high average effect on cell viability."),
									 	tags$li("The compound is predicted to show a wide range of effect on cell viability."),
									 	tags$li("The compound is predicted to vary from the CCLE average effect.")
									 	
									 ),
									 
									 p("These categories are not mutually exclusive, so it's possible that a single compound will be present in 
										multiple compound sets. Otherwise, the results are displayed as a set of small multiple graphs with the 
										compound name in the title section."),
						)
					),
					
					fluidRow(
						column(6,
									 h3("Minimal Predicted Effect"),plotOutput("minimal_eff_preds")),
						column(6,
									 h3("High Predicted Effect"),plotOutput("high_eff_preds"))
					),
					
					fluidRow(
						column(6,
									 h3("High Range of Predicted Effect"),plotOutput("high_range_preds")),
						column(6,
									 h3("Large Difference with CCLE Lines"),plotOutput("ccle_diff_preds"))
					),
					
					hr(),
					
					fluidRow(
						downloadButton("model_predictions_download", label = "Download Model Predictions"),
						downloadButton("predictions_summary_docx_download", label = "Download Predictions Report - DOCX Format")
						# downloadButton("predictions_summary_pdf_download", label = "Download Predictions Summary - PDF Format")
					)
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
	
	#Hide download buttons until after model has run
	shinyjs::hide("results")
	
	global_data <- reactiveValues(RNAseq = NULL,
																model_predictions = NULL,
																model_id = NULL,
																model_pred_summary = NULL)
	
	##############################################################################
	# RNAseq Input Processing
	##############################################################################
	
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
	
	##############################################################################
	# Model Running
	##############################################################################
	
	observeEvent(global_data$RNAseq, {
		global_data$model_id = substr(digest(global_data$RNAseq), 1, 6)
	})
	
	observeEvent(global_data$model_id, {
		run_model()
	})
	
	observeEvent(global_data$model_predictions, {
		global_data$model_pred_summary = global_data$model_predictions %>% 
			group_by(drug) %>%
			summarise(mean_via = mean(predicted_viability),
								CCLE_diff = abs(mean(predicted_viability - mean_via)),
								range_via = max(predicted_viability) - min(predicted_viability))
		
		global_data$model_pred_with_CCLE = global_data$model_predictions %>%
			left_join(CCLE_preds)
	})
	
	output$RNAseq_qc_text <- renderText({
		if (is.null(global_data$model_predictions)) return()
		
		if (dim(global_data$RNAseq)[1] == 110) {
			return(paste0("Your dataset contains ", dim(global_data$RNAseq)[1], ' of the 110 genes included in the model.'))
		} else {
			return(paste0("Your dataset contains ", dim(global_data$RNAseq)[1], ' of the 110 genes included in the model. 
					 For every gene missing, the average value from the original model data has been substituted.'))
		}
	})
	
	run_model <- reactive({
		if (is.null(global_data$model_id)) return()
		
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(2/3, detail = "Loading Model/Making Predictions")
		
		global_data$model_predictions = make_predictions(global_data$RNAseq)
		
		progress$close()
		
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(3/3, detail = "Building Results Report")
		
		render('build_inhibitor_overview.Rmd', 
					 output_file = here('www/',paste0("kinase_inhibitor_summary_",global_data$model_id,".docx")), 
					 params = list(predictions = global_data$model_predictions, RNAseq_data = global_data$RNAseq, model_id = global_data$model_id))
		
		shinyjs::show("results")
	})
	
	##############################################################################
	# Model Prediction Plotting
	##############################################################################	
	
	output$minimal_eff_preds <- renderPlot({
		if (is.null(global_data$model_pred_with_CCLE)) return()
		
		low_eff_drugs = global_data$model_pred_summary %>% arrange(desc(mean_via)) %>% slice(1:5) %>% pull(drug)
		
		global_data$model_pred_with_CCLE %>%
			filter(drug %in% low_eff_drugs) %>%
			mutate(drug = fct_relevel(drug, low_eff_drugs)) %>%
			plot_pred_set() + theme(text = element_text(size=16))
	})
	
	output$high_eff_preds <- renderPlot({
		if (is.null(global_data$model_pred_with_CCLE)) return()
		
		high_eff_drugs = global_data$model_pred_summary %>% arrange(mean_via) %>% slice(1:5) %>% pull(drug)
		
		global_data$model_pred_with_CCLE %>%
			filter(drug %in% high_eff_drugs) %>%
			mutate(drug = fct_relevel(drug, high_eff_drugs)) %>%
			plot_pred_set() + theme(text = element_text(size=16))
	})
	
	output$high_range_preds <- renderPlot({
		if (is.null(global_data$model_pred_with_CCLE)) return()
		
		high_range_drugs = global_data$model_pred_summary %>% arrange(desc(range_via)) %>% slice(1:5) %>% pull(drug)
		
		global_data$model_pred_with_CCLE %>%
			filter(drug %in% high_range_drugs) %>%
			mutate(drug = fct_relevel(drug, high_range_drugs)) %>%
			plot_pred_set() + theme(text = element_text(size=16))
	})
	
	output$ccle_diff_preds <- renderPlot({
		if (is.null(global_data$model_pred_with_CCLE)) return()
		
		ccle_diff_drugs = global_data$model_pred_summary %>% arrange(desc(CCLE_diff), desc(range_via)) %>% slice(1:5) %>% pull(drug)
		
		global_data$model_pred_with_CCLE %>%
			filter(drug %in% ccle_diff_drugs) %>%
			mutate(drug = fct_relevel(drug, ccle_diff_drugs)) %>%
			plot_pred_set() + theme(text = element_text(size=16))
	})
	
	##############################################################################
	# Model Results Download
	##############################################################################
	
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
}

# Run the application 
shinyApp(ui = ui, server = server)