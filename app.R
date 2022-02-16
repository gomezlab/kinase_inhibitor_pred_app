#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(dqshiny)
library(here)
library(tidymodels)
library(reactable)
library(shinyjs)

options(shiny.maxRequestSize=30*1024^2)

all_geo_archs_ids = read_rds(here('data/ARCHS_GEO_IDs.rds'))

convert_salmon_to_HGNC_TPM <- function(transript_data) {
	
	#look for ENST in the name column
	if (mean(str_detect(transript_data$Name,"ENST")) > 0.95) {
		#Check for the version dot in the Name column, if there, remove it with separate
		if (any(str_detect(transript_data$Name, "\\."))) {
			transript_data = transript_data %>% separate(Name, into = c("Name",NA), sep = "\\.")
		}
		
		ENST_to_hgnc = read_csv(here('data/model_expression_genes.csv'))
		
		transript_data = transript_data %>% 
			filter(Name %in% ENST_to_hgnc$ensembl_transcript_id) %>% 
			left_join(ENST_to_hgnc, by = c("Name"="ensembl_transcript_id")) %>% 
			group_by(hgnc_symbol) %>% 
			summarise(TPM = sum(TPM))
		
		return(transript_data)
	} else {
		
	}
}

make_predictions <- function(processed_RNAseq) {
	progress <- shiny::Progress$new()
	# Make sure it closes when we exit this reactive, even if there's an error
	on.exit(progress$close())

	average_exp_vals = read_rds(here('data/average_model_exp_vals.rds'))
	
	klaeger_wide = read_rds(here('data/klaeger_wide.rds')) %>%
		filter(concentration_M != 0)
	
	rand_forest_model = read_rds(here('data/final_model_500feat_100trees.rds'))
	
	model_data = processed_RNAseq %>% 
		mutate(model_feature = paste0("exp_",hgnc_symbol), 
					 trans_TPM = log2(TPM + 1)) %>% 
		select(-hgnc_symbol,-TPM)

	missing_RNAseq_vals = average_exp_vals %>% 
		filter(! model_feature %in% model_data$model_feature)
	
	model_data = model_data %>%
		bind_rows(missing_RNAseq_vals) %>%
		pivot_wider(names_from = model_feature,values_from = trans_TPM) %>%
		slice(rep(1:n(), each = dim(klaeger_wide)[1])) %>%
		bind_cols(klaeger_wide)
	
	progress$inc(3/3, detail = "Making Model Predictions")
	
	model_predictions = predict(rand_forest_model, model_data %>% mutate(klaeger_conc = NA, imputed_viability = NA, depmap_id = NA))
	model_predictions$drug = model_data$drug
	model_predictions$concentration_M = model_data$concentration_M
	
	model_predictions = model_predictions %>%
		rename(predicted_viability = .pred) %>%
		select(drug, concentration_M, everything())
	
	rm(rand_forest_model)
	rm(klaeger_wide)
	rm(average_exp_vals)
	gc()
	
	return(model_predictions)
}

# Define UI for application that draws a histogram
ui <- fluidPage(
	
	shinyjs::useShinyjs(), 
	
	# Application title
	titlePanel("Kinase Inhibitor Cell Viability Prediction"),
	
	tags$p(""),
	
	# Sidebar with a slider input for number of bins 
	sidebarLayout(
		sidebarPanel(
			fileInput("RNAseq_file", "Choose RNAseq Results File",
								multiple = TRUE),
			tags$hr(),
			
			autocomplete_input("GEO_ARCHS_ID", 
												 h2("GEO ARCHS ID"), 
												 all_geo_archs_ids,
												 placeholder = "Start Typing to Find Your GEO ID",
												 max_options = 100),
			
			tags$hr(),
			
		),

		# Show a plot of the generated distribution
		mainPanel(
			tags$h1("RNAseq Data Checks"),
			
			textOutput("RNAseq_qc_text"),
			fluidRow(
				column(6,tableOutput("RNAseq_sample")),
				
				column(6,tableOutput("prediction_sample"))
			),
			
			downloadButton("model_predictions_download", label = "Download Model Predictions")
		)
	)
)

# Define server logic required to draw a histogram
server <- function(input, output) {
	
	
	RNAseq_data <- reactive({
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())

		progress$inc(1/3, detail = "Recieving Data")
		
		read_delim(input$RNAseq_file$datapath, delim = "\t") %>%
			convert_salmon_to_HGNC_TPM()
	})
	
	model_predictions <- reactive({
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(3/3, detail = "Making Model Predictions")
		
		prediction_results = make_predictions(RNAseq_data())
		
		shinyjs::show("model_predictions_download")
		
		prediction_results
	})
	
	output$RNAseq_qc_text <- renderText({
		req(input$RNAseq_file)
		
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(2/3, detail = "Processing Data")
		
		paste0("Your file contains ", dim(RNAseq_data())[1], '/110 genes.')
		
	})
	
	output$RNAseq_sample <- renderTable({
		
		req(input$RNAseq_file)
		
		
		return(head(RNAseq_data()))
	})
	
	output$prediction_sample <- renderTable({
		
		req(input$RNAseq_file)

		return(head(model_predictions()))
	})
	
	output$model_predictions_download <- downloadHandler(
		filename = function() {
			paste0("kinase_inhbitor_model_predictions.csv")
		}, 
		content = function(file) {
			write_csv(model_predictions(), file)
		})
	
	shinyjs::hide("model_predictions_download")
}

# Run the application 
shinyApp(ui = ui, server = server)
