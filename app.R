library(shiny)
library(tidyverse)
library(dqshiny)
library(here)
library(tidymodels)
library(reactable)
library(shinyjs)
library(digest)

options(shiny.maxRequestSize=30*1024^2)

dir.create(here('data/uploads'), showWarnings = F)

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
	
	progress$inc(2/3, message = NULL, detail = "Loading Model")
	
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
		mutate(predicted_viability = signif(predicted_viability,3)) %>%
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
	
	tags$p("Hello and welcome to the kinase inhibitor cell viability prediction web server. This website is a companion 
				 to a forthcoming publication concerning prediction of cell viability after exposure to kinase inhibitors. 
				 The primary model developed in this paper uses RNAseq gene expression and the kinase inhibitor target profiles
				 to make these predictions. We've made this server available to allow interested biologists to submit gene expression
				 data they have gathered where there is some interest in how a set of kinase inhbitors would affect their model 
				 system."),
	
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
server <- function(input, output) {
	RNAseq_data <- reactive({
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(1/3, detail = "Processing Data")
		
		hide("instructions")
		
		TPM_data = read_delim(input$RNAseq_file$datapath, delim = "\t") %>%
			convert_salmon_to_HGNC_TPM()
		
		file.copy(input$RNAseq_file$datapath, here('data/uploads',paste0(substr(digest(TPM_data), 1, 6))))
		TPM_data
	})
	
	model_predictions <- reactive({
		prediction_results = make_predictions(RNAseq_data())
		
		shinyjs::show("results")
		shinyjs::show("model_predictions_download")
		
		prediction_results
	})
	
	output$RNAseq_qc_text <- renderText({
		req(input$RNAseq_file)
		
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
			
			paste0("kinase_inhbitor_model_predictions_",paste0(substr(digest(RNAseq_data()), 1, 6)),".csv")
		}, 
		content = function(file) {
			write_csv(model_predictions(), file)
		})
	
	shinyjs::hide("model_predictions_download")
	shinyjs::hide("results")
}

# Run the application 
shinyApp(ui = ui, server = server)