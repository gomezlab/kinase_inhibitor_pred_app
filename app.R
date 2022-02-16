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

options(shiny.maxRequestSize=30*1024^2)

all_geo_archs_ids = read_rds(here('data/ARCHS_GEO_IDs.rds'))

ENST_to_hgnc = read_csv(here('data/model_expression_genes.csv'))

convert_salmon_to_HGNC_TPM <- function(transript_data) {
	#look for ENST in the name column
	if (mean(str_detect(transript_data$Name,"ENST")) > 0.95) {
		#Check for the version dot in the Name column, if there, remove it with separate
		if (any(str_detect(transript_data$Name, "\\."))) {
			transript_data = transript_data %>% separate(Name, into = c("Name",NA), sep = "\\.")
		}
		
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
	klaeger_wide = read_rds(here('data/klaeger_wide.rds'))
	
	rand_forest_model = read_rds(here('data/final_model_500feat_100trees.rds'))
	
	model_data = processed_RNAseq %>% 
		mutate(model_feature = paste0("exp_",hgnc_symbol), 
					 trans_TPM = log2(TPM + 1)) %>% 
		select(-hgnc_symbol,-TPM) %>% 
		pivot_wider(names_from = model_feature,values_from = trans_TPM) %>% 
		slice(rep(1:n(), each = dim(klaeger_wide)[1])) %>%
		bind_cols(klaeger_wide)
	
	model_predictions = predict(rand_forest_model, model_data %>% mutate(klaeger_conc = NA, imputed_viability = NA, depmap_id = NA))
	model_predictions$drug = model_data$drug
	model_predictions$concentration_M = model_data$concentration_M
	
	return(model_predictions)
}

# Define UI for application that draws a histogram
ui <- fluidPage(

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
        	
        	tableOutput("contents"),
        	
        	tableOutput("pred")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
	
	RNAseq_data <- reactive({
		read_delim(input$RNAseq_file$datapath, delim = "\t") %>%
			convert_salmon_to_HGNC_TPM()
	})

	output$RNAseq_qc_text <- renderText({
		req(input$RNAseq_file)
		
		paste0("Your file contains ", dim(RNAseq_data())[1], '/110 genes.')
		
	})
	
	output$contents <- renderTable({
		
		req(input$RNAseq_file)
		
		return(RNAseq_data())
	})
	
	output$pred <- renderTable({
		
		req(input$RNAseq_file)
		
		return(make_predictions(RNAseq_data()))
	})
}

# Run the application 
shinyApp(ui = ui, server = server)
