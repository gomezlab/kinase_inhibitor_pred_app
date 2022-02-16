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

options(shiny.maxRequestSize=30*1024^2)

all_geo_archs_ids = read_rds(here('data/ARCHS_GEO_IDs.rds'))

ENST_to_hgnc = read_csv(here('data/model_expression_genes.csv'))

convert_salmon_to_HGNC_TPM <- function(transript_data) {
	#look for ENST in the name column
	if (mean(str_detect(transript_data$Name,"ENST")) > 0.95) {
		if (any(str_detect(transcript_test$Name, "\\."))) {
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
        	
        	tableOutput("contents")
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
}

# Run the application 
shinyApp(ui = ui, server = server)
