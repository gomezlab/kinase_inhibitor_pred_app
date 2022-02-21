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
	
	if (shiny::isRunning()) {
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		
		progress$inc(2/3, message = NULL, detail = "Loading Model")
	}
	average_exp_vals = read_rds(here('data/average_model_exp_vals.rds'))
	
	klaeger_wide = read_rds(here('data/klaeger_wide.rds')) %>%
		filter(concentration_M != 0)
	
	rand_forest_model = read_rds(here('data/final_model_500feat_100trees.rds'))
	
	if (shiny::isRunning()) {
		progress$inc(3/3, detail = "Making Model Predictions")
	}
	
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
	
	model_predictions = predict(rand_forest_model, 
															model_data %>% mutate(klaeger_conc = NA, imputed_viability = NA, depmap_id = NA))
	model_predictions$drug = model_data$drug
	model_predictions$concentration_M = model_data$concentration_M
	
	model_predictions = model_predictions %>%
		rename(predicted_viability = .pred) %>%
		mutate(predicted_viability = signif(predicted_viability,3)) %>%
		select(drug, concentration_M, everything())
	
	# rm(rand_forest_model)
	
	rm(klaeger_wide)
	rm(average_exp_vals)
	gc()
	
	return(model_predictions)
}
