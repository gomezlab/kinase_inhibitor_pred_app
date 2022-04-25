convert_salmon_to_HGNC_TPM <- function(transript_data) {
	
	transcript_data_TPM = c()
	
	#look for ENST in the name column
	if (mean(str_detect(transript_data$Name,"ENST")) > 0.95) {
		#Check for the version dot in the Name column, if there, remove it with separate
		if (any(str_detect(transript_data$Name, "\\."))) {
			transript_data = transript_data %>% separate(Name, into = c("Name",NA), sep = "\\.")
		}
		
		ENST_to_hgnc = read_csv(here('data/model_expression_genes_ENST.csv'))
		
		transcript_data_TPM = transript_data %>% 
			filter(Name %in% ENST_to_hgnc$ensembl_transcript_id) %>% 
			left_join(ENST_to_hgnc, by = c("Name"="ensembl_transcript_id")) %>% 
			group_by(hgnc_symbol) %>% 
			summarise(TPM = sum(TPM))
	} else if (mean(str_detect(transript_data$Name, "^NM_")) > 0.5) {
		RefSeq_to_hgnc = read_csv(here('data/model_expression_genes_refseq.csv'))
		
		transcript_data_TPM = transript_data %>% 
			filter(Name %in% RefSeq_to_hgnc$refseq_mrna) %>% 
			left_join(RefSeq_to_hgnc, by = c("Name"="refseq_mrna")) %>% 
			group_by(hgnc_symbol) %>% 
			summarise(TPM = sum(TPM))
	}
	
	return(transcript_data_TPM)
}

plot_pred_set <- function(data_set) {
	ggplot(data_set, aes(x=log10(concentration_M),y=mean_via)) +
		geom_ribbon(aes(ymin=lower_bound,ymax=upper_bound), alpha=0.25) +
		geom_line(alpha=0.75) +
		geom_line(aes(x=log10(concentration_M),y=predicted_viability), color = "blue") +
		geom_point(aes(x=log10(concentration_M),y=predicted_viability), color = "blue") +
		ylim(c(0,NA)) +
		labs(x = "Compound Concentration (Log 10 M)",y = "Predicted Viability") +
		BerginskiRMisc::theme_berginski() +
		facet_grid(cols = vars(drug))
}

make_predictions <- function(processed_RNAseq,prediction_model) {
	average_exp_vals = read_rds(here('data/average_model_exp_vals.rds'))
	
	klaeger_wide = read_rds(here('data/klaeger_wide.rds')) %>%
		filter(concentration_M != 0)
	
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
	
	model_predictions = predict(prediction_model, 
															model_data %>% mutate(klaeger_conc = NA, imputed_viability = NA, depmap_id = NA))
	model_predictions$drug = model_data$drug
	model_predictions$concentration_M = model_data$concentration_M
	
	model_predictions = model_predictions %>%
		rename(predicted_viability = .pred) %>%
		mutate(predicted_viability = signif(predicted_viability,3)) %>%
		select(drug, concentration_M, everything())
	
	rm(klaeger_wide)
	rm(average_exp_vals)
	gc()
	
	return(model_predictions)
}
