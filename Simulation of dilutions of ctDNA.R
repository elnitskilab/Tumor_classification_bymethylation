######################################################
######################################################
######################################################
# An R script that applies logistic regression models 
# to pooled methylation data for several cancer types at 
# three probes with 10 fold cross validation to find a 
# common threshold and applies the threshold to the same
# data to test its accuracy.


# Author: Karen Funderburk
# Created: June 19, 2019
# Modified: May 4, 2022

######################################################

 rm(list = ls(all = TRUE))  # resets R to fresh

 library('pROC')
 library('ResourceSelection')
 library('caret')

##################################################
# Reading in input data


types = c("BLCA", "BRCA", "COAD")

probes_list = c('cg03502002', 'cg14861089', 'cg21790626')


load("Blood_EXAMPLE.RData")
load("TCGA_EXAMPLE_FILE.RData")

# dilution percentages 1-100%
frac = c(0.01,0.1,0.25,0.5, 0.75, 0.9, 0.99, 0.999, 1)

#setting up cross validation
train.control <- trainControl(method = "cv", number = 10)

results = data.frame()
tumor_df = data.frame()

#constructing methylation data frame with TCGA data for selected tumor types
 for (type in types){
 	
   tmp_df = TCGA_EXAMPLE[TCGA_EXAMPLE$type==paste0(type),]#selecting tumor type
 	tmp_df = tmp_df[tmp_df$sample_type==1,] #Selecting only tumor samples
 	colnames(tmp_df)[1] = "tum_status"
  tmp_df$tum_status = factor(1)
  tumor_df = rbind(tumor_df, tmp_df)
 }
     	
#constructing normal (blood) data
set.seed(2023)      	
 blood_df = data.frame(Blood_EXAMPLE[sample(nrow(Blood_EXAMPLE), size = nrow(tumor_df), replace = T),])
 blood_df$tum_status = factor(0)
 blood_df$type = 'blood'

for (tdna_frac in frac){
  
  #adjusting tumor methylation level using blood data (simulated dilutions)
  adj_tumor_df = tumor_df
  adj_tumor_df[, probes_list] = (tdna_frac * tumor_df[, probes_list]) + 
  								((1-tdna_frac) * blood_df[, probes_list])
  
  methylation_df = rbind(blood_df, adj_tumor_df)
  
		
	for (probe in probes_list){

		# gathering training data
		TrainingData <- data.frame(cpg = methylation_df[,probe], 
								   status = methylation_df$tum_status)
		
		#building logistic regression model						   
		single_lr = train(status ~ cpg, 
						  data = TrainingData, 
						  method = 'glm', 
						  trControl = train.control)

		#generating roc object from regression probabilities 
		single_roc = roc(TrainingData$status ~ single_lr$finalModel$fitted.values, 
						 quiet = T, 
						 plot = F)
					 
		#storing best threshold, specificity, sensitivity, and AUC. Youden's statistic = max(sens + spec)
		single_results = coords(single_roc, 
			                    x = 'best', 
			                    best.method = 'youden', 
			                    ret = c('threshold','sens', 'spec','youden',
			                        	'tp','tpr','fp', 'fpr','tn', 'tnr','fn', 'fnr',
			                        	'precision','recall','accuracy'), 
			                    transpose = F)
			
		#choosing threshold with most tp if there are multiple 'bests'
		if (length(single_results$threshold) >1){single_results = single_results[single_results$tp == max(single_results$tp),]}
				

		# running validations on each type 
		threshold = single_results$threshold
		
		tum_validations = numeric()
		
		for (type in c(types, 'blood')){
					
			if(type == 'blood'){
			
				TestingData <- data.frame(cpg = blood_df[,probe])
				predictions = predict(single_lr, newdata = TestingData, type = 'prob')['1']

				tn = sum(predictions < threshold)
	 	    	fp = sum(predictions > threshold)
	 	    	specificity = tn/(tn + fp)
	 	    	
	 	    }else{
	 	    
				TestingData <- data.frame(cpg = adj_tumor_df[adj_tumor_df$type == type,probe])
				predictions = predict(single_lr, newdata = TestingData, type = 'prob')['1']

				tp = sum(predictions > threshold)
	     		fn = sum(predictions < threshold)
	     		sensitivity = tp/(tp + fn)
	     		
	     		tmp = c(tp = tp, fn = fn, sens = sensitivity)
	     		names(tmp) = paste0(type, '_', names(tmp))
				tum_validations = c(tum_validations, tmp)
			}}

 		#storing validation results
 		validations = data.frame(cbind( blood_tn = tn, 
 										blood_fp = fp, 
 										blood_spec= specificity, 
 										t(tum_validations)))

		
		
		
		#storing all results
		results = rbind(results, 
					    data.frame(row.names = NULL,
				                   cbind(probe, 
				                   		 combo_type = 'Single',
				                         tdna_frac, 
				                         single_results,
				                         auc = single_roc$auc,
				                         cv_accuracy = single_lr$results$Accuracy, 
				                         cv_kappa = single_lr$results$Kappa, 
				                         validations))
								)		

		
		
	}
	



	#pairwise probe roc calculations

	combos = combn(probes_list, 2)

	for(i in 1:ncol(combos)){
	
		probe1 = combos[1,i]
		probe2 = combos[2,i]
		

        TrainingData <- data.frame(cpg1 = methylation_df[,probe1], 
        							   cpg2 = methylation_df[,probe2], 
        							   status = methylation_df$tum_status)
        							   
        							   
        pairwise_lr = train(status ~ cpg1 + cpg2, 
        						data = TrainingData, 
        						method = 'glm', 
        						trControl = train.control)

	    #generating roc object for pairwise combinations of probes
        pw_roc = roc(TrainingData$status ~ pairwise_lr$finalModel$fitted.values, 
                     		 quiet = T, 
	                     	 plot = F) 
	        
		pw_results = coords(pw_roc, 
							x = 'best', 
								best.method = 'youden', 
								ret = c('threshold','sens', 'spec','youden',
											'tp','tpr','fp', 'fpr','tn', 'tnr','fn', 'fnr',
											'precision','recall','accuracy'), 
								transpose = F)   
	        
		#selects result with the most true positives if multiple 'bests' exist                
		if (length(pw_results$threshold) >1){pw_results = pw_results[pw_results$tp == max(pw_results$tp),]}



		# running validations on each type
		threshold = pw_results$threshold
		tum_validations = numeric()
				
		for (type in c(types, 'blood')){
					
			if(type == 'blood'){
			
				TestingData <- data.frame(cpg1 = blood_df[,probe1],
										  cpg2 = blood_df[,probe2])
				predictions = predict(pairwise_lr, newdata = TestingData, type = 'prob')['1']

				tn = sum(predictions < threshold)
	 	    	fp = sum(predictions > threshold)
	 	    	specificity = tn/(tn + fp)
	 	    	
	 	    }else{
	 	    
				TestingData <- data.frame(cpg1 = adj_tumor_df[adj_tumor_df$type == type,probe1],
										  cpg2 = adj_tumor_df[adj_tumor_df$type == type,probe2])
				predictions = predict(pairwise_lr, newdata = TestingData, type = 'prob')['1']

				tp = sum(predictions > threshold)
	     		fn = sum(predictions < threshold)
	     		sensitivity = tp/(tp + fn)
	     		
	     		tmp = c(tp = tp, fn = fn, sens = sensitivity)
	     		names(tmp) = paste0(type, '_', names(tmp))
				tum_validations = c(tum_validations, tmp)
			}}

 		#storing validation results
 		validations = data.frame(cbind( blood_tn = tn, 
 										blood_fp = fp, 
 										blood_spec= specificity, 
 										t(tum_validations)))

		
		
		
		#storing all results
		results = rbind(results, 
					    data.frame(row.names = NULL,
				                   cbind(probe = paste(probe1, 'vs', probe2), 
				                   		 combo_type = 'Pairwise',
				                         tdna_frac, 
				                         pw_results,
				                         auc = pw_roc$auc,
				                         cv_accuracy = pairwise_lr$results$Accuracy, 
				                         cv_kappa = pairwise_lr$results$Kappa, 
				                         validations))
								)		


		
	}
	

plot.new()

	#triplicate probe roc calculations

	combos = combn(probes_list, 3)

	for (i in 1:ncol(combos)){
		probe1 = combos[1,i]
		probe2 = combos[2,i]
		probe3 = combos[3,i]

		TrainingData <- data.frame(cpg1 = methylation_df[,probe1], 
									   cpg2 = methylation_df[,probe2], 
									   cpg3 = methylation_df[,probe3], 
									   status = methylation_df$tum_status)
									   
		tripl_lr = train(status ~ cpg1 + cpg2 + cpg3, 
							 data = TrainingData, 
							 method = 'glm', 
							 trControl = train.control)

		#generates roc object 
		tp_roc = roc(TrainingData$status ~ tripl_lr$finalModel$fitted.values, 
					 quiet = T, 
					 plot = F) 

				
		tp_results = coords(tp_roc, 
							x = 'best', 
							best.method = 'youden', 
							ret = c('threshold','sens', 'spec','youden',
									'tp','tpr','fp', 'fpr','tn', 'tnr','fn', 'fnr',
									'precision','recall','accuracy'), 
							transpose = F)
									
		#selects result with the most true positives if multiple 'bests' exist    
		if (length(tp_results$threshold) >1){tp_results = tp_results[tp_results$tp == max(tp_results$tp),]}
			

		# running validations on each type
		threshold = tp_results$threshold
		tum_validations = numeric()
		
		for (type in c(types, 'blood')){
					
			if(type == 'blood'){
			
				TestingData <- data.frame(cpg1 = blood_df[,probe1],
										  cpg2 = blood_df[,probe2],
										  cpg3 = blood_df[,probe3])
				predictions = predict(tripl_lr, newdata = TestingData, type = 'prob')['1']

				tn = sum(predictions < threshold)
	 	    	fp = sum(predictions > threshold)
	 	    	specificity = tn/(tn + fp)
	 	    	
	 	    }else{
	 	    
				TestingData <- data.frame(cpg1 = adj_tumor_df[adj_tumor_df$type == type,probe1],
										  cpg2 = adj_tumor_df[adj_tumor_df$type == type,probe2],
										  cpg3 = adj_tumor_df[adj_tumor_df$type == type,probe3])
				predictions = predict(tripl_lr, newdata = TestingData, type = 'prob')['1']

				tp = sum(predictions > threshold)
	     		fn = sum(predictions < threshold)
	     		sensitivity = tp/(tp + fn)
	     		
	     		tmp = c(tp = tp, fn = fn, sens = sensitivity)
	     		names(tmp) = paste0(type, '_', names(tmp))
				tum_validations = c(tum_validations, tmp)
			}}

 		#storing validation results
 		validations = data.frame(cbind( blood_tn = tn, 
 										blood_fp = fp, 
 										blood_spec= specificity, 
 										t(tum_validations)))

		
		
		
		#storing all results
		results = rbind(results, 
					    data.frame(row.names = NULL,
				                   cbind(probe = paste(probe1, 'vs', probe2, 'vs', probe3), 
				                   		 combo_type = 'Triplicate',
				                         tdna_frac, 
				                         tp_results,
				                         auc = tp_roc$auc,
				                         cv_accuracy = tripl_lr$results$Accuracy, 
				                         cv_kappa = tripl_lr$results$Kappa, 
				                         validations))
								)		

	
		
	}


}



write.table(results, 'SIMULATED DILUTION_EXAMPLE FILE.txt', quote = FALSE, row.names = FALSE, sep = '\t')


