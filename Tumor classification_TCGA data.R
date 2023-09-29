######################################################
######################################################
######################################################
# An R script that demonstrates how methylation data from TCGA (beta values) can be used in logistic regression models to classify tumor
# versus normal, followed by validation with inclusion of independent datasets. 

# Author: Karen Funderburk 
# Created: June 19, 2019
# Modified: November 1, 2022

######################################################

 library('pROC')
 library('car')
 library('ResourceSelection')

##################################################

types = c("BLCA", "BRCA", "COAD")

probes_list = c('cg03502002', 'cg14861089', 'cg21790626')


# reading in independent methylation datasets

ind_types = c('ind_BRCA')

load("TCGA_EXAMPLE_FILE.RData")
load("validation_EXAMPLE_FILE.RData")
ind_df = validation_EXAMPLE

results = data.frame()

for (type in types){
  
  df=TCGA_EXAMPLE[TCGA_EXAMPLE$type==paste0(type),]#selects TCGA data for specified type
  
	for (probe in probes_list){
	  
	  #SINGLE PROBE ROC CALCULATIONS
	  
	  # running logistic regression
	  training_data = data.frame(cpg = df[[probe]],
	                             status = df[['sample_type']])
	  
	  single_lr = glm(status ~ cpg, 
	                  data = training_data,
	                  family = 'binomial')
	  
	  #feeding in fitted values to ROC function
	  single_roc = roc(df$sample_type ~ single_lr$fitted.values, 
	                   quiet = T, 
	                   plot = F)
	              
	  
	  #stores best threshold, specificity, and sensitivity. Youden's statistic = max(sens + spec)
	  single_results = coords(single_roc, 
	                          x = 'best', 
	                          best.method = 'youden', 
	                          ret = c('threshold','sens', 'spec','youden',
	                                  'tp','tpr','fp', 'fpr','tn', 'tnr','fn', 'fnr',
	                                  'precision','recall','accuracy'), 
	                          transpose = 'F')
	  
	  #choosing threshold with most tp if there are multiple 'bests'
	  if (length(single_results$threshold) >1){single_results = single_results[single_results$tp == max(single_results$tp),]}
	  
	  #validations
    validation_df = data.frame(0)
    
    for (vtype in ind_types){ 
  	  

  	    new_data = ind_df[ind_df$type == vtype,]
  	    test_data = data.frame(cpg = new_data[[probe]],
  	                           status = new_data$tum_status)
      
      predictions = predict.glm(single_lr, newdata = test_data, type = 'response')
  	  
  	  tn = sum(predictions[test_data$status == 0] < single_results$threshold, na.rm = T)
  	  tp = sum(predictions[test_data$status == 1] > single_results$threshold, na.rm = T)
  	  fn = sum(predictions[test_data$status == 1] < single_results$threshold, na.rm = T)
  	  fp = sum(predictions[test_data$status == 0] > single_results$threshold, na.rm = T)
  	  specificity = tn/(tn + fp)
  	  sensitivity = tp/(tp + fn)
  	  
  	  tmp = c(tp = tp, tn = tn, fn = fn, fp = fp, sens = sensitivity, spec = specificity)
  	  names(tmp) = paste0(vtype, '_', names(tmp))
  	  validation_df=cbind(validation_df, t(tmp))
  	}  
	  
	  results = rbind(results, 
	                  data.frame(row.names = NULL,
	                             cbind(probe = probe, 
	                                   type = type,
	                                   combo_type = 'Single',
	                                   single_results,
	                                   auc = single_roc$auc,
	                                   validation_df[,-1])))
	  
	}


  
  
  
	#PAIRWISE PROBE ROC CALCULATIONS

	combos = combn(probes_list, 2)
	
	for (i in 1:ncol(combos)){
	
		probe1 = combos[1,i]
		probe2 = combos[2,i]

		#uses logistic regression estimates to get linear combination of 2 predictors 
		training_data = data.frame(cpg1 = df[[probe1]],
		                           cpg2 = df[[probe2]],
		                           status = df$sample_type)
		
		pairwise_lr = glm(status ~ cpg1 + cpg2, 
		                  data = training_data,
		                  family = 'binomial')

		#generates roc object for combinations of probes
		pw_roc = roc(df$sample_type ~ pairwise_lr$fitted.values, 
		             quiet = T, 
		             plot = F)
		
		#stores best threshold, specificity, and sensitivity. Youden's statistic = max(sens + spec)
		pw_results = coords(pw_roc, 
		                    x = 'best', 
		                    best.method = 'youden', 
		                    ret = c('threshold','sens', 'spec','youden',
		                            'tp','tpr','fp', 'fpr','tn', 'tnr','fn', 'fnr',
		                            'precision','recall','accuracy'),
		                    transpose = 'F')	
		
		#choosing threshold with most tp if there are multiple 'bests'
		if (length(pw_results$threshold) >1){pw_results = pw_results[pw_results$tp == max(pw_results$tp),]}
		
		#validations
		validation_df = data.frame(0)
		for (vtype in ind_types){ 
		  
		    new_data = ind_df[ind_df$type == vtype,]
		    test_data = data.frame(cpg1 = new_data[[probe1]],
		                           cpg2 = new_data[[probe2]],
		                           status = new_data$tum_status)
		

		  predictions = predict.glm(pairwise_lr, 
		                        newdata = test_data, 
		                        type = 'response')
		  
		  tn = sum(predictions[test_data$status == 0] < pw_results$threshold, na.rm = T)
		  tp = sum(predictions[test_data$status == 1] > pw_results$threshold, na.rm = T)
		  fn = sum(predictions[test_data$status == 1] < pw_results$threshold, na.rm = T)
		  fp = sum(predictions[test_data$status == 0] > pw_results$threshold, na.rm = T)
		  specificity = tn/(tn + fp)
		  sensitivity = tp/(tp + fn)
		  
		  tmp = c(tp = tp, tn = tn, fn = fn, fp = fp, sens = sensitivity, spec = specificity)
		  names(tmp) = paste0(vtype, '_', names(tmp))
		  validation_df=cbind(validation_df, t(tmp))
		}
		
		
		#storing results
		results = rbind(results, 
		                data.frame(row.names = NULL,
		                           cbind(probe = paste0(probe1, '_vs_', probe2), 
		                                 type = type,
		                                 combo_type = 'Pairwise',
		                                 pw_results,
		                                 auc = pw_roc$auc,
		                                 validation_df[,-1])))
	} 		


	#TRIPLICATE PROBE ROC CALCULATIONS

	combos = combn(probes_list, 3)

	for (i in 1:ncol(combos)){
		probe1 = combos[1,i]
		probe2 = combos[2,i]
		probe3 = combos[3,i]
	
		#logistic regression of 3 predictors to get linear combination for use in roc
    training_data = data.frame(cpg1 = df[[probe1]],
                               cpg2 = df[[probe2]],
                               cpg3 = df[[probe3]],
                               status = df$sample_type)
    
		tripl_lr = glm(status ~ cpg1 + cpg2 + cpg3, 
		               data = training_data,
		               family = 'binomial')	
		
		#generates roc object for triplicate combinations of probes
		tp_roc = roc(df$sample_type ~ tripl_lr$fitted.values, 
		             quiet = T, 
		             plot = F)
		
		#stores threshold, specificity, and sensitivity
		tp_results = coords(tp_roc, 
		                    x = 'best', 
		                    best.method = 'youden', 
		                    ret = c('threshold','sens', 'spec','youden',
		                            'tp','tpr','fp', 'fpr','tn', 'tnr','fn', 'fnr',
		                            'precision','recall','accuracy'),
		                    transpose = 'F')	
		
		#choosing threshold with most tp if there are multiple 'bests'
		if (length(tp_results$threshold) >1){tp_results = tp_results[tp_results$tp == max(tp_results$tp),]}
		
		
		#validations
		validation_df = data.frame(0)
		
		for (vtype in ind_types){ 
		  
		    new_data = ind_df[ind_df$type == vtype,]
		    test_data = data.frame(cpg1 = new_data[[probe1]],
		                           cpg2 = new_data[[probe2]],
		                           cpg3 = new_data[[probe3]],
		                           status = new_data$tum_status)

		  
		  predictions = predict.glm(tripl_lr, 
		                        newdata = test_data, 
		                        type = 'response')
		  
		  tn = sum(predictions[test_data$status == 0] < tp_results$threshold, na.rm = T)
		  tp = sum(predictions[test_data$status == 1] > tp_results$threshold, na.rm = T)
		  fn = sum(predictions[test_data$status == 1] < tp_results$threshold, na.rm = T)
		  fp = sum(predictions[test_data$status == 0] > tp_results$threshold, na.rm = T)
		  specificity = tn/(tn + fp)
		  sensitivity = tp/(tp + fn)
		  
		  tmp = c(tp = tp, tn = tn, fn = fn, fp = fp, sens = sensitivity, spec = specificity)
		  names(tmp) = paste0(vtype, '_', names(tmp))
		  validation_df=cbind(validation_df, t(tmp))
		}  
		
		
		#storing results
		results = rbind(results, 
		                data.frame(row.names = NULL,
		                           cbind(probe = paste0(probe1, 'vs', probe2, 'vs', probe3), 
		                                 type = type,
		                                 combo_type = 'Triplicate',
		                                 tp_results,
		                                 auc = tp_roc$auc,
		                                 validation_df[,-1])))
		
		
	}
	

}


write.csv(results, file = 'ANALYSIS_OF_TCGA_DATA_EXAMPLE_FILE.csv', row.names = F)


