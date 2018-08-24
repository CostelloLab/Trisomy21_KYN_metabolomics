require(glmnet, quietly = T)
require(plyr, quietly = T)
require(RColorBrewer, quietly = T)
library(ROCR)
library(pROC)
cols = brewer.pal(11, 'Spectral')

LassoTally = function(train.dataset, pos.class.label, n.rounds = NULL, subset.train = F,
                      n.perfect = NULL, test.dataset = NULL, train.percent = .8,
                      good.acc = .9, n.folds = 5, model.note = ''){
  
  # train.dataset and test.dataset should have samples as rows and features as
  # the columns
  
  # Make sure that either n.perfect or n.rounds is defined
  if (is.null(n.rounds) & is.null(n.perfect) | !is.null(n.rounds) & !is.null(n.perfect)){
    stop('Please define either the n.rounds or n.perfect argument\n')
  }
  
  # Choose the iterator value, when to stop iterating
  if (!is.null(n.rounds)){
    iterator.max = n.rounds
  } else if (!is.null(n.perfect)){
    iterator.max = n.perfect
  }
  
  # Track these
  i = 0                         # Iterator - either the n.rounds or n.perfect
  n.models = 0                  # Number of models fit total
  n.perfect.models = 0          # Number of 100% accurate models
  perfect.features = c()        # Features (proteins) used in 100% accurate models
  perfect.train.samples = c()   # Samples in the training set when the model was 100% accurate
  perfect.test.samples = c()    # Samples in the testing set when the model was 100% accurate
  good.features = c()           # Features used in >= good.acc % accurate models
  all.train.samples = c()       # Samples in the training set
  all.test.samples = c()        # Samples in the testing set
  all.acc = c()                 # All accuracy values acheived
  best.acc = 0
  mod.specificities = c()
  mod.sensitivities = c()
  roc.auc = c()
  
  while (i < iterator.max){
    
    # Define training and testing sets
    if (is.null(test.dataset)){
      data.for.model = train.dataset
      sample.ids = row.names(data.for.model)
      inTrain = sample.ids %in% sample(sample.ids, floor(train.percent*length(sample.ids)))
      train.data = data.for.model[inTrain,]
      test.data = data.for.model[!inTrain,]
      
      # Format for glmnet - exclude SampleGroup from the feature matrices
      feature.matrix = data.matrix(train.data[,-ncol(train.data)])
      response.vector = factor(train.data[,ncol(train.data)])
      testing.matrix = data.matrix(test.data[,-ncol(test.data)])
      testing.response = factor(test.data[,ncol(test.data)])
      
    } else if (!is.null(test.dataset)){
      if (subset.train){
        inTrain = sample(row.names(train.dataset),
                         floor(train.percent*nrow(train.dataset)))
        
        feature.matrix = data.matrix(train.dataset[inTrain,-ncol(train.dataset)])
        response.vector = factor(train.dataset[inTrain,ncol(train.dataset)])
        testing.matrix = data.matrix(test.dataset[,-ncol(test.dataset)])
        testing.response = factor(test.dataset[,ncol(test.dataset)])
      } else{
        feature.matrix = data.matrix(train.dataset[,-ncol(train.dataset)])
        response.vector = factor(train.dataset[,ncol(train.dataset)])
        testing.matrix = data.matrix(test.dataset[,-ncol(test.dataset)])
        testing.response = factor(test.dataset[,ncol(test.dataset)])
      }
    }
    
    #cat(table(response.vector), '\n')
    
    # Don't use if a class has only 1 observation
    if (length(table(response.vector)) != 2 | 
        any(as.numeric(table(response.vector)) < 2)){
      #cat('wrong!')
      next
    }
    
    # Train lasso model
    n.models = n.models + 1
    lasso.fit = cv.glmnet(x = feature.matrix,
                          y = response.vector,
                          family = 'binomial',
                          nfolds = n.folds,
#                          type.measure = 'auc',
                          intercept = F,
                          alpha = 1)
    
    # Store feature coefficients for each model
    coeffs = coef(lasso.fit, s='lambda.min', exact=TRUE)
    idx = which(coeffs[,1] !=0)
    variables = row.names(coeffs)[idx]
    df = data.frame(Original_ID = row.names(coeffs)[idx],
                    coeffs@x)
    
    # Predict on testing
    lasso.preds = predict(lasso.fit, newx = testing.matrix, type = 'class')
    lasso.preds2 = predict(lasso.fit, newx = testing.matrix)
    pr = roc(as.character(testing.response), lasso.preds2[,1])
    mod.sensitivities = c(mod.sensitivities, pr$sensitivities)
    mod.specificities = c(mod.specificities, pr$specificities)
    roc.auc = c(roc.auc, pr$auc)

    metrics = ReportPerfMetrics(lasso.preds[,ncol(lasso.preds)], 
                                testing.response, pos.class = pos.class.label)
    all.acc = c(all.acc, metrics$Accuracy)
    
    # Store all model performance metrics
    if (exists('all.perf.metrics')){
      suppressMessages(all.perf.metrics <- rbind(all.perf.metrics, as.data.frame(metrics)))
    } else{
      all.perf.metrics = as.data.frame(metrics)
    }
    
    # Store all model coefficients
    names(df)[2] = get('n.models')
    if (exists('all.coeffs.df')){
      suppressMessages(all.coeffs.df <- join(all.coeffs.df, df, type = 'full'))
    } else{
      all.coeffs.df = df
    }
    
    # Good coefficients
    if (metrics$Accuracy >= good.acc){
      if (exists('good.coeffs.df')){
        suppressMessages(good.coeffs.df <- join(good.coeffs.df, df, type = 'full'))
      } else{
        good.coeffs.df = df
      }
    }
    
    # Save best lasso to use for later prediction outside function
    if (metrics$Accuracy > best.acc){
      best.acc = metrics$Accuracy
      best.lasso = lasso.fit
    }
    
    if (metrics$Accuracy == 1){
      # If we're interested in finding perfect models, add 1 to the iterator
      if (!is.null(n.perfect)){
        i = i+1
      }
      
      # If the model had perfect prediction, track the features used
      n.perfect.models = n.perfect.models + 1
      perfect.features = c(perfect.features, as.character(df$Original_ID))
      
      # Print progress
      if (n.perfect.models %% 10 == 0){ cat('Perfect models found:', n.perfect.models, '\n')}
      
      # Also track the feature coefficients used
      if (exists('perfect.coeffs.df')){
        suppressMessages(perfect.coeffs.df <- join(perfect.coeffs.df, df, type = 'full'))
      } else{
        perfect.coeffs.df = df
      }
      
      # And the samples in the training and testing groups
      perfect.train.samples = c(perfect.train.samples, row.names(feature.matrix))
      perfect.test.samples = c(perfect.test.samples, row.names(testing.matrix))
      
    } 
    
    if (metrics$Accuracy >= good.acc){
      # Track features used in 90% accurate models
      good.features = c(good.features, as.character(df$Original_ID))
    }
    
    # If we're just tracking total number of models, add one to iterator
    if (!is.null(n.rounds)){
      i = i + 1
      # Print progress
      if (i %% 10 == 0) { cat('Rounds completed:', i, '\n') }
    }
  } 
  
  if (!exists('perfect.coeffs.df')){
    perfect.coeffs.df = list()
  }
  if (!exists('good.coeffs.df')){
    good.coeffs.df = list()
  }
  
  # Save and return 
  final.results = list(ModelNote = model.note,
                       TrainPercent = train.percent,
                       PerfectCoeffsDF = perfect.coeffs.df,
                       AllCoeffsDF = all.coeffs.df,
                       GoodCoeffsDF = good.coeffs.df,
                       Features = perfect.features,
                       NumModels = n.models,
                       NumPerfect = n.perfect.models,
                       Accuracy = all.acc,
                       AllPerfMetrics = all.perf.metrics, 
                       GoodAcc = good.acc,
                       PerfectFeatures = perfect.features,
                       GoodFeatures = good.features,
                       PerfectTrainSamples = perfect.train.samples,
                       PerfectTestSamples = perfect.test.samples,
                       AllTrainSamples = all.train.samples,
                       AllTestSamples = all.test.samples,
                       LassoFit = best.lasso,
                       Preds = lasso.preds,
                       ModSensitivities = mod.sensitivities,
                       ModSpecificities = mod.specificities,
                       AUC = roc.auc)
  return(final.results)
}

CoeffsDFBoxplot = function(coeffs.df, plot.group = 'both', plot.subset = NA, 
                           plot.names = T, main = '', ylab = '', col1 = NULL,
                           col2 = NULL){
  
  # Group to plot can be 'all', 'positive' or 'negative'
  # Subset can be NA for no subset, or a number for showing only the top n
  
  # Colors
  if (is.null(col1)){
    col1 = cols[1]
    col2 = cols[10]
  }
  
  # Format data
  model.coeffs = coeffs.df
  row.names(model.coeffs) = model.coeffs[,1]
  model.coeffs[,1] = NULL
  model.coeffs = t(model.coeffs)
  
  # Order by coefficient value
  prot.avg = apply(model.coeffs, 2, mean, na.rm = T)
  model.coeffs = model.coeffs[,order(prot.avg)]
  prot.avg = prot.avg[order(prot.avg)]
  
  # Subset
  if (is.na(plot.subset)){
    if (plot.group == 'positive'){
      sub = length(prot.avg[prot.avg > 0])
    }
    if (plot.group == 'negative'){
      sub = length(prot.avg[prot.avg < 0])
    }
    if (plot.group == 'both'){
      sub = length(prot.avg)
    }
  } else{
    sub = plot.subset
  }
  
  if (plot.group == 'positive'){
    
    if (sub > nrow(model.coeffs[, which(prot.avg > 0)])) {
      cat('Plot subset is larger than the number of positive coefficients. Plotting all positive coefficients instead.')
      sub = nrow(model.coeffs[, which(prot.avg > 0)])
    }
    
    # If there is only one positive feature, plot it
    #if (is.null(dim(model.coeffs[,which(prot.avg > 0)]))){
    #  coeffs.to.plot = model.coeffs[, which(prot.avg > 0)][1:sub]
    # Otherwise plot multiple positive features
    #} else{
    coeffs.to.plot = model.coeffs[, which(prot.avg > 0)]
    coeffs.to.plot = coeffs.to.plot[,(ncol(coeffs.to.plot)-(sub-1)):ncol(coeffs.to.plot)]
    #}
  } 
  if (plot.group == 'negative'){
    if (sub > nrow(model.coeffs[, which(prot.avg < 0)])) {
      cat('Plot subset is larger than the number of negative coefficients. Plotting all negative coefficients instead.')
      sub = nrow(model.coeffs[, which(prot.avg < 0)])
    }
    coeffs.to.plot = model.coeffs[, which(prot.avg < 0)][,1:sub]
  } 
  if (plot.group == 'both'){
    if (sub == ncol(model.coeffs)){
      coeffs.to.plot = model.coeffs
    } else{
      subsub = floor(sub/2)
      coeffs.to.plot = model.coeffs[,c(1:subsub, (ncol(model.coeffs)-(subsub-1)):ncol(model.coeffs))]
    }
  }
  prot.avg = prot.avg[colnames(coeffs.to.plot)]
  
  #if (ncol(coeffs.to.plot) <= 30){
  #  names.to.plot = colnames(coeffs.to.plot)
  #} else{
  #  names.to.plot = rep('', ncol(coeffs.to.plot))
  #}
  if (plot.names){
    names.to.plot = colnames(coeffs.to.plot)
  } else{
    names.to.plot = rep('', ncol(coeffs.to.plot))
  }
  
  par(las = 1)
  boxplot(coeffs.to.plot,
          horizontal = T,
          main = main,
          ylab = ylab,
          xlab = 'Model coefficient',
          names = names.to.plot,
          col = ifelse(as.logical(prot.avg > 0), col1, col2))
  abline(v = 0, col = 'grey')
}

AccuracyCDF = function(acc, plot.cols = NULL, plot.verticals = F){
  
  # acc should be a list of named accuracy vectors
  
  if (is.null(plot.cols)){
    plot.cols = sample(brewer.pal(11, 'Spectral'), length(acc))
  }
  
  plot(ecdf(acc[[1]]),
       verticals = plot.verticals,
       col = plot.cols[1],
       main = 'Accuracy CDF',
       xlab = 'Accuracy',
       xlim = c(0,1))
  
  if (length(acc) > 1){
    for (i in 2:length(acc)){
      plot(ecdf(acc[[i]]),
           verticals = plot.verticals,
           col = plot.cols[i],
           add = T)
    }
  }
  legend(x = 0, y = 1, names(acc), col = plot.cols, pch = 19, bty = 'n')
}

ReportPerfMetrics = function(predicted.labels, true.labels, pos.class){
  # Calculate the accuracy, precision and recall for two-class prediction
  tp = sum(true.labels == pos.class & predicted.labels == pos.class)
  fp = sum(true.labels != pos.class & predicted.labels == pos.class)
  tn = sum(true.labels != pos.class & predicted.labels != pos.class)
  fn = sum(true.labels == pos.class & predicted.labels != pos.class)
  n = tp + fp + tn + fn
  
  accuracy = (tp + tn)/n
  precision = tp/(tp + fp)
  recall = tp/(tp + fn)
  
  return(list(TP = tp,
              TN = tn,
              FP = fp,
              FN = fn,
              Accuracy = accuracy, 
              Precision = precision, 
              Recall = recall))
}