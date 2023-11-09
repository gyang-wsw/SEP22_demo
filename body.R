library(R6)
library(jsonlite)
library(rstan)

metric_cols <- c('GRPs', 'Clicks', 'Impressions')

#SEP22
#model_params <- list(
#  search=list(dep_var='Search_GM_Clicks', ethnicity=c('GM'), exclude_media='Search_GM_Clicks'),
#  others_under150=list(dep_var='enrolls_others_under150', ethnicity=c('GM')),
#  AA_under150=list(dep_var='enrolls_AA_under150', ethnicity=c('GM', 'AA'))
#)

#SEP23
model_params <- list(
  search=list(dep_var='Search_GM_Clicks', ethnicity=c('GM')),
  others=list(dep_var='enrolls_others', ethnicity=c('GM')),
  AA=list(dep_var='enrolls_AA', ethnicity=c('GM', 'AA')),
  latino=list(dep_var='enrolls_latino', ethnicity=c('GM', 'Latino'))
)

modelClass <- R6Class(
  "modelClass", 
  list(
    model_args = NULL,
    df = NULL,
    df_transformed = NULL,
    model_dict = NULL,
    combine_map = NULL,
    predictors = NULL,
    demo_cols = NULL,
    media_vars = NULL,
    adstocks_init = NULL,
    adstocks = c(),
    predictors_media = NULL,
    predictors_media_vars = NULL,
    reg = NULL,
    model = list(),
    coef_mult = NULL,
    components = NULL,
    effects = NULL,
    effects_sum = NULL,
    preds = NULL,
    initialize = function(data, model_args) {
      self$model_args <- model_args
      self$df <- data$df
      
      #use transformed df from now on
      self$df_transformed <- self$df
      
      self$demo_cols <- data$demo_cols
      #base_ethnicity <- media_names_df %>% filter(Ethnicity==self$model_args$ethnicity) %>% pull(varname)
      self$media_vars <- data$media_vars #[names(data$media_vars) %in% base_ethnicity]
      
      self$adstocks_init <- setNames(data$adstocks_init[[self$model_args$dep_var]], data$adstocks_init$var)
      
      path <- paste(model_dir, self$model_args$model, self$model_args$model_file, sep='/')
      self$model_dict <- read_json(path)
    }
  )
)

modelClass$set("public", "addUnits", function(name, values) {
  for (unit in metric_cols) {
    try({
      # Attempt to combine the specified columns
      to_combine <- paste0(values, "_", unit)
      combined <- self$df_transformed[to_combine] %>% rowSums()
      name_unit <- paste0(name, "_", unit)
      
      # Break out of the loop if successful
      break
    }, silent=TRUE)
  }
  lst(name_unit, to_combine)
})

#assume combine_dict is valid i.e. doesn't contain vars not in ethnicity
modelClass$set("public", "combine", function() {
  combine_dict <- self$model_dict$combine
  
  self$combine_map <- list()
  combined <- list()
  for (k in names(combine_dict)) {
    
    addunits <- self$addUnits(k, combine_dict[[k]])
    to_combine_unit_k <- addunits$to_combine
    
    #make combine_map
    n <- length(combine_dict[[k]])
    combine_map_k <- rep(k, n) %>% as.list()
    names(combine_map_k) <- to_combine_unit_k
    self$combine_map <- c(self$combine_map, combine_map_k)
    
    self$df_transformed[addunits$name_unit] <- self$df_transformed[to_combine_unit_k] %>% rowSums()
    
    #in case we want to transform (adstock) after combining
    self$df[addunits$name_unit] <- self$df[to_combine_unit_k] %>% rowSums()
    
    combined <- append(combined, addunits$name_unit)
  }
  
  to_combine <- names(self$combine_map) %>% to_base_names()
  singles <- setdiff(self$media_vars %>% to_base_names(), to_combine)
  self$predictors_media <- c(self$media_vars[singles], combined)
  self$predictors_media <- self$predictors_media[!(self$predictors_media %in% self$model_dict$exclude)]
  self$predictors_media <- self$predictors_media %>% str_subset(paste(self$model_args$ethnicity, collapse='|'))
  
  self$predictors <- c(self$demo_cols, self$predictors_media, self$model_args$add_cols)
  
  invisible(self)
})

modelClass$set("public", "adstock", function(adstock_ths) {
  #adstock_ths is a named vector
  #it gets replaced so can call multiple times without stacking
  
  for (n in names(adstock_ths)) {
    if (adstock_ths[[n]] == 0) {
      self$df_transformed[[n]] <- self$df[[n]]
    }
    else {
      self$df_transformed[[n]] <- adstock_group(self$df, n, adstock_ths[n])
    }
  }
  
  #update current adstock
  #doesn't include adstocks_post
  self$adstocks[names(adstock_ths)] <- adstock_ths
  
  invisible(self)
})

modelClass$set("public", "dimret", function(steep) {
  #can only call once
  
  for (n in names(steep)) {
    self$df_transformed[[n]] <- hill(self$df_transformed[[n]], 0.5, steep)
  }
  
  invisible(self)
})

modelClass$set("public", "transform", function(pre=TRUE) {
  if (pre) {
    #initialize adstocks
    #don't pre-adstock to_combine that will be adstocked after combining
    to_combine_post <- self$model_dict$combine[self$model_dict$adstocks_post %>% names() %>% to_base_names()] %>% unlist()
    adstocks_init <- self$adstocks_init[!(names(self$adstocks_init) %in% to_combine_post)]
    self$adstock(adstocks_init)
    
    #override with defined adstocks
    self$adstock(self$model_dict$adstocks_pre %>% unlist())
  
    #dimret transform is after hill, make sure to only transform once
    self$dimret(self$model_dict$dimret_pre %>% unlist())
  }
  else {
    self$adstock(self$model_dict$adstocks_post %>% unlist())
    self$dimret(self$model_dict$dimret_post %>% unlist())
  }
  
  invisible(self)
})

modelClass$set("public", "runOLS", function(prep=TRUE, ...) {
  #do pre-reg steps outside of function for more flexibility
  #transform (has units) before combine (need to add units)
  #as a result units chosen by combine may be different from those defined by transform
  
  if (prep) {
    model1$transform(pre=TRUE) #for to_combine
    model1$combine() #transform before and after combine
    model1$transform(pre=FALSE) #for combined
  }
  
  formula_predictors <- paste(self$predictors, collapse="+")
  
  #other transforms besides for media
  for (n in names(self$model_dict$transforms)) {
    formula_predictors <- gsub(n, self$model_dict$transforms[[n]], formula_predictors)
  }
  
  if (self$model_args$log) {
    formula <- sprintf('log(%s+1) ~ %s', self$model_args$dep_var, formula_predictors)
  }
  else {
    formula <- sprintf('%s / %s ~ %s', self$model_args$dep_var, as.character(self$model_args$scale), formula_predictors)
  }
  
  self$reg <- lm(as.formula(formula), data=self$df_transformed, x=TRUE, y=TRUE)
  self$model$X <- self$reg$x
  self$model$y <- self$reg$y
  self$model$coefficients <- self$reg$coefficients
  
  invisible(self)
})

modelClass$set("public", "runRidge", function(...) {
  #use matrix from lm formula for dummies
  formula_predictors <- paste(self$predictors, collapse="+")
  formula <- sprintf('%s / %s ~ %s', self$model_args$dep_var, as.character(self$model_args$scale), formula_predictors)
  
  reg_lm <- lm(as.formula(formula), data=self$df_transformed, x=TRUE, y=TRUE)
  
  X <- reg_lm$x[,-1]
  y <- reg_lm$y
  
  #param constraints
  xl <- rep(-Inf, ncol(X))
  names(xl) <- colnames(X)
  xl[self$predictors_media %>% unlist()] <- 0
  
  #intercept=TRUE gives diff answer than intercept=FALSE and column of 1s
  self$reg <- glmnet::glmnet(X, y, 'gaussian', intercept=TRUE, standardize=FALSE, lower.limits=xl, alpha=0, ...)
  self$model$X <- cbind(`(Intercept)`=1, X)
  self$model$y <- y
  self$model$coefficients <- coef(self$reg)[,1]
  
  invisible(self)
})

modelClass$set("public", "prepareBayesian", function(psbeta2, psbeta3) {
  #use matrix from lm formula for dummies
  formula_predictors <- paste(self$predictors, collapse="+")
  formula <- sprintf('%s / %s ~ %s', self$model_args$dep_var, as.character(self$model_args$scale), formula_predictors)
  
  reg_lm <- lm(as.formula(formula), data=self$df_transformed, x=TRUE, y=TRUE)
  
  X <- reg_lm$x[,-1] #remove intercept, contained in stan
  y <- reg_lm$y
  
  exclude_scale <- grep(paste(
    c(self$df %>% select(where(is.character)) %>% colnames(), 'Week_num'),
    collapse='|'), colnames(X))
  
  #scale to mean=1, need to exclude dummy variables
  #X[,exclude_scale] <- apply(X[,exclude_scale], 2, function(x) x/mean(x))
  #y <- y/mean(y)
  
  self$model$X <- X
  self$model$y <- y
  
  predictors_media <- self$predictors_media %>% unlist()
  predictors_base <- setdiff(colnames(X), predictors_media)
  
  data_stan <-list(
    N = nrow(X),
    K1 = length(predictors_base),
    K2 = length(model1_coefs_found_zero),
    K3 = length(model1_coefs_found_weak),
    y = y,
    X1 = X[, predictors_base],
    X2 = X[, names(model1_coefs_found_zero)],
    X3 = X[, names(model1_coefs_found_weak)],
    psbeta2 = rep(psbeta2, length(model1_coefs_found_zero)),
    pmubeta3 = model1_coefs_found_weak,
    psbeta3 = psbeta3
  )
  
  data_stan
})

modelClass$set("public", "runBayesian", function(stan_file, data_stan) {
  self$reg <- stan(file = stan_file, data = data_stan, verbose=TRUE, open_progress=FALSE)
  
  invisible(self)
})

modelClass$set("public", "runModel", function(model='ols', ...) {
  if (model == 'ols') {
    self$runOLS(...)
  }
  
  #note: both glmnet and stan doesn't support categorical so need to create dummies
  #also need to add log support
  
  else if (model == 'ridge') {
    self$runRidge(...)
  }
  
  else if (model == 'bayesian') {
    self$runBayesian(...)
  }
  
  invisible(self)
})

modelClass$set("public", "evaluate", function(adstock_ths) {
  #returns r2- a scalar
  
  self$adstock(adstock_ths)
  self$runRidge(lambda=0)
  
  self$reg$dev.ratio
})

modelClass$set("public", "breakout_combined", function(coefs_no_units, no_breakout) {
  #breakouts combined by mult its coef with to_combine, then cbinds singles
  
  combine_map <- self$combine_map[!(self$combine_map %in% no_breakout)]
  
  combined_coefs <- data.frame(combined=do.call(rbind, combine_map))
  combined_coefs$name <- rownames(combined_coefs)
  combined_coefs <- merge(combined_coefs, data.frame(coefs=coefs_no_units), by.x='combined', by.y='row.names')
  
  if (nrow(combined_coefs) != 0) {
    coef_mult_combined <- create_coef_mult(combined_coefs$coefs, self$df_transformed[, combined_coefs$name])
    colnames(coef_mult_combined) <- to_base_names(colnames(coef_mult_combined))
  } else {
    coef_mult_combined <- numeric(0)
  }
  
  remove <- combined_coefs$combined
  
  self$predictors_media_vars <- to_base_names(self$predictors_media)
  self$predictors_media_vars <- c(self$predictors_media_vars[!(self$predictors_media_vars %in% remove)], colnames(coef_mult_combined))
  
  self$coef_mult <- cbind(self$coef_mult[, !(colnames(self$coef_mult) %in% remove)], coef_mult_combined)
  
  invisible(self)
})

modelClass$set("public", "create_results_model", function(no_breakout) {
  coefs_no_units <- setNames(self$model$coefficients, to_base_names(names(self$model$coefficients)))
  
  self$coef_mult <- create_coef_mult(coefs_no_units, self$model$X)
  
  self$breakout_combined(coefs_no_units, no_breakout)
  
  invisible(self)
})

modelClass$set("public", "create_results", function(groups=list(), breakout=list(), no_breakout=c()) {
  #groups: named list of vectors
  
  individual <- length(groups) == 0
  
  coef_mult_breakouts <- list()
  for (m in names(breakout)) {
    c <- self$model$coefficients[breakout[[m]]] * models[[m]]$coef_mult
    #assume same index as self$df
    c <- c[, models[[m]]$predictors_media_vars]
    
    coef_mult_breakouts[[m]] <- c
  }
  
  self$create_results_model(no_breakout)
  
  coef_mult <- self$coef_mult
  
  if (self$model_args$log) {
    inverseTransform <- exp(self$reg$fitted.values)-1
    aggFunction <- exp_prod
  } 
  else {
    inverseTransform <- self$df_transformed[[self$model_args$scale]] * self$reg$fitted.values
    aggFunction <- sum
  }
  
  preds <- data.frame(DMA=self$df_transformed$DMA, Week=self$df_transformed$Week_num,
                      y=self$df_transformed[[self$model_args$dep_var]], y_fitted = inverseTransform)
  preds$res <- preds$y - preds$y_fitted
  
  self$components <- data.frame(
    baseline_main=aggFunction(coef_mult[,setdiff(colnames(coef_mult), self$predictors_media_vars)])
  )
  
  if(length(breakout) != 0) {
    coef_mult_media <- do.call(
      data.frame, c(list(main=coef_mult[, self$predictors_media_vars]), coef_mult_breakouts)
    )
  } 
  else { #to see variables effects w/o search
    coef_mult_media <- data.frame(main=coef_mult[, self$predictors_media_vars])
  }
  
  #split into groups based on value
  regex <- sprintf("(main|%s)\\.", paste(names(breakout), collapse="|"))
  coef_mult_media <- sapply(split.default(coef_mult_media, sub(regex, "", colnames(coef_mult_media))), rowSums)
  
  #effects of media types instead of individual vars
  if (!individual) {
    coef_mult_media_groups <- list()
    
    for (e in names(groups)) {
      cols_e <- colnames(coef_mult_media) %in% groups[[e]]
      if (sum(cols_e) > 0) {
        coef_mult_media_groups[[e]] <- coef_mult_media[,cols_e,drop=FALSE] %>% rowSums()
      }
    }
    
    coef_mult_media <- data.frame(coef_mult_media_groups) %>% as.matrix()
  }
  
  with_media <- preds$y_fitted
  
  if (self$model_args$log) {
    without_media <- with_media*exp(-coef_mult_media)
  }
  else {
    without_media <- with_media - coef_mult_media
  }
  
  self$effects <- data.frame(
    enrolled = self$df_transformed[[self$model_args$dep_var]],
    media_effect = create_media_effect(with_media, self$components$baseline_main),
    
    #multiply each col by fitted
    create_media_effect(with_media, without_media)
  )
  
  effects_sum <- colSums(self$effects)
  self$effects_sum <- data.frame(name=names(effects_sum), effect=effects_sum)
  
  self$preds <- preds
  preds
})