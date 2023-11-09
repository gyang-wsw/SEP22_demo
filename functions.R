to_base_names <- function(vect) {
  str_replace(vect, '_(Clicks|Impressions|GRPs)(_scaled)?(_0.[0-9])?$', '')
}

dim_ret_function <- function(variable, mu) {
  atan(variable/mu)
}

exp_prod <- function(mat) {
  apply(exp(mat), 1, prod)
}

create_coef_mult <- function(beta, data) {
  t(apply(data, 1, function(x) beta * x))
}

create_media_effect <- function(with_media, without_media) {
  with_media - without_media
}

create_effects_individual <- function(model_y, coef_mult) {
  model_y * (1-exp(-coef_mult))
}

adstock <- function(x, th) {
  adstocked <- vector(length=length(x))
  adstocked[1] <- x[1]
  for (t in 2:length(x)) {
    adstocked[t] <- x[t] + th*adstocked[t-1]
  }
  adstocked
}

adstock_group <- function(df, c, th)
  #c is name of column
  suppressMessages({
    df %>% group_by(DMA) %>% select(all_of(c)) %>% mutate(across(everything(), function(x) adstock(x, th))) %>% pull(c)
  })

#steeper with lower a
hill <- function(x, g, a) {
  1 / (1 + (x/g)^-a)
}

calc_spends_effects <- function(spends_sum, effects_sum, normalize=TRUE) {
  #hacky sol for TV
  tv_effect <- effects_sum[grep('^TV', names(effects_sum))] %>% sum()
  effects_sum <- effects_sum[grep('^TV', names(effects_sum), invert=TRUE)]
  effects_sum['TV'] <- tv_effect
  
  spends_effects <- merge(data.frame(effects_sum), spends_sum, by='row.names', all.x=TRUE)
  
  spends_effects$effects_sum <- spends_effects$effects_sum / sum(spends_effects$effects_sum)
  spends_effects$spends <- spends_effects$spends / sum(spends_effects$spends)
  
  spends_effects
}

#tools for param search

#apply to to_combine (pre)
evaluate_tstat <- function(var, adstock_ths, model) {
  #assumes var is combined and adjusts those params
  combined_units <- model1$addUnits(var, model1$model_dict$combine[[var]])
  th_named <- setNames(adstock_ths, combined_units$to_combine)
  model$adstock(th_named)
  model$combine()
  model$runOLS(prep=FALSE)
  summary(model$reg)$coefficients[combined_units$name_unit, c('Estimate', 'Pr(>|t|)')]
}

max_tstat <- function(var, model, iters=100L) {
  ng <- import('nevergrad', delay_load=TRUE)
  
  combined_units <- model1$addUnits(var, model1$model_dict$combine[[var]])
  p <- length(combined_units$to_combine)
  
  instrumentation <- ng$p$Array(shape = tuple(p), lower=0, upper=0.8)
  optim <- ng$optimizers$registry[['CMA']](parametrization=instrumentation, budget=iters)
  
  for (i in 1:optim$budget) {
    th <- optim$ask()
    
    tstat <- evaluate_tstat(var, th$value, model1)['Estimate']
    
    optim$tell(th, -1 * tstat)
  }
  
  t_max <- optim$provide_recommendation()$value
  evaluate_tstat(var, t_max, model1) %>% print()
  
  list(optim=optim, params=setNames(t_max, combined_units$to_combine))
}

#apply to singles
#warning: make sure after search for first var, model and data is same as before search
adstock_range <- seq(0,0.8,0.1)
params_search <- function(vars, model) {
  transforms_results <- list()
  for (x in vars) {
    xt <- matrix(NA, nrow=length(adstock_range), ncol=2, dimnames=list(adstock_range, c('Estimate', 'Pr(>|t|)')))
    
    v <- model$df_transformed[[x]]
    #adstock_current <- model$adstocks[[x]]
    
    for (i in 1:length(adstock_range)) {
      th <- adstock_range[i]
      th_named <- setNames(c(th), x)
      model$adstock(th_named)
      model$runOLS(prep=FALSE)
      xt[i,] <- summary(model$reg)$coefficients[x, c('Estimate', 'Pr(>|t|)')]
    }
    
    #reset adstocks to defined adstocks, if not then 0
    model$df_transformed[[x]] <- v
    #model$adstocks[[x]] <- adstock_current
    #all.equal(model1$df, model1$df_transformed)
    
    transforms_results[[x]] <- xt
  }
  transforms_results
}

dimret_range <- c(seq(0.1, 0.9, 0.1), seq(1, 5, 1))
params_search_dimret_adstock <- function(x, model) {
  results <- list()
  n1 <- length(adstock_range)
  n2 <- length(dimret_range)
  results$grid_est <- array(NA, dim=c(n1,n2), dimnames=list(adstock_range, dimret_range))
  results$grid_pval <- array(NA, dim=c(n1,n2), dimnames=list(adstock_range, dimret_range))
  
  v <- model$df_transformed[[x]]
  
  for (i in 1:n1) {
    for (j in 1:n2) {
      th <- adstock_range[i]
      steep <- dimret_range[j]
      model$adstock(setNames(c(th), x))
      model$dimret(setNames(c(steep), x))
      
      model$runOLS(prep=FALSE)
      results$grid_est[i,j] <- summary(model$reg)$coefficients[x, 'Estimate']
      results$grid_pval[i,j] <- summary(model$reg)$coefficients[x, 'Pr(>|t|)']
      
      model$df_transformed[[x]] <- v #dimret can only be run once, this is to reset
    }
  }
  
  results
}

params_search_2d <- function(a, b, model) {
  results <- list()
  n <- length(adstock_range)
  results$grid_est <- array(NA, dim=c(n,n,2), dimnames=list(adstock_range, adstock_range, NULL))
  results$grid_pval <- array(NA, dim=c(n,n,2), dimnames=list(adstock_range, adstock_range, NULL))
  
  v1 <- model$df_transformed[[a]]
  v2 <- model$df_transformed[[b]]
  
  for (i in 1:n) {
    for (j in 1:n) {
      th_a <- adstock_range[i]
      th_b <- adstock_range[j]
      th_named <- setNames(c(th_a, th_b), c(a,b))
      model$adstock(th_named)
      model$runOLS(prep=FALSE)
      results$grid_est[i,j,] <- summary(model$reg)$coefficients[c(a,b), 'Estimate']
      results$grid_pval[i,j,] <- summary(model$reg)$coefficients[c(a,b), 'Pr(>|t|)']
    }
  }
  
  model$df_transformed[[a]] <- v1
  model$df_transformed[[b]] <- v2
  
  results
}