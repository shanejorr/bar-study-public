# custom functions that are used in doctrinal-required
req_class_cum_rank <- function(df, req_class_remove) {
  
  # this function calculates the cumulative GPA of all required classes except the one
  # referenced in the parameter 'req_class_remove' and LMC and BEF
  
  # parameters:
  #  df: data frame containing class grade information
  #  req_class_remove: the required class that will be removed prior to 
  #         calculating total GPA for all required classes except the one entered ehre 
  
  keep_cols <- colnames(df)[str_detect(colnames(df),  "cum_gpa_.*|credits_.*")]
  
  calc_total_gpa <- df %>%
    select(keep_cols) %>%
    # remove standardized grade columns
    select(-contains('_ctr')) %>%
    # remove the specific class that we are comparing total GPA against
    select(-contains(req_class_remove), -contains("lmc"))
  
  # divide total number of column by two, so we can separate grade and class rank columns
  total_cols <- ncol(calc_total_gpa)
  divider <- total_cols / 2
  start_credits <- divider+1
  
  # multiply class grade by number of hours for all classes
  cum_gpa <- calc_total_gpa[1:divider] * calc_total_gpa[start_credits:total_cols]
  
  # sum rows for total cumulative GPA in required courses
  df['grade_sum'] <- rowSums(cum_gpa, na.rm=T)
  
  # calculate percentile rank by class
  df <- df %>%
    group_by(class) %>%
    # calculate standardized GPA
    mutate(grade_sum = scale(grade_sum)) %>%
    ungroup()
  
  return(df$grade_sum)
}

mod_with_total_grade <- function(df, req_class_mod) {
  # creates linear regression model with total cumulative grades minus one class as predictors
  
  df_mod <- df
  
  # create GPA column name for class
  req_class_name <- glue("cum_gpa_{req_class_mod}_ctr")
  
  # model formula
  req_class_form <- glue("result ~ cum_gpa_doc + {req_class_name} + state_rate_std + (1 | class)")
  
  # create dataset that has the total class rank for all required classes except included class
  df_mod$cum_gpa_doc <- req_class_cum_rank(df, req_class_mod)
  
  # prior probabilities for the regression coefficients in models with cumulative gpa
  coef_prior_req <- normal(location = c(1, 0, .3),
                           scale = c(2.5, 1, 2.5))
  
  # create model
  req_mod <- single_class_model(df_mod, req_class_form, coef_prior_req, weight = F)
  
  return(req_mod)
  
}

# get order for plot axis by ranking HDI of coefficients
create_axis_orders <- function(posterior_draws) {
  
  posterior_draws %>%
    group_by(class_name) %>%
    median_hdi(estimate) %>%
    arrange(estimate) %>%
    select(class_name) %>%
    .[[1]]
}

create_cum_gpa_draws <- function(mod, class_col) {
  # function that produces posterior for overall cum gpa and individual gpa coefficients
  
  # required so we can use object name in function
  class_name <- as.name(class_col)
  
  mod %>%
    spread_draws(`cum_gpa_.*`, regex = T) %>%
    select(contains("cum_gpa")) %>%
    rename(estimate = class_col) %>%
    mutate(class_name = !!class_col)
}

# create posterior predictions for fitted draws of individual 1L classes
make_fitted_draws <- function(df, class_col, fit_mod) {
  
  # standardized class GPA values to use for prediction
  # first, set length of range of values, for use in creating prediction dataset
  # msut be even number since we will divide it by two
  n <- 20
  class_pred <- round(seq(-.8, .8, length.out = n), 3)
  
  # create predictions for students at the bottom and top 33% of cum gpa distribution
  # value represent standard deviations for these percentiles
  cum_gpa_values <- round(qnorm(c(.3333, .6666), mean = 0, sd = 1), 3)
  
  tibble(
    jur = rep('NC_ube', n),
    class = rep(2019, n),
    cum_gpa_doc = rep(cum_gpa_values, n/2),
    !!class_col := class_pred,
    state_rate_std = rep(.2, n)
  )  %>%
    add_fitted_draws(fit_mod, n = 300, seed = 456)%>%
    mutate(class_name = str_remove(!!class_col, " ")) %>%
    ungroup() %>%
    rename(class_grade = !!class_col) %>%
    select(class_name, class_grade, cum_gpa_doc, .row, .draw, .value) %>%
    mutate(cum_gpa_doc = as.character(cum_gpa_doc),
           cum_gpa_doc = str_replace(cum_gpa_doc, as.character(cum_gpa_values[1]), "Bottom 33% GPA"),
           cum_gpa_doc = str_replace(cum_gpa_doc, as.character(cum_gpa_values[2]), "Top 33% GPA"))

}
