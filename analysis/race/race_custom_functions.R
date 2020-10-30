#####################################################################
#
# This script contains custom functions that are used in race.Rmd
#
#####################################################################

# function to produce line plots for racial differences in the raw data
plot_race_diff <- function(df, x, y, c, title_label, x_label, y_label, add_perc) {
  
  x <- enquo(x)
  y <- enquo(y)
  c <- enquo(c)
  
  plot_race <- ggplot(df, aes(!!x, !!y, color = !!c)) +
    geom_point() +
    geom_line(aes(group = !!c)) +
    scale_colour_brewer(type = "qual", palette = "Set2", aesthetics = c("colour", "fill")) +
    theme_minimal() +
    labs(title = title_label,
         x = x_label,
         y = y_label,
         color = NULL,
         linetype = NULL) +
    theme(plot.title = element_text(size = rel(1.1)), legend.position = "bottom")
  
  if (add_perc) {
    plot_race <- plot_race +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))
  }
  
  return(plot_race)
}

# function to calculate posterior, which will be used for this and the GPA / race interaction model
create_interaction_posterior <- function(model) {
  
  # extract the posterior of interaction terms from a model that includes race and year interacting
  # parameters:
  #    model: the model with interactions terms
  
  year_minority_post <- model %>%
    # extract all coefficients that pertain to White status
    spread_draws(`.*White|.*White.*`, regex = T) %>%
    select(-.chain, -.iteration, -contains('sigma')) %>%
    # add fixed Black coefficient to all Black / class interaction terms
    mutate_at(vars(contains("classificationWhite:year")), ~(classificationWhite + .))
  
  # make column names only the years (example: 2016-2017)
  year_col_names <- str_extract(colnames(year_minority_post), "[0-9]{4}-[0-9]{4}")
  year_col_names <- year_col_names[!is.na(year_col_names)]
  colnames(year_minority_post) <- c(".draw", "2017D-2019", year_col_names)
  
  # currently, each column is a year / minority status coefficient
  # convert to long form where each row is a trace value and year / minority status
  year_minority_post <- year_minority_post %>%
    pivot_longer(cols = matches("^[0-9]"), names_to = 'class_name', values_to = 'estimate')
}

avg_plot_title <- function(avg_metric) {
  # create plot titles for plots that show average difference between races
  glue("Average {avg_metric} by Year and Race")
}

metric_averages <- function(df, metric_cols) {
  # create a table of average law school GPA, LSAT, and undergrad GPA, by race and graduating class
  # all values are scaled
  
  # parameters:
  #   df: dataframe with class and classification columns, and each metric in a different columns
  #   metric_cols: vector of column names for the metrics to calculate average of
  
  metric_name_recode <- c(cum_gpa_ctr = "LGPA", 
                          admission_idx = "Admissions Index",
                          lsat = "LSAT", 
                          u_gpa = "UGPA")
  
  df %>%
    group_by(class, classification) %>%
    summarize_at(vars(!!metric_cols), ~weighted.mean(x = ., w = weights, na.rm=T)) %>%
    pivot_longer(cols = !!metric_cols, names_to = 'metric', values_to = 'estimate') %>%
    mutate(metric = factor(metric, levels = names(metric_name_recode)),
           metric = recode(metric, !!!metric_name_recode))
}

calc_bar_pass <- function(weighted_df, matching_group) {
  
  # this function calcualtes the pass passage rates from an matched and weighted dataframe
  # the dataframe is created with `metric_averages`
  
  # parameters:
  #   weighted_df = weighted dataframe created with `metric_averages`
  #   matching_group = string with a description of the matched dataset
  
  weighted_df %>%
    group_by(classification, result) %>%
    summarize(sum_weight = sum(weights), # creates total counts based on weights
              count = n()) %>%
    ungroup() %>%
    mutate(result = recode(result, `0` = "Fail", `1` = "Pass")) %>%
    pivot_wider(id_cols = "classification", names_from = "result", values_from = "sum_weight") %>%
    mutate(pass_rate = Pass / (Pass + Fail),
           matching_group = !!matching_group) 
  
}

# the function creates class rank quantile
create_quantile <- function(perc_rank_col, num_quantiles) {
  
  num_quant <- num_quantiles
  quantile_labels <- seq(1, num_quant)
  
  quant <- cut_number(perc_rank_col, n = num_quant, labels = quantile_labels)
  
  return(quant)
  
}
