################################################################
#
# This script contains functions used throughout the analysis files
#
################################################################

# models -----------------------------------------------------------------

# function to run models with the formual and the input
single_class_model <- function(df, model_formula, coef_prior, weight = F) {
  
  # for the intercept prior, assume 70% pass rate with all explanatory variables at 0,
  # which equates to a .85 log-odds
  intercept_prior <- student_t(df = 5, location = .85, scale = 3, autoscale = F)
  
  stan_glmer(formula(model_formula),
             data = df, adapt_delta = .99,
             prior_intercept = intercept_prior,
             prior = coef_prior,
             family = binomial("logit"),
             weights = if (weight == F) NULL else if (weight == T) weights,
             iter = 4000,
             seed = 12345)
}

# create a dataset with propensity score matches based on classes
matching_class_dataset <- function(df, class, match_formula, match_method) {
  
  # create matching object
  match_mod <- matchit(formula(match_formula),
                       data = df, method = match_method, replace = F)
  
  # create dataset of matched observations
  match_df <- match.data(match_mod)
  
  # add class name to end of dataset
  match_df$class_name <- class
  
  return(match_df)
  
}

# create a model based on a matched dataset of whether someone took an elective class
matching_class_model <- function(df, class, match_formula, match_method, model_formula, coef_prior) {

  # create matched dataset
  matched_data <- matching_class_dataset(df, class, match_formula, match_method)
  
  # create bayesian hierarchical model from matched dataset
  mod <- single_class_model(matched_data, model_formula, coef_prior, weight = T)
  
  return(mod)
  
}

# manipulate posterior draws ---------------------------------------------------

# function to extract HPI
hdi_single_class <- function(single_class_model, class) {
  
  # required so we can use object name in function
  class_name <- as.name(class)
  
  # remove '_rank' from end of class name
  class <- str_remove(class, "_rank")
  
  df <- single_class_model %>%
    spread_draws(!!class_name, perc_rank_center) %>%
    # calculate median highest posterior density interval
    median_hdi(!!class_name, perc_rank_center) %>%
    # add the class name, so that when we iterate through all models we can keep track of the class
    mutate(name = !!class) %>%
    # don't need these columns, because they are the same for every iteration
    select(-.width, -.point, -.interval)
  
  # make standard column names, so all class dataframes have the same column names
  colnames(df) <- c('class', 'class.lower', 'class.upper',
                    'perc_rank_center', 'perc_rank_center.lower', 'perc_rank_center.upper',
                    'class_name')
  
  return(df)
  
}

# create posterior prediction of elective classes
# used in doctrinal-electives.Rmd.
make_fitted_draws_electives <- function(class_name, fit_mod, perc_rank_col) {
  
  pred_num_perc_rank <- 50 
  
  prediction_grid <- data.frame(
    # create prediction dataset
    # perc_rank_center = rep(c(.25, .5, .75), each = 2),
    perc_rank_ctr = seq_range(perc_rank_col, n = pred_num_perc_rank),
    jur = "NC_ube",
    class = 2019,
    state_rate_std = .2,
    class_take = rep(c(0, 1), pred_num_perc_rank/2),
    stringsAsFactors = FALSE
    )

  # change column name to the name of the class, so that it matches the column name used in the model
  colnames(prediction_grid)[length(prediction_grid)] <- class_name
  
  prediction_grid <- prediction_grid %>%
    # create predictions from prediction dataset with model
    add_fitted_draws(fit_mod, n = 500) %>%
    ungroup() %>%
    select(perc_rank, took_class = class_name, .value, .row) %>%
    mutate(class_name = class_name)
  
  return(prediction_grid)
  
}

create_class_draws <- function(mod, class_col) {
  
  # create posterior draws of class coefficient from regression models
  
  # required so we can use object name in function
  class_name <- as.name(class_col)
  
  mod %>%
    spread_draws(!!class_name) %>%
    rename(estimate = class_col) %>%
    mutate(class_name = !!class_col) %>%
    select(estimate, class_name)
  
}

post_point_prob <- function(posterior_values, group_term, expression, pretty_output = F) {
  
  # calculates the percentage of the posterior distribution that is 
  # above or below a given threshold
  
  # parameters
  #   poseterior_values = dataframe of the posterior distribution
  #     can be created with `create_class_draws`
  #   group_term = column name that identifies groups
  #   expression = boolean expression to calculate with, should be wrapped with `expr`
  #     this will include the following: 
  #        column names that contains values, greater than or less than sign, and threshold value
  #     example: `expr(estimate > 1)`
  #   pretty_output = T or F, whether to print the probability as a percentage ("99%"): use TRUE
  #      or just decimals (.99): use FALSE
  
  # produce error message if pretty_output does not equal T or F
  if (!(pretty_output %in% c(T, F))) {
    stop("`pretty_output` must be either TRUE or FALSE")
  }
  
  group_term <- enquo(group_term)
  
  posterior_values <- posterior_values %>%
    group_by(!!group_term) %>%
    mutate(exceed_threshold = ifelse(!!expression, T, F)) %>%
    summarize(exceed_threshold = sum(exceed_threshold) / n(),
              .groups = 'drop') %>%
    mutate(exceed_threshold = ifelse(exceed_threshold > .99, .99, exceed_threshold),
           exceed_threshold = ifelse(exceed_threshold < .01, .01, exceed_threshold),
           exceed_threshold = round(exceed_threshold, 2))
  
  if (pretty_output) {
    posterior_values <- posterior_values %>%
      mutate(exceed_threshold = glue("{exceed_threshold*100}%"))
  }
  
  return(posterior_values)
  
}

interaction_trace <- function(mod, dummy_demo, group_demo) {
  
  # create dataframe of trace of interaction terms for a demgoraphic
  
  # pull out trace of model and only the coefficients for the interaction terms and main effect
  df <- mod %>%
    spread_draws(`^treatment_catTreatment.*`, regex=T) %>%
    select(starts_with("treatment"))
  
  # add the main effect to the interaction term for each interaction term
  for (i in seq_along(df)) {
    if (i == 1) next
    df[[i]] <- df[[1]] + df[[i]]
  }
  
  # convert from wide form to long form and relabel
  df <- df %>% 
    pivot_longer(cols = everything(),
                 names_to = 'demographic', values_to = "value") %>%
    mutate(demographic = str_remove(demographic, "treatment_catTreatment:"),
           demographic = str_replace_all(demographic, "treatment_catTreatment", !!dummy_demo),
           group = !!group_demo)
  
  return(df)
  
}

# create a dataset for making predictions with a Bayesian model, for whether person took class
class_prediction_dataset <- function(df, class_name, class_years, jur) {
  
  # parameters:
  #  class_name: the name of the class column in the original dataset
  #  class_years: years of graduating classes to include
  #  jur: bar jurisdictions to include
  
  df %>%
    data_grid(perc_rank_center = seq_range(perc_rank_center, 30, trim = .25, pretty = T),
              !!class_name := c(0, 1),
              class = !!class_years,
              jur = !!jur,
              state_rate_std = 0)
  
}

# function to find the probability that a coefficient (ex: year) is the highest or lowest
prob_rank <- function(post_dist, rank_num, pretty_print = F) {
  
  # paramerters:
  #   post_dist = posterior distribution of coefficients in long form
  #   rank_num = the rank order that we want to use to calculate the probability that coef. has that rank order
  #   pretty_print = converts percentages to printable format
  
  post_dist <- post_dist %>%
    group_by(.draw) %>%
    mutate(draw_rank = min_rank(estimate)) %>%
    # only keep the highest ranked year
    filter(draw_rank == !!rank_num) %>%
    group_by(class_name) %>%
    count() %>%
    ungroup() %>%
    mutate(perc_highest = n / sum(n))
  
  if (pretty_print) {
    post_dist <- post_dist %>%
      select(class_name, perc_highest) %>%
      mutate(perc_highest = round(perc_highest*100, 0),
             perc_highest = glue("{perc_highest}%"))
  }
  
  return(post_dist)
}

# visualizations -------------------------------------------------------------

# create point plot of regression coefficients
point_plot <- function(df, x, y, reorder_fctr, x_max, x_min, color_hex = "#9B2805", 
                       plot_title = NULL, x_title = NULL, y_title = NULL) {
  
  x <- enquo(x)
  y <- enquo(y)
  x_max <- enquo(x_max) 
  x_min <- enquo(x_min)
  reorder_fctr <- enquo(reorder_fctr)
  
  ggplot(df, aes(!!x, fct_reorder(!!y, !!reorder_fctr))) +
    geom_point(color = color_hex) +
    # bayesian 90% HPI
    geom_errorbarh(aes(xmax = !!x_max, xmin = !!x_min), color = color_hex) +
    labs(title = plot_title,
         x=x_title,
         y=y_title) +
    theme_minimal() + 
    theme(legend.position="none",
          plot.title = element_text(size = rel(1.1)))
  
}

plot_took_class <- function(df, x_value, y_value = class, color_value, x_label, x_limits = c(0, 1), plot_title,
                            y_sorting_order,
                            color_sorting_order = c("Took class", "Did not take class")) {
  
  # this function creates a point plot for classes on the y axis and
  # either bar passage or class rank on the x axis, grouped by students who
  # took or did not take course
  
  x_value <- enquo(x_value)
  y_value <- enquo(y_value)
  color_value <- enquo(color_value)
  
  ggplot(df, aes(!!x_value, factor(!!y_value, levels = y_sorting_order), 
                 color = factor(!!color_value, levels = color_sorting_order))) +
    geom_point() +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = x_limits) +
    scale_colour_brewer(palette = "Set2") +
    labs(title = plot_title,
         x = x_label,
         y = NULL,
         color = NULL) +
    theme_minimal() +
    theme(legend.position="bottom",
          plot.title = element_text(size = rel(1.1)))
  
}

# rank coefficients
rank_coef <- function(gathered_coef, ranks) {
  
  gathered_coef %>%
    group_by(.draw) %>%
    # rank coefficients within draws
    mutate(rank = min_rank(.value)) %>%
    filter(rank %in% !!ranks) %>%
    group_by(.variable) %>%
    count() %>%
    mutate(percentage = n / 8000) %>%
    arrange(desc(percentage)) %>%
    ungroup() %>%
    mutate(.variable = str_remove(.variable, "_rank"),
           .variable = recode(.variable, !!!map_names),
           percentage = round(percentage, 2) *100,
           percentage = glue("{percentage}%")) %>%
    select(.variable, percentage)
}

# plot of coefficients with point estimate, 50% CI, and 95% CI
class_coeff_plot <- function(class_draws, title_label, x_label, y_label, axis_order, x_limit = NULL, horz_line = 0, sub_label = NULL) {
  
  # create point plot with 95 and 50% intervals of coefficients
  
  coef_plot <- class_draws %>%
    group_by(class_name) %>%
    ggplot(aes(y = factor(class_name, levels = axis_order), x = estimate)) +
    stat_pointinterval(color = "#2960AC") +
    geom_vline(xintercept = horz_line, alpha = .6) +
    labs(title = title_label,
         subtitle = sub_label,
         x = x_label,
         y = y_label) +
    theme_minimal() +
    theme(plot.title = element_text(size = rel(1.1)))
  
  # add x limits if needed
  if (!is.null(x_limit)) {
    coef_plot <- coef_plot +
      lims(x = x_limit)
  }
  
  return(coef_plot)
  
}

class_bar_pass_prob <- function(df, lab_title, lab_subtitle) {
  
  # create line plot showing probability of bar passage based on whether
  # person took class, at all class ranks
  
  df %>%
    ggplot(aes(x, predicted, colour = group)) + 
    geom_line() +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill=group), alpha = .2, size = 0) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .20), labels = scales::percent) +
    scale_x_continuous(limits = c(.2, .8), breaks = seq(.2, .8, .2), labels = scales::percent) +
    scale_color_brewer(type = "Qual", palette = "Set2") +
    scale_fill_brewer(type = "Qual", palette = "Set2") +
    labs(title = lab_title,
         subtitle = lab_subtitle,
         x = "Cumulative Law School Class Percentile Rank",
         y = "Bar Passage Probability",
         color = NULL,
         fill = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(size = rel(1.1)))
}

# create a plot of predictions by law school class ranks from a model
# models are for whether student took class
plot_predictions <- function(data, x_var, class_name, x_axis_limit, plot_title, 
                             x_lab = "Law school class percentile rank", y_lab = "Bar passage probability", use_perc = T) {
  
  # parameters:
  #  data: dataframe from predictions (comes from add_fitted_draws)
  #  class_name: name of the class
  #  facet_plot: whether to facet the plot (probably by years)
  #  year_order: order of the years for the facet
  #  plot_title: title of the plot
  #  use_perc: convert the x and y axis to percentages
  
  class_name <- enquo(class_name)
  x_var <- enquo(x_var)
  
  pred_plot <- data %>%
    ggplot(aes(x = !!x_var, y = .value, color = !!class_name, fill = !!class_name)) +
    geom_line(alpha = .8, size = .8) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), 
                alpha = .1, linetype = 0, size = 0, show.legend = FALSE) +
    scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set2") +
    labs(title = plot_title,
         subtitle = "Bands are 90% credible intervals",
         x = x_lab,
         y = y_lab,
         color = NULL) +
    theme_minimal() +
    theme(plot.title = element_text(size = rel(1.1)),
          legend.position="bottom")
  
  if (use_perc == T) {
    pred_plot <- pred_plot +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .20), labels = scales::percent) +
      scale_x_continuous(limits = x_axis_limit, breaks = seq(x_axis_limit[1], x_axis_limit[2], .1), 
                         labels = scales::percent_format(accuracy = 1))
  }
  
  return(pred_plot)
  
}

# lollipop plot with a single group
single_group_point_plot <- function(df, x, y, color_var, plot_title, x_label, perc_label, x_offset = 1.07, sub_label = NULL) {
  
  x <- enquo(x)
  y <- enquo(y)
  color_var <- enquo(color_var)
  
  plt <- ggplot(df, aes(!!x, !!y, color = !!color_var)) + 
    geom_point() +
    geom_segment(aes(x=0, xend=!!x, y=!!y, yend=!!y)) +
    geom_text(aes(label = if (perc_label == T) {sprintf("%1.0f%%", !!x*100)} else round(!!x, 2), 
                  y = !!y, x = ifelse(!!x > 0, !!x + x_offset, !!x - x_offset)),
              color = "#333333", size = 3.3) +
    scale_color_brewer(type = "Qual", palette = "Set2") +
    labs(title = plot_title,
         subtitle = sub_label,
         x = x_label,
         y = NULL,
         color = NULL,
         fill = NULL) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.1)))
  
  if (perc_label == T) {
    plt <- plt +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))
  }
  
  return(plt)
}

# misc. analysis ----------------------------------------------------------------

# calculate average correlation for doctrinal courses
class_corr <- function(df, class_name, class_col, cor_col) {
  
  class_col <- enquo(class_col)
  cor_col <- enquo(cor_col)
  
  df %>%
    filter(str_detect(!!class_col, !!class_name)) %>%
    summarize(corr = round(mean(!!cor_col), 2)) %>%
              mutate(class = !!class_name) %>%
    select(class, corr)
  
}

calc_pass_rates <- function(df, grouping_terms, pivot_id_cols) {
  
  # this function calculates the bar passage rates of user defined groups
  # and outputs the rates in a table
  
  # note: dataframe must contain a column with bar results called 'result' and 0 msut be fail and 1 msut be pass
  
  # parameters:
  #  df: dataframe with the groups and bar passage data
  #  grouping_terms: terms to group on in calculating bar passage rates
  #                  must include the column that signifies pass or fail, and this column must be called 'result'
  #  pivot_id_cols: the id_cols parameter in pivot_wider, these will be the unique rows in the final table
  
  df %>%
    drop_na(result) %>%
    group_by_at(grouping_terms) %>%
    count() %>%
    ungroup() %>%
    mutate(result = recode(result, `0` = "Fail", `1` = "Pass")) %>%
    pivot_wider(id_cols = pivot_id_cols, names_from = "result", values_from = "n") %>%
    replace_na(list(`Pass` = 0, `Fail` = 0)) %>%
    mutate(pass_rate = Pass / (Pass + Fail))
  
}

# clean data --------------------------------

# calculate difference in values (bar pass rate) from matched dataset
diff_matched <- function(df, match_col, matching_form) {
  
  # parameters:
  #   df: single dataframe for values to match
  #   match_col: the column to match on
  #   matching_form: formula for matching
  matching_class_dataset(df = df, class = match_col, match_formula = matching_form, match_method = 'exact') %>%
    select(class_name, took_class = !!match_col, result, weights) %>%
    # calculate weighted pass rate for each class
    mutate(took_class = recode(took_class, `0` = "not_take", `1` = "took")) %>%
    group_by(took_class) %>%
    summarize(pass_rate = weighted.mean(result, weights, na.rm = T)) %>%
    mutate(class_name = !!match_col) %>%
    # # pivot wider so we can calculate difference
    pivot_wider(id_cols = 'class_name', names_from = 'took_class', values_from = 'pass_rate') %>%
    mutate(pass_diff = took - not_take,
          pass_diff = round(pass_diff*100, 0))
  
}

# import race data
import_race_data <- function(file_path) {
  
  # recode race classification
  race_recode <- c(
    BA = 'Black',
    AI = 'American Indian or Alaskan Native',
    AS = 'Asian or Pacific Islander',
    CW = 'White',
    HL = 'Hispanic/Latinx',
    HS = 'Hispanic/Latinx',
    DI = 'Did not indicate',
    NH = 'Native Hawaiian or Other Pacific Islander',
    OT = 'Other',
    AP = 'Asian or Pacific Islander',
    NO = 'Other',
    OTHER = 'Other',
    `HISP AM` = 'Hispanic/Latinx',
    NI = "Other",
    `AS/BA` = 'Other',
    DIN = 'Other',
    NO = 'Other',
    PR = 'Hispanic/Latinx',
    `A/PI` = 'Asian or Pacific Islander',
    AA = 'Black',
    `AS/CW` = 'Other',
    `CS` = 'White',
    No = 'Did not indicate'
  )
  
  df <- read_csv(file_path) %>%
    mutate(class = str_replace(class, "2016-D", "2017"),
           class = str_replace(class, "2015-D", "2016")) %>%
    # make Feb 19 NC bar takers and later their own jurisdiction since this is when NC transitioned to
    # UBE and bar rates went up a lot
    mutate(jur = ifelse((date %in% c("19F", "19J", "20F")) & (jur == "NC"), "NC_ube", jur),
           race = recode(race, !!!race_recode, .missing = 'Did not indicate'),
           classification = ifelse(race %in% c('White', 'Black'), race, 'Other')) %>%
    # scale and center lsat, u_gpa
    mutate_at(vars(lsat, u_gpa), ~scale(.)) %>%
    mutate(admission_idx = scale(lsat + u_gpa)) %>%
    # scale and center the law school grade columns, but group by year first
    group_by(class) %>%
    mutate_at(vars(cum_gpa, onel_gpa, first_sem_gpa), .funs = list(ctr = ~scale(.))) %>%
    ungroup() %>%
    # center percent rank columns
    mutate_at(vars(contains("_rank")), ~scale(., center = T, scale = F)) %>%
    # remove one visitor
    filter(class != "Visitor")
  
}

# import skill-electives data
import_skills_data <- function(file_path) {
  
  read_csv(file_path) %>%
    # skills elective classes don't start until 2014
    filter(class >= 2014) %>%
    # make missing state rate value (there is only one) 0, which is average
    replace_na(list(state_rate_std = 0)) %>%
    # convert 2016-D class to 2017 and 2015-D to 2015 since these students entered the same year
    mutate(class = str_replace(class, "2016-D", "2017"),
           class = str_replace(class, "2015-D", "2016")) %>%
    # convert class grades to either 0 if missing value (did not take class)
    # or 1 if a grade is present (took class)
    mutate_at(vars(essays:mla), ~as.character(.)) %>%
    mutate_at(vars(essays:mla), ~replace_na(., "Did not take")) %>%
    mutate_at(vars(essays:mla), ~ifelse(. == "Did not take", 0, 1)) %>%
    mutate_at(vars(essays:mla), ~as.numeric(.)) %>%
    # standardize cumulative GPA
    group_by(class) %>%
    mutate(cum_gpa_ctr = scale(cum_gpa)) %>%
    ungroup() %>%
    # duplicate MBE column, so we have an NC only column to iterate through
    mutate(mbeNC = mbe,
           # make Feb 19 and 20 NC bar takers their own jurisdiction since this is when NC transitioned to
           # UBE and bar rates went up a lot
           jur = ifelse((date %in% c("19F", "19J", "20F")) & (jur == "NC"), "NC_ube", jur))
  
}

make_perc_pretty <- function(percentage_col) {
  
  # take a raw percentage (.77) and make it pretty formatted for a table (77%)
  
  # percentage_col = column containing percentages that we want to make pretty
  
  perc <- round(percentage_col, 2) * 100
  perc <- glue("{perc}%")
  
  return(perc)
  
}
