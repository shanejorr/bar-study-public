# function that knits the entire report together from various Rmarkdown files
bookdown::render_book('index.Rmd', new_session = FALSE)
