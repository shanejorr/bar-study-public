# Elon School of Law Bar Passage Evaluation

## Project structure

This project contains all code used for an evaluation of Elon School of Law's bar passage rates. All analyses were completed in Rmarkdown files, with one file correpsonding to one section of the report. The Rmarkdown file for each section can be found in the following locations:

- *Race and Bar Passage:*  `analysis/race/01-race.Rmd`
- *Required Bar-related Courses and Bar Passage:* `analysis/doctrinal-required/02-doctrinal-required.Rmd`
- *Skills Electives and Bar Passage:* `analysis/skills-electives/03-skills-elective.Rmd`
- *Doctrinal Electives and Bar Passage:* `analysis/doctrinal-elective/04-skills-elective.Rmd`
- *Commercial Prep Courses and Bar Passage:* `analysis/commercial-prep/05-commercial-prep.Rmd`

Custom functions are used throughout the project. The file `analysis/analysis_functions.R` contains functions used in all sections. Additionally, two sections have custom functions used only in that section. These include:

- *Race and Bar Passage:*  `analysis/race/race_custom_functions.R`
- *Required Bar-related Courses and Bar Passage:* `analysis/doctrinal-required/doctrinal-required-functions.R`

## Purpose of this repo

This repo does not contain any data due to student privacy concerns. Therefore, users will not be able to run the code. But, the code allows others to inspect the models used in the report and provides an example of creating an entire report in `Rmarkdown` and `bookdown`.