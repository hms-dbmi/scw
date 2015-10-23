#library(knitr)
#knitr::knit2html('analysis_of_heterogeneity_and_subpopulations.Rmd')
# knitr still uses R Markdown v1, which is uglier
# use V2 to integrate Pandoc and Bootstrap themes
library(rmarkdown)
rmarkdown::render('analysis_of_heterogeneity_and_subpopulations.Rmd', 'html_document')