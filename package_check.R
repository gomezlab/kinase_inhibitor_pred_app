#!/usr/bin/env Rscript

# to find all the library commands run:
#  grep -Rh library * | sort | uniq
# 
# then reformat the library calls to use p_load as below, plus dealing with the github only packages

if("pacman" %in% rownames(installed.packages()) == FALSE) {
	install.packages("pacman")
}

library(pacman)

p_load(shiny)
p_load(tidyverse)
p_load(here)
p_load(tidymodels)
p_load(shinyjs)
p_load(digest)
p_load(rmarkdown)
p_load(rhdf5)

p_load_gh('mbergins/BerginskiRMisc')
p_load_gh("daqana/dqshiny")
