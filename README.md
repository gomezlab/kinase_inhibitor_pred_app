# Kinase Inhibitor Cell Viability Prediction Application

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6835479.svg)](https://doi.org/10.5281/zenodo.6835479)

This repository contains the code related to the Kinase Inhibitor Prediction Application. The application was built using the R shiny framework and uses a model developed in a seperate [project](https://github.com/gomezlab/PRISM_perturbations). You can access the application [here](https://kipred.gomezlab.org/).

![screenshot of kinase inhibitor pred app](application_screenshot.png)

## Running the Application Locally

In addition to the above website, the application can also be run on your personal computer. There are a few files that are too large for Github, so I've made them available through [Zenodo](https://dx.doi.org/10.5281/zenodo.6835479). 

* model_only_500feat_100trees.rds: This is the model file used by the system to make predictions. It should be placed in data/model_only_500feat_100trees.rds (~350 MB)
* matt_model_matrix.h5: This is a subset of the data available through [ARCHS](https://maayanlab.cloud/archs4/) customized to only contain the RNAseq data needed for the model. It should be put in "data/ARCHS_subset/matt_model_matrix.h5", this file is unnecessary if you only plan on running your own RNAseq data through the system. (~160 MB)

After setting up these additional files, the easiest way to start the application is through R Studio by openning the .Rproj file and then clicking the "Run App" button. Alternatively, you can also start the application with this R command:

```{r}
shiny::runApp()
```

I've also made a [`short script`](package_check.R) to check for any missing packages and install them using the pacman library.

### Using the Full 500 Trees Model

The web app uses the 100 tree version of the model, which has nearly identical performance to the 500 tree model, but has 1/5 the RAM requirements. If you are interested in running the full 500 trees model, we provide the trained model file from the above zenodo link. In order to use this file instead of the 100 trees model, you will need to modify the source code. Search the code in the repository for "model_only_500feat_100trees.rds" and replace with 
"model_only_500feat_500trees.rds". You should only find a single instance of the model file in the [`functions.r`](functions.R) file.

## Reproducibility

Since this work is largely based around making an applicaion available to use a previously described model, there isn't much to say here. We do describe producing a reduced model to make it feasible to make the service available and the methods to reproduce this is described in another [github repository](https://github.com/gomezlab/PRISM_perturbations).
