
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DSFworld

## Overview

DSFworld is a web application and set of supporting scripts for the
analysis and interactive modeling of Differential Scanning Fluorimetry
data. The full DSFworld web application can be access at the url:
<https://gestwickilab.shinyapps.io/dsfworld/>.

This repository contains the full code for the DSFworld web application,
as well as the associated scripts, and modularized mini web application
for each of the individual tasks handled at DSFworld (uploading,
layouts, plotting, analysis, and downloads.)

  - 0_full_app contains the full DSFworld application
  - 1_upload_raw contains the data uploading the pre-formatting module 
  - 2_add_layouts contains the data labeling module
  - 3_plot contains the interactive plotting module
  - 4_analyze contains the Tma analyses (dRFU, fits 1-4, and BIC) module 
  - 5_download contains the data and results downloading module
  - The interactive DSF data model can be found in the 0_full_app > scripts > data_modeling.R

The DSFworld website is presented as a part of a larger publication
entitled “Three Essential Resources to Improve Differential Scanning
Flourimetry (DSF) Experiments”. The full text of this paper and its
supplementary information can be downloaded from the DSFworld website,
as well as from the 0\_full\_app folder of this repository.

## Usage

DSFworld and the associated scripts use the following packages:

  - [shiny](https://shiny.rstudio.com/), for building a web
    interface.

  - [shinyBS](https://cran.r-project.org/web/packages/shinyBS/shinyBS.pdf),
    for drop-down
    panels.

  - [shinyalert](https://cran.r-project.org/web/packages/shinyalert/shinyalert.pdf),
    for pop-up
    messages.

  - [shinycssloaders](https://cran.r-project.org/web/packages/shinycssloaders/shinycssloaders.pdf),
    for busy
    spinners.

  - [rhandsontable](https://cran.r-project.org/web/packages/rhandsontable/rhandsontable.pdf),
    for editable tables.

  - [tidvyerse](https://www.tidyverse.org/), for data handling and
    visualization.

  - [modelr](https://cran.r-project.org/web/packages/modelr/modelr.pdf),
    for data and model
    handling.

  - [minpack.lm](https://cran.r-project.org/web/packages/minpack.lm/minpack.lm.pdf),
    for fitting.

  - [signal](https://cran.r-project.org/web/packages/signal/signal.pdf),
    for Savitsky-Golay
    filtering.

  - [quantmod](https://cran.r-project.org/web/packages/quantmod/quantmod.pdf),
    for assistance starting parameter
    estimates.

  - [SciViews](https://cran.r-project.org/web/packages/SciViews/SciViews.pdf),
    for natural logarithms.

## Contact

You can contact us at dsfworlducsf (at) gmail (dot) com.
