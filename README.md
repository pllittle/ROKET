# ROKET

## Introduction

This package is designed to perform optimal transport and hypothesis testing on kernel matrices when performing kernel regression. The software is optimized for calculating distance between pairs of samples based on the pairs of mutated gene statuses.

## Installation

```R
# Dependencies
req_packs = c("devtools","Rcpp",
	"RcppArmadillo","MiRKAT","ggplot2")
all_packs = as.character(installed.packages()[,1])

for(pack in req_packs){
	if( pack %in% all_packs ) next
	stop(sprintf("Install R package = %s",pack))
}

# Other dependencies
if( !("smartr" %in% all_packs) )
	devtools::install_github("pllittle/smartr")

# Install
if( !("ROKET" %in% all_packs) )
	devtools::install_github("pllittle/ROKET")
```

By default, the software runs a single thread and loops through all pairs of samples for distance calculations. However if OpenMP is installed, the user can make use of multi-threaded calculations.

## Vignette

```R
# An Introduction
vignette(topic = "intro",package = "ROKET")
```

## Citation
Little, P., Hsu, L., Sun, W. (2021). ROKET: Associating Somatic Mutation with Clinical Outcomes through Kernel Regression and Optimal Transport. *bioRxiv*. [[HTML](https://www.biorxiv.org/content/10.1101/2021.12.23.474064v1), [PDF](https://www.biorxiv.org/content/10.1101/2021.12.23.474064v1.full.pdf)]

## Workflow

R package and code to perform the manuscript's workflow are provided [here](https://github.com/pllittle/ROKETworkflow).
