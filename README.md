# ROKET

## Introduction

This package is designed to perform optimal transport and hypothesis testing on kernel matrices when performing kernel regression. The software is optimized for calculating distance between pairs of samples based on the pairs of mutated gene statuses.

## Installation

```R
# Dependencies
req_packs = c("devtools","Rcpp","RcppArmadillo","reshape2",
	"smartr","MiRKAT","ggplot2","ROKET")
all_packs = as.character(installed.packages()[,1])
rerun = 0

for(pack in req_packs){
	if( pack %in% all_packs ){
		library(package = pack,character.only = TRUE)
		next
	}

	if( pack %in% c("smartr","ROKET") ){
		repo = sprintf("pllittle/%s",pack)
		bb = tryCatch(devtools::install_github(repo = repo,
			build_vignettes = TRUE),
			error = function(ee){"error"})
	} else {
		bb = tryCatch(install.packages(pack),
			error = function(ee){"error"})
	}
	
	if( !is.null(bb) && bb == "error" )
		stop(sprintf("Error for package = %s",pack))
	rerun = 1
}

if( rerun == 1 ) stop("Re-run above code")
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
