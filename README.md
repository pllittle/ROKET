<div align="left">
<a href=""><img src="https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink&label=ROKET" height="80" /></a>
</div>

<!-- badges: start -->
![C++](https://img.shields.io/badge/C++-%2300599C.svg?style=square&logo=c%2B%2B&logoColor=gold)
![R](https://img.shields.io/badge/R-%23276DC3.svg?style=square&logo=r&logoColor=pink)
![CRAN status](https://www.r-pkg.org/badges/version/ROKET)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2021.12.23.474064.svg)](https://doi.org/10.1101/2021.12.23.474064)
<!-- badges: end -->

This package is designed to perform optimal transport and hypothesis 
testing on kernel matrices when performing kernel regression. The 
software is optimized for calculating distance between pairs of 
samples based on the pairs of mutated gene statuses.

## Installation

<details>

<summary>Click to expand!</summary>

```R
# Dependencies
req_packs = c("devtools","Rcpp","RcppArmadillo","reshape2",
	"ggdendro","smarter","MiRKAT","ggplot2","ROKET")
all_packs = as.character(installed.packages()[,1])
rerun = 0
build_vign = ifelse(Sys.getenv("RSTUDIO_PANDOC") == "",FALSE,TRUE)

for(pack in req_packs){
	if( pack %in% all_packs ){
		library(package = pack,character.only = TRUE)
		next
	}
	
	bb = NULL
	if( pack %in% c("smarter","ROKET") ){
		repo = sprintf("pllittle/%s",pack)
		bb = tryCatch(devtools::install_github(repo = repo,
			build_vignettes = build_vign,
			dependencies = TRUE),
			error = function(ee){"error"})
	} else {
		bb = tryCatch(install.packages(pkgs = pack,
			dependencies = TRUE),
			error = function(ee){"error"})
	}
	
	if( !is.null(bb) && bb == "error" )
		stop(sprintf("Error for package = %s",pack))
	rerun = 1
}

if( rerun == 1 ) stop("Re-run above code")
```

By default, the software runs a single thread and loops through all 
pairs of samples for distance calculations. However if OpenMP is 
installed, the user can make use of multi-threaded calculations.

</details>

## Vignette

```R
# An Introduction
vignette(topic = "intro",package = "ROKET")
```

## Citation
Little, P., [Hsu, L.](https://www.fredhutch.org/en/faculty-lab-directory/hsu-li.html), 
[Sun, W.](https://github.com/sunway1999) (2021). ROKET: Associating 
Somatic Mutation with Clinical Outcomes through Kernel Regression 
and Optimal Transport. *bioRxiv*. 
[[HTML](https://www.biorxiv.org/content/10.1101/2021.12.23.474064v1), 
[PDF](https://www.biorxiv.org/content/10.1101/2021.12.23.474064v1.full.pdf)]

## Workflow

R package and code to perform the manuscript's workflow are 
provided [here](https://github.com/pllittle/ROKETworkflow).
