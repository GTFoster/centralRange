`newProject`


**Authors**:
Grant Foster, Cleber Ten Caten, Tad A. Dallas

**Title**: 
The impacts of geographic range size and niche breadth on species centrality in floral visitation networks


**Idea**:
Originally inspired by [this paper](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.13986) by Moulatlet et al (2023), this project looks to see if the breadth of geographic and climatic environments experienced by a species correlate to its role within local interaction networks or larger-scale metawebs.  

**Abstract**:
Understanding the factors that influence species' roles within ecological networks can provide insights into the rules governing the assembly of interaction networks. Centrality, a suite of measures describing how connected a given node is to the rest of a network, is often used to identify species whose loss may result in disproportionate effects on biotic communities. Understanding the drivers of species' centrality in floral visitation networks is crucial, as declines in pollinator abundance and diversity threaten to reduce pollination services globally. Here, we test how geographic range size and niche breadth influence species centrality in floral visitation networks across spatial scales. Using a globally distributed dataset of floral visitation networks comprising 31,677 interactions among 2,703 plant and 4,305 floral visitor species across 99 sites, we model species’ geographic and environmental ranges from large-scale occurrence data to evaluate how these traits shape species’ roles in both local networks and global metawebs. In global metawebs, geographic range size was positively related to closeness centrality. However, after accounting for spatial turnover in interactions, this effect largely disappeared or reversed. Instead, niche breadth showed a positive relationship with the standardized centrality of flowering plant species in both network types. In contrast, in local-scale insect-floral visitor networks, we find significant interactions between geographic range size and niche breadth in the absence of main effects, suggesting multiple ecological pathways through which a species may become highly central in local networks. Together, these results emphasize the need for a multiscale perspective to understand the drivers and consequences of changes in ecological network properties across environments and scales. 

### How to run the analyses

code organization and details on how to reproduce the analyses.

### Repository Structure

```{bash}
├── README.md. 
├── analysis
│   ├── Clean_Analysis.Rmd. This is the main analysis and supplement as presented in the main text of the manuscript. should recreate the entire analysis start to finish provide you install the appropriate packages (easiest way to do this is to set initial chunk eval=TRUE). All visualizations and statistical analyses are created by this script.
manuscript
│   ├── RangeSizing_Function.R This script creates estimates the geographic range size and niche breadth of species in our dataset. Uses species occurence files created by GBIFdb_Querying.R
│   ├── RangeSizing_FillIn_Function.R A modified version of the above code that checks for species we already have estimates for excludes them from new rangesizing runs. (Helpful in batch processing)
│   ├── GBIFdb_Querying.R Script for querying a locally downloaded copy of GBIF to create species occurence record files in a format taken by RangeSizing_Function.R. Requires a local GBIF copy (see the vignette for the package `GBIFdb`)
│   ├── studykey.csv Small csv file used for adding metaData to Web of Life pollination networks
│   ├── backbone_matches.csv Small csv file containing GBIF taxonomy data for species used in our study. Used for GBIFdb querying. 
├── ├── CleanData Intermediate data products created by Clean analysis.RMD
├── ├── CleanFigs Data visualizations output by Clean analysis.RMD. Final publication versions of these are save as .svg files (visualization polishing was done using InkScape)
├── ├── data Input data for our analysis
├── ├── ├──  GBIF folder containing gbif occurence points used in the analysis.
├── ├── worldclim Environmental data used in the calculation of species' niche breadths
├── analysis
```


### Funding

None. Cheap science rules!





