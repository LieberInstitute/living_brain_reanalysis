# living_brain_reanalysis

[![DOI](https://zenodo.org/badge/670288774.svg)](https://zenodo.org/doi/10.5281/zenodo.10010510)

Re-analysis of the data from Liharska et al. [10.1101/2023.04.21.23288916](https://doi.org/10.1101/2023.04.21.23288916). 

> A study of gene expression in the living human brain
> Lora E. Liharska, You Jeong Park, Kimia Ziafat, Lillian Wilkins, Hannah Silk, Lisa M. Linares, Ryan C. Thompson, Eric Vornholt, Brendan Sullivan, Vanessa Cohen, Prashant Kota, Claudia Feng, Esther Cheng, Jessica S. Johnson, Marysia-Kolbe Rieder, Jia Huang, Joseph Scarpa, Jairo Polanco, Emily Moya, Alice Hashemi, Matthew A. Levin, Girish N. Nadkarni, Robert Sebra, John Crary, Eric E. Schadt, Noam D. Beckmann, Brian H. Kopell, Alexander W. Charney
> medRxiv 2023.04.21.23288916; doi: https://doi.org/10.1101/2023.04.21.23288916

The results published here are in whole or in part based on data obtained from the [AD Knowledge Portal](https://adknowledgeportal.org/). Specifically, we downloaded the data (FASTQ raw files and phenotype data tables) from [Synapse (syn26337520)](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520).

# Supplementary Tables

See [processed-data/SupplementaryTables](processed-data/SupplementaryTables/) for all supplementary tables. The descriptions include GitHub permalinks to the scripts that made them.

# Code structure

In this re-analysis we used [`SPEAQeasy`](https://doi.org/10.1186/s12859-021-04142-3) to re-align the FASTQ files to GRCh38. The specific version `SPEAQeasy` we used is available at [LieberInstitute/SPEAQeasy/tree/living_brain_reanalysis](https://github.com/LieberInstitute/SPEAQeasy/tree/living_brain_reanalysis).

Files are organized following the structure from [LieberInstitute/template_project](https://github.com/LieberInstitute/template_project). Log files include the corresponding R session information with details about version numbers of the packages we used.

# Citation

We hope that this repository will be useful for your research. Please use the following [BibTeX](https://en.wikipedia.org/wiki/BibTeX) information to cite this code repository as well as the commentary paper. Thank you!

```
@article {TODO
}
```

Please note that this project was only made possible thanks to many other R and bioinformatics software authors, which are cited either in the vignettes and/or the paper(s) describing this package.

# Internal

JHPCE location: `/dcs05/lieber/lcolladotor/living_brain_LIBD001/living_brain_reanalysis`
