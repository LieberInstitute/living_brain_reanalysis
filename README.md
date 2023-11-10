# living_brain_reanalysis

[![DOI](https://zenodo.org/badge/670288774.svg)](https://zenodo.org/doi/10.5281/zenodo.10010510)

Re-analysis of the data from Liharska et al. [10.1101/2023.04.21.23288916](https://doi.org/10.1101/2023.04.21.23288916). 

> A study of gene expression in the living human brain
> Lora E. Liharska, You Jeong Park, Kimia Ziafat, Lillian Wilkins, Hannah Silk, Lisa M. Linares, Ryan C. Thompson, Eric Vornholt, Brendan Sullivan, Vanessa Cohen, Prashant Kota, Claudia Feng, Esther Cheng, Jessica S. Johnson, Marysia-Kolbe Rieder, Jia Huang, Joseph Scarpa, Jairo Polanco, Emily Moya, Alice Hashemi, Matthew A. Levin, Girish N. Nadkarni, Robert Sebra, John Crary, Eric E. Schadt, Noam D. Beckmann, Brian H. Kopell, Alexander W. Charney
> medRxiv 2023.04.21.23288916; doi: https://doi.org/10.1101/2023.04.21.23288916

The results published here are in whole or in part based on data obtained from the [AD Knowledge Portal](https://adknowledgeportal.org/). Specifically, we downloaded the data (FASTQ raw files and phenotype data tables) from [Synapse (syn26337520)](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520).

# Processed-data

* The gene level LBP dataset re-processed by _SPEAQeasy_ is available at [processed-data/02_SPEAQeasy/pipeline_output/count_objects/rse_gene_living_brain_reanalysis_n516.Rdata](processed-data/02_SPEAQeasy/pipeline_output/count_objects/rse_gene_living_brain_reanalysis_n516.Rdata). It is available as a [_SummarizedExperiment_](https://bioconductor.org/packages/SummarizedExperiment/) R/Bioconductor object. You will need the sample demographic information that is access-controlled and available from [Synapse (syn26337520)](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520) to fully reproduce our analyses.
* The transcript level LBP dataset re-processed by _SPEAQeasy_ is available at [processed-data/02_SPEAQeasy/pipeline_output/count_objects/rse_tx_living_brain_reanalysis_n516.Rdata](processed-data/02_SPEAQeasy/pipeline_output/count_objects/rse_tx_living_brain_reanalysis_n516.Rdata). Given that this is a 807 MB file, it is available with [GIT LFS](https://git-lfs.com/). You will need the sample demographic information that is access-controlled and available from [Synapse (syn26337520)](https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520) to fully reproduce our analyses.
* The transcript level RNA degradation dataset from Jaffe et al., _PNAS_, 2017 processed by _recount3_ is available at [processed-data/05_check_degradation/rse_degrade_tx.Rds](processed-data/05_check_degradation/rse_degrade_tx.Rds). It was created [here](https://github.com/LieberInstitute/living_brain_reanalysis/blob/00e469cb1e029717cd6ca407ea6f06e7d0bb8100/code/05_check_degradation/01_combine_recount3_salmon_output.R#L74-L82).

# Supplementary Tables

See [processed-data/SupplementaryTables](processed-data/SupplementaryTables/) for all supplementary tables. The descriptions include GitHub permalinks to the scripts that made them.

# Code structure

In this re-analysis we used [`SPEAQeasy`](https://doi.org/10.1186/s12859-021-04142-3) to re-align the FASTQ files to GRCh38. The specific version `SPEAQeasy` we used is available at [LieberInstitute/SPEAQeasy/tree/living_brain_reanalysis](https://github.com/LieberInstitute/SPEAQeasy/tree/living_brain_reanalysis).

Files are organized following the structure from [LieberInstitute/template_project](https://github.com/LieberInstitute/template_project). Log files include the corresponding R session information with details about version numbers of the packages we used.

# Citation

We hope that this repository will be useful for your research. Please use the following [BibTeX](https://en.wikipedia.org/wiki/BibTeX) information to cite this code repository as well as our re-analysis of the Liharska et al Living Brain Study dataset. Thank you!

> **Comparison of gene expression in living and postmortem human brain**
> Leonardo Collado-Torres, Lambertus Klei, Chunyu Liu, Joel E Kleinman, Thomas M Hyde, Daniel H Geschwind, Michael J Gandal, Bernie Devlin, Daniel R Weinberger
> medRxiv 2023.11.08.23298172; doi: <https://doi.org/10.1101/2023.11.08.23298172>

```
@article {Collado-Torres2023.11.08.23298172,
	author = {Leonardo Collado-Torres and Lambertus Klei and Chunyu Liu and Joel E Kleinman and Thomas M Hyde and Daniel H Geschwind and Michael J Gandal and Bernie Devlin and Daniel R Weinberger},
	title = {Comparison of gene expression in living and postmortem human brain},
	elocation-id = {2023.11.08.23298172},
	year = {2023},
	doi = {10.1101/2023.11.08.23298172},
	publisher = {Cold Spring Harbor Laboratory Press},
	URL = {https://www.medrxiv.org/content/early/2023/11/09/2023.11.08.23298172},
	eprint = {https://www.medrxiv.org/content/early/2023/11/09/2023.11.08.23298172.full.pdf},
	journal = {medRxiv}
}
```

Please note that this project was only made possible thanks to many other R and bioinformatics software authors, which are cited either in the vignettes and/or the paper(s) describing this package.

# Internal

JHPCE location: `/dcs05/lieber/lcolladotor/living_brain_LIBD001/living_brain_reanalysis`
