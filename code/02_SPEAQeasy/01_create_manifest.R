library("here")
library("sessioninfo")

## Read in information from the Synapse download query
synapse <- read.csv(here("raw-data/psychENCODE/v2023-07-14/RNAseq/SYNAPSE_TABLE_QUERY_127308191.csv"))
## Subset to just RNAseq data
synapse <- subset(synapse, fileFormat == "fastq")
dim(synapse)
# [1] 7838   46

## Check that individualID and specimenID do not overlap
table(sapply(with(synapse, split(individualID, specimenID)), function(x) length(unique(x))))
#   1
# 516

## Check that all fastq.gz files were downloaded
synapse$JHPCE_path <- file.path(here("raw-data/psychENCODE/v2023-07-14/RNAseq"), synapse$name)
synapse$JHPCE_file_exists <- file.exists(synapse$JHPCE_path)
table(synapse$JHPCE_file_exists)
# FALSE  TRUE
#   227  7611
synapse_missing <- subset(synapse, !JHPCE_file_exists)
synapse_present <- subset(synapse, JHPCE_file_exists)

unique(synapse_missing$specimenID)
#  [1] "LBPSEMA4BRAIN639" "LBPSEMA4BRAIN730" "LBPSEMA4BRAIN176" "LBPSEMA4BRAIN285"
#  [5] "LBPSEMA4BRAIN469" "LBPSEMA4BRAIN043" "LBPSEMA4BRAIN141" "LBPSEMA4BRAIN085"
#  [9] "LBPSEMA4BRAIN607" "LBPSEMA4BRAIN401" "LBPSEMA4BRAIN348" "LBPSEMA4BRAIN478"
# [13] "LBPSEMA4BRAIN402" "LBPSEMA4BRAIN778" "LBPSEMA4BRAIN623"
unique(synapse_missing$specimenID) %in% synapse_present$specimenID
#  [1]  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# [13] FALSE FALSE FALSE

## Official syntax from https://adknowledgeportal.synapse.org/Explore/Studies/DetailsPage/StudyData?Study=syn26337520
# synapse get -q "SELECT * FROM syn11346063.39 WHERE ( \"study\" HAS ( 'LBP' ) ) AND ( \"specimenID\" LIKE '%LBPSEMA4BRAIN639%' OR \"specimenID\" LIKE '%LBPSEMA4BRAIN730%' )"
paste(paste0("\"specimenID\" LIKE '%", unique(synapse_missing$specimenID), "' OR"), collapse = "")
# "\"specimenID\" LIKE '%LBPSEMA4BRAIN639' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN730' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN176' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN285' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN469' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN043' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN141' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN085' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN607' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN401' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN348' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN478' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN402' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN778' OR\"specimenID\" LIKE '%LBPSEMA4BRAIN623' OR"

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
