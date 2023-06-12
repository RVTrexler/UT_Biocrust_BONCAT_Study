# BONCAT-FACS-Seq reveals the active fraction of a biocrust community undergoing a wet-up event
Data and analysis files from a study using BONCAT-FACS-Seq to identify active microbes in a biocrust microbial community

Trexler RV, Van Goethem MW, Goudeau D, Nath N, Malmstrom RR, Northen TR, Couradeau E (2023) BONCAT-FACS-Seq reveals the active fraction of a biocrust community undergoing a wet-up event. *Front. Microbiol.*, 14:1176751. https://doi.org/10.3389/fmicb.2023.1176751

## Files

### R analyses
This directory contains all R based analyses

* **Bracken_Analysis.R:** analysis of taxonomic data from metagenomes produced via BONCAT-FACS-Seq.
* **COG_Analysis.R:** analysis of COG data extracted from IMG/M of the metagenomes produced via BONCAT-FACS-Seq.

### Data
This directory contains important data not included in the manuscript supplementary material. The raw sequencing data can be found on the NCBI SRA under BioProject PRJNA938738.

* **Metadata_key.txt:** metadata table describing samples; used in both Bracken_Analysis.R and COG_Analysis.R.
* **bracken_combined_v2.txt:** taxa table output from bracken; analyzed in Bracken_Analysis.R.
* **COG_Data.zip:** zipped directory containing the COG data from each metagenome sample; analyzed in COG_Analysis.R.
* **COG_ID_Annotation.txt:** data table linking COG IDs to associated annotation data; used in COG_Analysis.R.
* **COG_Categories_Expanded.txt:** data table COG Category IDs to associated annotation data; used in COG_Analysis.R.
