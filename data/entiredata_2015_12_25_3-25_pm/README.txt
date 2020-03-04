This folder contains files downloaded from the CisBP-RNA database, which
contains information on RNA binding proteins (RBPs), and their motifs.
The folder itself is named after the subset of data you chose to download,
with the date and time it was downloaded appended at the end. Depending on what 
you chose to download, this directory contains a subset of the following files:


logos_all_motifs - directory containing .png images in the form of sequence
logos, which visually depict the sequence binding preferences of a given RBP
obtained from a single experimental source.  Each file is named by its
Motif_ID (see description of the text files below).


pwms_all_motifs - directory containing text files giving the frequency of each
base at each position in the motif.  Each file is named by its Motif_ID (see
below).


Escores/Zscores.txt - these text files are matrices providing the E- and
Z-scores from each RNAcompete assay, for all 16,382 possible 7-mers.  E-scores
range from -0.5 to +0.5; 7-mers with E-scores >= 0.45 are considered "strongly
bound" by the associated RBP.  Z-scores can have arbitrary ranges, with higher
scored indicating stronger binding.  See Ray et al. Nature 2013 for additional
details.


RBP_Information.txt and RBP_Information_all_motifs.txt -
These files contain information on the RBPs.  'RBP_Information.txt' contains,
for each RBP, all directly determined motifs (see below). If an RBP does not 
have a directly determined motif, this file will also include its best
inferred motif.  'Best' is defined as the motif(s) obtained from the most
similar RBP (based on the %ID in the amino acids of its RBP) that has a directly 
determined motif.

RBP_Information_all_motifs.txt is a superset of 'RBP_Information.txt'.  It
also includes any motif that can be inferred for a given RBP, up to an
RBD %ID of 50%.  Note that %IDs below 70% often have similar motifs, but there
are a greater number of cases where their motifs begin to deviate. See Ray et al. 
Nature 2013 for additional information. 

The files have identical formats. Each file contains one row per RBP.  RBPs with 
multiple motifs have one line per motif. The columns in the file are as follows:

RBP_ID	Internal CisBP-RNA ID for the RNA binding protein.  Each gene has a unique RBP_ID.
Family_ID	Internal CisBP-RNA ID for the RBP family.  A family is the unique set of RNA binding domains (RBDs) present in the protein.
RSource_ID	Internal CisBP-RNA ID for the source of the RNA binding protein (i.e. where its genome sequence was obtained).
Motif_ID	Internal CisBP-RNA ID for the associated motif.
MSource_ID	Internal CisBP-RNA ID for the source of the motif (i.e. which database or study it came from)
DBID	External ID of the RBP (e.g., Ensembl ID)
RBP_Name	Name of the RBP
RBP_Species	Species of the RBP
RBP_Status	Motif status of the RBP: 'D' stands for directly determined motif.  'I' indicates that the motif is inferred from another RBP, based on RBD similarity (see Ray et al. Nature 2013 for details). 'N' means no motif is available.
Family_Name	Name of the RBP's family
RBDs	The unique set of RBDs (Pfam names) present in the RBP
RBD_Count	Number of unique RBDs in the RBP
Cutoff	Cutoff used to infer motifs for the RBP family
DBID	Motif ID from the associated database or study
Motif_Type	Experimental assay used to determine the motif
MSource_Identifier	ID for the source of the motif (i.e., its project name)
MSource_Type	Internal CisBP-RNA ID for the motif category
MSource_Author	First author for the source of the motif
MSource_Year	Year of publication of the motif source
PMID	Pubmed ID of the motif source
MSource_Version	Version of the source (i.e. database build)
RBPSource_Name	Source of the RBP (i.e. where did the genome build come from?)
RBPSource_URL	URL of the RBP source
RBPSource_Year	Year the genome data was downloaded
RBPSource_Month	Month the genome data was downloaded
RBPSource_Day	Day the genome data was downloaded

Questions can be directed to Matt Weirauch: mattweirauch@gmail.com
