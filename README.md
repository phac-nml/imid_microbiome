
## 16S rRNA gene targeted-amplicon data processing ##
Input: 16S rRNA sequence data<br />
Scripts:<br />
 
  * sbatch_make_contigs.sh   <br />
  * make_contigs.sh  <br />
  * merge_fasta.sh  <br />
  * groups_20dec2017.sh   <br />
  * analysis_20dec2017.sh  <br />
  * summary.single.sh  <br />
  

## Post-OTU analyses: normalization and machine learning classification##
Inputs: An OTU count table with taxanomic information for the OTUs and sample information.

Scripts:<br />
  * microbiome_1stStep.r: data wrangling  <br />
  * microbiome_2ndStep.r: data filtering, normalization and machine learning classification  <br />
  * ***_batch.r files: sourced in the previous two files to do analyses in batch  <br />

