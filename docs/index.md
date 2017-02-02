# Before you start
The following tools are required for running VarSelect, please read carefully and install all the following tools for full benefits. Please note that the current version of VarSelect supports only human genome reference hg19/GRCh37, and stay in touch for GRCH38 update. 

## Gemini

Gemini is a framework for human genetic variations analysis by the SQLite database engine, and can be downloaded and installed by following the instruction at https://gemini.readthedocs.io/en/latest/content/installation.html . Once installed, please follow the commands at https://gemini.readthedocs.io/en/latest/content/installation.html#updating-the-gemini-executables-and-annotations to download the data file of the additional annotations, such as the GERP and CADD scores.
 
## VEP
VEP is part of the Ensembl project and is a comprehensive variants annotation tool. To install VEP, please follow the instruction at http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html. Please note that the pre-computed cache files are required to speeding up the annotation process. Please follow the instruction step by step at http://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html to download the VEP of GRCh37 version.

## VEP plugins
VEP supports plugin modules to incorporate annotations from the external datasets. VarSelect requires the dbNSFP plugin, which are available at https://github.com/Ensembl/VEP_plugins/blob/release/86/dbNSFP.pm, respectively. 

## dbNSFP
dbNSFP annotates the non-synonymous single nucleotide variations (nsSNVs). The data file of dbNSFP is available at https://sites.google.com/site/jpopgen/dbNSFP .

## ANNOVAR
ANNOVAR  is a variant annotation tool with high efficiency to a variety of annotation databases, and is available at http://annovar.openbioinformatics.org/en/latest/user-guide/download/. Please note that a license is required. Please follow the instruction at http://annovar.openbioinformatics.org/en/latest/user-guide/startup/ to install scripts into proper directories when all the required packages are downloaded. Databases will be automatically installed by VarSelect installation script.

## snpEff
snpEff annotates and predicts the impact of genetic variants, and is available at http://snpeff.sourceforge.net/download.html. After downloading the software, the pre-built snpEff database for human genome of GRCh37 is needed. Please download it with following command:

```
java –jar /path/to/your/snpEff.jar download –v GRCh37.75
```
## vcftools
vcftools are a set of tools for manipulating genetic variants in the VCF-formatted files, and are available at https://vcftools.github.io/index.html. Please follow the instruction to install vcftools at https://vcftools.github.io/examples.html.

## bcftools, bgzip, tabix
bcftools, bgzip and tabix are tools to compress, index and manipulate VCF files. bcftools are available at http://www.htslib.org/download/, and include the bgzip and the tabix tools in the release.

# Download VarSelect
The latest version of VarSelect package can be downloaded from https://github.com/VarSelect/VarSelect

# Install VarSelect
Please make sure that you have already downloaded all the required packages and resources for running VarSelect. When you are all set, please run the following command to decompress the VarSelect file.

```
tar zxvf VarSelect-latest.tar.gz
```

After extracting the package, run the VarSelect installation script

```
/path/to/your/VarSelect/install_VarSelect.pl
```

Add the VarSelect path to your system's $PATH settings

```
export PATH=/path/to/your/VarSelect/dir:$PATH
```


# Quick start

VarSelect script is executable only on command line. Please use -h flag for the basic usage information. 

```
varselect.pl –h
```

VarSelect annotates and analyzes sequence variants, and compares the results from different primary analyses. To start using VarSelect, please use "annotate" to process your vcf file(s) of interests.

```
varselect.pl annotate -v /path/to/vcf_files_list \
                      -p /path/to/ped/file \
                      -m workflow_mode
```

The annotation function combines the variants from different samples and comprehensively annotates all reported variants. The -v option specifies the file recording the links between samples and corresponding variant files. The links is specified by a comma separator, one file per line as the following format:

```
sample1,/path/to/vcf/file1
sample1,/path/to/vcf/file2
sample2,/path/to/vcf/file3 
```

The gender and phenotype information is also required, and is specify by -p option followed by a PED file containing the pedigree information. An example of the PED file is available at http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped

Once annotated, VarSelect creates a SQLite database file by using the GEMINI framework. The database can be queried and filtered by using the GEMINI commands to extract the annotation information.

```
gemini query —header -q 'select * from variants limit 10' \
/path/to/varselect.db
```

Besides annotation, VarSelect also integrates copy number information to the variants. The option –c specifies the relationship between sample and copy number variation call from CNVkit. The relationship is specified in the file by comma separated pair, one file per line as following format:

```
sample1,/path/to/cns/file1
sample2,/path/to/cns/file2
```

VarSelect will annotate with tags: cnv_samples, cnv_log2. Quantitative change by log base 2 between the paired samples is computed and annotated on the cnv_fc_samples and cnv_foldchange_log2 tags. 

VarSelect also incorporate user's gene expression profile to annotate variants and supports the output files from DEXseq, sailfish or featureCounts. VarSelect annotates the variants on tags: xpr_samples, xpr_readcounts, xpr_fpkm, xpr_tpm.

```
sample1,readcount,/path/to/file/of/reads_counts
sample1,exoncount,/path/to/file/of/dexseq_generated_exoncounts
sample1,fpkm,/path/to/combined.fpkm
sample1,tpm,/path/to/combined.gene.sf.tpm
sample2,readcount,/path/to/file/of/reads_counts
sample2,exoncount,/path/to/file/of/dexseq_generated_exoncounts
sample2,fpkm,/path/to/combined.fpkm
sample2,tpm,/path/to/combined.gene.sf.tpm
```

Quantitative change between the paired sample is pairwise computed and annotated on tags: xpr_fc_samples, xpr_foldchange_readcount, xpr_foldchange_fpkm , xpr_foldchange_tpm

## Built-in analytic workflow for family and paired sample analysis

VarSelect provides two basic analytic workflows: 1) family analysis or 2) paired case-control analysis, by specifying the option “-m family” or “-m paired”, respectively.

For family analysis, VarSelect analyzes the genetic variants for five genetic models, including autosomal recessive (AR), compound heterozygosity (CH), de novo recessive (DNR), two-hit recessive (THR) and X-link recessive (XLR). The variants follow each of the model will be labelled on the tags: is_AR, is_CH, is_DR, is_TH and is_XL, respectively.

For paired case-control analysis, VarSelect labels the changes of nucleotides between the paired samples with 'loss of heterozygous (LOH)' or 'de novo somatic mutations'. The changes are labelled on tags: is_LOH and is_de_novo. 
The result of each analysis is recorded by a tag in_analysis_jobid, while the 'jobid' is the time stamp when the analysis is performed by VarSelect.

## Start from vcf files by multiple callers
Different variant callers often result in inconsistent variant calling reports while mostly are consistent with some inconsistent variants. VarSelect is engineered to deal with such situation by processing multiple VCF files by different variants callers in two ways: 1) unify all reported variants or 2) intersect variants reported by all callers. The option -k flag triggers this function, followed by option -u (union) or option -i (intersection), depending on the analytic purpose.

```
varselect.pl annotate -v /path/to/vcf/files/list \
                      -p /path/to/ped/file \
                      -m workflow_mode  \
                      -k -i

```

The format of the list of VCF files is different from single-caller mode since now a sample would have multiple VCF files. User must specify sample, variant caller, and the associated VCF file. An extra comma separated field with the name of variant caller is added to the end of each line as follows:

```
sample1,/path/to/vcf/file1,caller1
sample1,/path/to/vcf/file2,caller2
sample2,/path/to/vcf/file3,caller1
sample2,/path/to/vcf/file4,caller2
```

Please note that regardless filtering for the union or intersection, inconsistent calls among different callers are regarded ambiguity and are marked and removed from further analysis. The list of removed variants will be stored in the result directory.

## Re-analysis and update the analysis database 
Samples in a complex study design can be annotated together and analyzed according to various analytic purpose. The label of the phenotypic information (e.g. tumor and normal; affected and unaffected) of the samples can be rewritten and thereafter analyzed accordingly, according to the labels specified in the PED file. Samples begun with a hash character '#' in PED file are excluded from the downstream analysis. An example PED file is shown here, and the samples “uncle” and “aunt” are excluded from the downstream analysis

```
#family_id  sample_id   paternal_id maternal_id    sex     phenotype
family1     father       0           0         1       1
family1     mother      0           0         2       1
family1     daughter    father      mother      2       2
#family1    uncle       0           0         2       2
#family1    aunt        0           0         2       2
```

An initial full primary generates a VarSelect database, which is required for recording the results of the subsequent re-analysis. The -d option specifies the location of where the database file is kept. 

```
varselect.pl analysis -d /path/to/gemini/db \
                      -p /path/to/modified/ped/file \
                      -m workflow mode \
                      -k -u
```
After re-analysis, a new analysis directory will be created with new jobid. The results and logs will be stored in new directory. Filtered variants will be assigned new tag in_analysis_new_jobbid in the VarSelect database.

## Secondary analysis: compare results from any two primary analysis 

Analysis in full or in part of the annotation, single or multiple variant callers, and analytic workflow is termed ‘primary analysis’, while the results from different primary analysis can be compared for specific purposes and is termed ‘secondary analysis’. Comparison between any two primary analyses is specified by the options “–a” and “–b” as follows.

```
varselect.pl compare  -a /dir/to/analysis/of/datasetA \
                   -b /dir/to/analysis/of/datasetB \
                   -c [1-4] 
```

The secondary analysis includes four comparisons including: 1) the union of analysis A and B. 2) intersection of analysis A and B. 3) variants presented in the analysis A but not in the analysis B. 4) variants present in the analysis B but not in the analysis A. 

Results of new secondary analysis will be stored in a new directory. The filtered variants will be assigned a tag in_analysis with the Job ID.


# Description of VarSelect scripts
VarSelect is composed of many individual scripts. Bellows are the description of each VarSelect script. 
* varselect.pl is the main script of VarSelect. It provides three commands included in VarSelect: annotate (initial annotation), analysis (primary and re-analysis) and compare (secondary analysis).
  * Command “annotate” triggers the script vs_annotate.pl to annotate VCF files from scratch and the script vs_analysis.pl to analyze variants through workflow of choice. There are three required options: -v sample-vcf file list, -p PED file, -m workflow mode.
  * Command “analysis” triggers the script vs_analysis.pl. There are three required options: -d gemini db file, -p PED file, -m workflow mode. 
  * Command “compare” triggers the script vs_compare.pl to compare results between two primary analysis. There are three required options: “-a” and “-b” to specify the path of the analysis A and B. “-c” specify the method of comparison (1. union, 2. intersection, 3. A only, and 4. B only)

* vs_annotate.pl is called by varselect.pl for three functions. 
  * Firstly, it prepares VCF file for annotation. VCF files from same sample will be joined together by vcf-concat of VCFtools. VCF files of different samples are then merged into a VCF file by vcf-merge in the VCFtools. The genotypes of sex chromosome on merged VCF file are curated by vcf-fix-ploidy in the VCFtools.
  * Secondly, the script triggers VEP, snpEff and ANNOVAR for annotation
  * Thirdly, the script triggers gemini framework to generate a sqlite db for downstream analysis.


# Examples

## Examples 1 samples from a family 

In this example, we use variants from chromosome 22 of family NA12878,NA12891,NA12892 as an example for family-based analysis.

All required files are stored in directory varselect/examples/example1/ :

* NA12878_chr22.vcf.gz  : Variant calls from chromosome 22 of these genomes. 
* NA12891_chr22.vcf.gz 
* NA12892_chr22.vcf.gz 
* example1.txt : A comma separated file that describe the relationship between sample and corresponding VCF file.  
* example1.ped : A tab separated file that describe the relationship between each sample, the sex of samples and the affected status of samples.  In this example, we assign NA12878 as an affected sample to simulate family case with disease unaffected parents and a disease affected child. 


There are six mandatory columns in ped file:  FamilyID, SampleID, PaternalSampleID, MaternalSampleID, Sex and Phenotype (affected status). For detail of ped format, please check the page (http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped) 
User can start to run example of VarSelect by following command:

```
varselect.pl annotate -v example1.txt -p example1.ped -m family
```

VarSelect uses timestamp of job submit as  job id (ex: 20170105162548).  Directory VarSelectAnalysisResult_20170105162548/ will be created when VarSelect starts running. Log files and intermediate results will be stored in this directory.

After the process of VarSelect, file example1_varselect.db will be created at current directory. User can use gemini commands to query the database for further analysis. 

Please check which column is in the variants table first.

```
gemini db_info example1_varselect.db |grep '^variants'
```

In family mode, there are five specific columns come with your job id: is_AR_20170105162548 , is_CH_20170105162548, is_DR_20170105162548, is_TH_20170105162548 and  is_XL_20170105162548 which represent that the variant is presented in following genetic inheritance models: autosomal recessive, compound heterozygotes, de novo recessive, two hit and X-linked respectively.

You can run following command to filer out which variant fit the compound heterozygotes  model.

```
gemini query --header
       –q ‘select chrom,start,ref,alt,gene,gts from variants
           where is_CH_20170105162548 = 1’
       example1_varselect.db
```

Column “in_analysis_20170105162548” is a union set of five models. Variants that fit any of the five genetic models can be filtered out by following command:

```
gemini query --header
       –q ‘select chrom,start,ref,alt,is_AR_20170105162548 ,
           is_CH_20170105162548, is_DR_20170105162548, is_TH_20170105162548,
           is_XL_20170105162548 from variants where in_analysis_20170105162548 = 1'
       example1_varselect.db   
```

## Examples 2

## Examples 3

## Examples 4

# Publication & Contact
