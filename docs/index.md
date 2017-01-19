# Introduction

# Install

# Quick start

# Examples

## Examples 1 samples from a family 

In this example, we use variants from chromosome 22 of family NA12878,NA12891,NA12892 as an example for family-based analysis.

All required files are stored in directory varselect/examples/example1/ :

* NA12878_chr22.vcf.gz  : Variant calls from chromosome 22 of these genomes. *

* NA12891_chr22.vcf.gz *

* NA12892_chr22.vcf.gz *

* example1.txt : A comma separated file that describe the relationship between sample and corresponding VCF file.  *

* example1.ped  : A tab separated file that describe the relationship between each sample, the sex of samples and the affected status of samples.  In this example, we assign NA12878 as an affected sample to simulate family case with disease unaffected parents and a disease affected child. *


There are six mandatory columns in ped file:  FamilyID, SampleID, PaternalSampleID, MaternalSampleID, Sex and Phenotype (affected status). For detail of ped format, please check the page ﻿http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped﻿ 
User can start to run example of VarSelect by following command:
varselect.pl annotate -v example1.txt -p example1.ped -m family
VarSelect uses timestamp of job submit as  job id (ex: 20170105162548).  Directory VarSelectAnalysisResult_20170105162548/ will be created when VarSelect starts running. Log files and intermediate results will be stored in this directory.

After the process of VarSelect, file example1_varselect.db will be created at current directory. User can use gemini commands to query the database for further analysis. 

Please check which column is in the variants table first.
gemini db_info example1_varselect.db |grep '^variants'
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
