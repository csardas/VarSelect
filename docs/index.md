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
