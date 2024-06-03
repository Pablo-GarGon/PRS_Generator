INPUTS:

Run the bash script PRS_generator.sh. This will extract the genetic info and call the R script which will perform the PRS calculations.
This script is programmed to use PLINK2 pgen files to extract genetic information. In case your input data is in a different format, convert to pgen or edit the script.
**Make sure you have a pruned list of independent variants when using this software**

EXAMPLE FILES:

    List (.txt) of variants to use, corresponding with the IDs on your genetic file. Example: GWS_Intelligence_DAVIES_2018_variantsforPRS.txt
    Summary Statistics Table. Example: Summary_statistics_GWS_Intelligence_DAVIES_2018.txt

Columns needed in Summary Statistics file:

    Chromosome: Chromosome (numeric 1,2,...,22)
    Position (hg38)
    MarkerName (rs). NOT mandatory, but useful
    A1_summary: Effect allele; the allele whose effect is reported
    A2_summary: The other allele
    Beta: Effect size. It is important that this is a beta or log(OR). If it is in another unit, it must be transformed.
    SE: Standard deviation. NOT mandatory
    Pvalue: NOT mandatory
    Effect Allele Frequency: NOT mandatory
    ID: SNP ID in your genetic file

To see the options: ./PRS_generator.sh

To test different parameters of MAF and R2, or to play with the parameters in general, use the -g option to deactivate the extraction of genetic data.
