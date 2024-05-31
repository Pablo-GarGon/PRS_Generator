#!/bin/bash


########### Parameters

# S: Sumstats file path
# V: list of variants
# O: OUTPUT
# D: directorio output
# F: para no aplicar filtro R2 y MAF (TRUE by default)
# E: EADB correction (FALSE by default)
# R: R2 Filter
# M: MAF filter
# G: extract genotypes (TRUE by default)


print_usage() {
  printf "Usage:
  [-s] PATH to Sumstats FILE (REQUIRED)
  [-v] PATH to variants FILE (REQUIRED). INPUT must be a .txt list. Do NOT include the .txt termination in the argument ( as with PLINK :) )
  [-o] output prefix (REQUIRED)
  [-d] output directory (DEFAULT: CURRENT DIRECTORY)
  [-f] turn off R2 and MAF filtering
  [-e] activates EADB correction
  [-r] R2 Filter THRESHOLD (DEFAULT=0.3)
  [-m] MAF filter THRESHOLD (DEFAULT=0.01)
  [-g] Turn off dosage extraction

  "
}


# EXIT if no arguments
if [[ $# -eq 0 ]] ; then
    echo NO ARGUMENTS INTRODUCED !!!
        print_usage
    exit 0
fi


# DEFAULT values
folder="$PWD"/
f_flag='T'
e_flag='F'
g_flag=true
r_filter='0.3'
maf_filter='0.01'


while getopts 's:v:o:d:fr:m:eg' flag; do
  case "${flag}" in
    s) sumstats="${OPTARG}" ;;
    v) extract="${OPTARG}" ;;
    o) output="${OPTARG}" ;;
        d) folder="${OPTARG}"/ ;;


    f) f_flag='F' ;;
    e) e_flag='T' ;;
        g) g_flag=false ;;


    r) r_filter="${OPTARG}" ;;
    m) maf_filter="${OPTARG}" ;;

    *) print_usage
       exit 1 ;;
  esac
done


echo SUMSTATS FILE is: $sumstats
echo VARIANTS FILE is: $extract
echo OUTPUT PREFIX is: $output
echo OUTPUT DIRECTORY: is $folder
echo EADB correction: $e_flag
echo R2 FILTER set to: $r_filter
echo MAF FILTER set to: $maf_filter
echo FILTER BY R2 and MAF: $f_flag
echo dosage extraction: $g_flag



# Extract genotype information:


# Variables:

# 1 - Name of file with the list of variant IDs (don't include .txt termination)
# 2 - Summary Statistics File
# 3 - Output name

# info

# Edit - SET WORKING DIRECTORY FOLDER    <--------
folder=WDIR

# Edit - SET IMPUTED GENOMIC DATA DIRECTORY    <--------
foldergen=/nas/GRACE/PSP_DEGESCO/PSP_2021/plink_files/TopMED/
# Edit - SET IMPUTATION INFO FILE DIRECTORY    <--------
folderinfo=/nas/GRACE/PSP_DEGESCO/PSP_2021/output_imputation/TopMED/

echo $foldergen
echo $folderinfo

zcat ${folderinfo}chr19.info.gz | head -1 > ${folder}Info.txt

for i in {1..22}
do
zcat ${folderinfo}chr${i}.info.gz | grep -wFf ${extract}.txt  >> ${folder}Info.txt
done


# working directory

## Aditional edits <----------
### Edit here and in the plink command below: name of genetic files. In this case, the script is prepared for pgen files (plink2). YOu can adapt the commands for input vcf files if you wish <--------
### Also edit path to the plink2 binary file <--------

/nas/software/plink2 --pfile ${foldergen}chr16_PSP_2021 --from-bp 19266158 --to-bp 19266158 --chr 16 --export A --out ${folder}temp
awk '{print$2} ' ${folder}temp.raw > ${folder}dosages.raw

rm ${folder}temp.raw

for i in {1..22}
do

/nas/software/plink2 --pfile ${foldergen}chr${i}_PSP_2021 --extract ${extract}.txt --export A --out ${folder}temp
cut  -f 7-  ${folder}temp.raw > ${folder}temp1

paste --delimiters=' ' ${folder}dosages.raw ${folder}temp1 | sed 's/ /\t/g' > ${folder}temp

mv ${folder}temp ${folder}dosages.raw

rm ${folder}temp.raw
rm ${folder}temp1

done

fi
## Edit <- Add your pipeline directory here. 
Rscript PATH_TO_PIPELINE_DIRECTORY/PRS_CALCULATION.R $sumstats ${folder}Info.txt ${folder}dosages.raw $output $folder $f_flag $r_filter $maf_filter $e_flag
