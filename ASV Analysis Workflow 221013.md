
# ASV Analysis Workflow 221013



## 1. Fastq file pre-processing

### 1.1 fastqc

```
fastqc \
-t 10 \
-o ./FastQC_results \
./WT_CCl4_[1368]wk/*.fastq.gz
```

### 1.2 Trimmomatic
- Trim low quality sequences

```
Installed to /opt/Trimmomatic-0.39
executed via /usr/bin/Trimmomatic
(java jar trimmomatic-0.39.jar)
```

```
#! /usr/bin/Rscript

library(glue)
setwd('/home/hsy/Programs/Trimmomatic-0.39/')

path <- '/data/CHJ_hepatocyte_RNAseq_RAW/'
wks <- c('1','PA1')

for(i in 1:length(wks)){
    wk <- wks[i]
    dir.create(path = file.path('/data/CHJ_hepatocyte_RNAseq_RAW/Trimmed'))
    fastqs <- list.files(path = file.path(path), pattern = 'gz$')
    fastqs <- gsub(pattern='_.\\.fastq\\.gz$', replacement = '',x = fastqs)
    fastqs <- unique(fastqs)
    for(j in 1:length(fastqs)){
        fastq_1 <- glue('{fastqs[j]}_1.fastq.gz')
        fastq_2 <- glue('{fastqs[j]}_2.fastq.gz')
        fastq_1_trim <- glue('{fastqs[j]}_1_trim.fastq.gz')
        fastq_2_trim <- glue('{fastqs[j]}_2_trim.fastq.gz')
        fastq_1_unpaired <- glue('{fastqs[j]}_1_unpaired.fastq.gz')
        fastq_2_unpaired <- glue('{fastqs[j]}_2_unpaired.fastq.gz')
        cmd <- glue('Trimmomatic PE -threads 16 \\
        /data/CHJ_hepatocyte_RNAseq_RAW/{fastq_1} \\
        /data/CHJ_hepatocyte_RNAseq_RAW/{fastq_2} \\
        /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/{fastq_1_trim} \\
        /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/{fastq_1_unpaired} \\
        /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/{fastq_2_trim} \\
        /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/{fastq_2_unpaired} \\
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True TRAILING:3 MINLEN:36')
        system(cmd)
    }
}
```


### 1.2 salmon quantification

#### 1.2.1 Decoy sequence preparation

- Using Conda Environment `salmon`
- Using GRCm39 for mouse genome

```
grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt
```
```
cat gencode.vM29.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome.fa.gz
```
```
salmon index -t gentrome.fa.gz -d decoys.txt -p 16 -i salmon_index --gencode
```


#### 1.2.2 Quantification
```
salmon quant \
-i /data/HSY/Gencode_GRCm39/salmon_index \
-p 16 \
-l A \
--numGibbsSamples 30 \
--gcBias \
--validateMappings \
-o Ctrl_quant \
-1 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/1_1_trim.fastq.gz \
-2 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/1_2_trim.fastq.gz
```
_**salmon quant -o option 에  숫자 와 언더바 같이 들어가면 안됨...**_

```
----incorrect example-----
salmon quant \
-i /data/Gencode_GRCm39/salmon_index \
-p 10 \
-l A \
--numGibbsSamples 30 \
--gcBias \
--validateMappings \
-o WT-V_4-4_quant \
-1 /data/CCl4/Trimmed/WT_CCl4_8wk/WT-V_4-4_1_trim.fastq.gz \
-2 /data/CCl4/Trimmed/WT_CCl4_8wk/WT-V_4-4_2_trim.fastq.gz
```

```
salmon quant \
-i /data/Gencode_GRCm39/salmon_index \
-p 16 \
-l A \
--numGibbsSamples 30 \
--gcBias \
--validateMappings \
-o PA_quant \
-1 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/PA1_1_trim.fastq.gz \
-2 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/PA1_2_trim.fastq.gz
```

##### R script version
```
#! /usr/bin/Rscript

library(glue)
path <- '/data/CCl4/Trimmed'
working_path <- '/data/CCl4/Salmon_Quant'
setwd('/data/CCl4/Salmon_Quant')
wks <- list.files(path)


for(i in 1:length(wks)){
 wk_folder<- file.path(path, wks[i])       
 wk_files <- list.files(wk_folder, pattern = '_trim.fastq.gz')
 for(j in 1:(length(wk_files)/2)){
         sample_1 <- wk_files[(2*j-1)]
         sample_2 <- wk_files[(2*j)]
         sample <- gsub(pattern = '_1_trim.fastq.gz$', replacement = '',x=sample_1)
         print(glue('{wk_folder}: {sample} : {sample_1} : {sample_2}'))
         dir.create(glue('{sample}_quant'))
         setwd(file.path(working_path,glue('{sample}_quant')))
         cmd <- glue('salmon quant \\
         -i /data/Gencode_GRCm39/salmon_index \\
        -p 16 \\
        -l A \\
        --numGibbsSamples 30 \\
        --gcBias \\
        --validateMappings \\
        -o quant \\
        -1 {wk_folder}/{sample_1} \\
        -2 {wk_folder}/{sample_2}')
         print(cmd)
         system(cmd)
         setwd(working_path)
         }
}
```

### 1.3 STAR

#### 1.3.1 STAR 1st pass

```
#! /usr/bin/Rscript

library(glue)
path <- '/data/CCl4/Trimmed'
setwd('/data/CCl4/STAR/STAR_1pass/')
wks <- list.files(path)

for(i in 1:length(wks)){
 wk_folder<- file.path(path, wks[i])       
 wk_files <- list.files(wk_folder, pattern = '_trim.fastq.gz')
 for(j in 1:(length(wk_files)/2)){
         sample_1 <- wk_files[(2*j-1)]
         sample_2 <- wk_files[(2*j)]
         sample <- gsub(pattern = '_1_trim.fastq.gz$',replacement = '',x=sample_1)
         print(glue('{wk_folder}: {sample} : {sample_1} : {sample_2}'))
         cmd <- glue('STAR \\
        --outSAMattributes All \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --readFilesCommand zcat \\
        --runThreadN 16 \\
        --sjdbGTFfile /data/Gencode_GRCm39/gencode.vM29.primary_assembly.annotation.gtf \\
        --outReadsUnmapped Fastx \\
        --outMultimapperOrder Random \\
        --outWigType wiggle \\
        --genomeDir /data/Gencode_GRCm39/STAR/ \\
        --readFilesIn /{wk_folder}/{sample_1} /{wk_folder}/{sample_2} \\
        --outFileNamePrefix {sample}_')
         system(cmd)

         }

}

```

#### 1.3.2 STAR 2nd pass

```
#! /usr/bin/Rscript

library(glue)
path <- '/data/CCl4/Trimmed'
setwd('/data/CCl4/STAR/STAR_2pass/')
wks <- list.files(path)
sjdbfiles <- list.files(path='/data/CCl4/STAR/STAR_1pass/', pattern = "SJ.out.tab$")
sjdbfiles <- paste0('/data/CCl4/STAR/STAR_1pass/',sjdbfiles)
sjdbfiles_str <- paste0(sjdbfiles, collapse = " ")

for(i in 1:length(wks)){
 wk_folder<- file.path(path, wks[i])       
 wk_files <- list.files(wk_folder, pattern = '_trim.fastq.gz')
 for(j in 1:(length(wk_files)/2)){
         sample_1 <- wk_files[(2*j-1)]
         sample_2 <- wk_files[(2*j)]
         sample <- gsub(pattern = '_1_trim.fastq.gz$',replacement = '',x=sample_1)
         print(glue('{wk_folder}: {sample} : {sample_1} : {sample_2}'))
         cmd <- glue('STAR \\
        --outSAMattributes All \\
        --outSAMtype BAM SortedByCoordinate \\
        --quantMode GeneCounts \\
        --readFilesCommand zcat \\
        --runThreadN 16 \\
        --sjdbGTFfile /data/Gencode_GRCm39/gencode.vM29.primary_assembly.annotation.gtf \\
        --sjdbFileChrStartEnd {sjdbfiles_str} \\
        --outReadsUnmapped Fastx \\
        --outMultimapperOrder Random \\
        --outWigType wiggle \\
        --genomeDir /data/Gencode_GRCm39/STAR/ \\
        --readFilesIn {wk_folder}/{sample_1} {wk_folder}/{sample_2} \\
        --outFileNamePrefix {sample}_')
         system(cmd)

         }

}

```

### 1.4 HiSAT2 & StringTie
#### 1.4.1 HiSAT2

HiSAT2 quantification -- now using GRCm38/mm10
```
/home/hsy/Programs/hisat2-2.2.1/hisat2 \
-x grcm38_tran/genome_tran \
-p 16 \
--rna-strandness RF \
-1 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/PA1_1_trim.fastq.gz \
-2 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/PA1_2_trim.fastq.gz \
--dta \
-S PA_Hisat2.sam
```
```
/home/hsy/Programs/hisat2-2.2.1/hisat2 \
-x grcm38_tran/genome_tran \
-p 16 \
--rna-strandness RF \
-1 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/1_1_trim.fastq.gz \
-2 /data/CHJ_hepatocyte_RNAseq_RAW/Trimmed/1_2_trim.fastq.gz \
--dta \
-S Ctrl_Hisat2.sam
--new-summary
````

##### Sorting hisat2 output
```
samtools view -bS PA_Hisat2.sam > PA_Hisat2.bam
samtools sort -@ 16 PA_Hisat2.bam -o PA_Hisat2.sorted.bam

samtools view -@ 16 -bS Ctrl_Hisat2.sam > Ctrl_Hisat2.bam
samtools sort -@ 16 Ctrl_Hisat2.bam -o Ctrl_Hisat2.sorted.bam
```

#### 1.4.2 StringTie

1. Stringtie quantification using original gtf & hisat sorted bam file

```
stringtie \
-p 16 \
-G /data/CHJ_hepatocyte_RNAseq_RAW/stringtie/mm10_NCBI_108.gtf \
-o Ctrl_StringTie.gtf \
--rf  \
/data/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/Ctrl_Hisat2.sorted.bam
```
```
stringtie \
-p 16 \
-G /data/CHJ_hepatocyte_RNAseq_RAW/stringtie/mm10_NCBI_108.gtf \
-o PA_StringTie.gtf \
--rf  \
/data/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/PA_Hisat2.sorted.bam

```
2. Stringtie merge all samples using GTFs
```
stringtie \
--merge \
-p 16 \
-G /data/CHJ_hepatocyte_RNAseq_RAW/stringtie/mm10_NCBI_108.gtf \
-o Merged_StringTie.gtf \
Ctrl_StringTie.gtf \
PA_StringTie.gtf
```
3. Gffcompare - compare merged transcriptome with the reference
```
gffcompare \
-r mm10_NCBI_108.gtf \
Merged_StringTie.gtf
```


4. stringtie quantification using merged GTF
```
stringtie \
-p 16 \
-eB \
-G Merged_StringTie.gtf \
-o ./Ctrl/Ctrl_Merged.gtf \
-A ./Ctrl/Ctrl_Merged.tab \
/data/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/Ctrl_Hisat_mm10.sorted.bam

stringtie \
-p 16 \
-eB \
-G Merged_StringTie.gtf \
-o ./PA/PA_Merged.gtf \
-A PA_Merged.tab \
/data/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/PA_Hisat_mm10.sorted.bam
```


## 2. Downstream Analysis

### 2.1 Gene Level Analysis


### 2.2 Isoform Level Analysis

#### 2.2.1 IsoformSwitchAnalyzeR



#### 2.2.2 DEXSeq

- but can substitute with isoformSwitchAnalyzeR

```
python dexseq_count.py \
-p yes \
-f bam \
-r pos \
/data/Gencode_GRCm39/gencode.vM29.primary_assembly.annotation_DEXseq.gff \
/data/CHJ_hepatocyte_RNAseq_RAW/STAR/STAR_2nd_pass/CtrlAligned.sortedByCoord.out.bam \
Ctrl_dexseq_count.txt


python dexseq_count.py \
-p yes \
-f bam \
-r pos \
/data/Gencode_GRCm39/gencode.vM29.primary_assembly.annotation_DEXseq.gff \
/data/CHJ_hepatocyte_RNAseq_RAW/STAR/STAR_2nd_pass/PAAligned.sortedByCoord.out.bam \
PA_dexseq_count.txt

https://support.bioconductor.org/p/9143537/#9143566
sed 's/\"//g' Ctrl_dexseq_count.txt > Ctrl_dexseq_count_clean.txt
sed 's/\"//g' PA_dexseq_count.txt > PA_dexseq_count_clean.txt
```



## 3. Sequence analysis  

#### Gffread transcript fasta file
```
gffread -O -W -w merged.stringtie.transcripts.exons.all.fa -g /data/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/mm10.fa Merged_StringTie.gtf

gffread -O -W -w merged.macrogen.stringtie.transcripts.exons.all.fa -g /data/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/mm10.fa /data/CHJ_hepatocyte_RNAseq_RAW/result_RNAseq/Novel_transcript_analysis/StringTie/Annotation/merged.gtf

gffread -O -W -w 221115.fixed.merged.stringtie.transcripts.exons.all.fa -g /data/HSY/CHJ_hepatocyte_RNAseq_RAW/HiSAT2/mm10.fa 221115_fixed_merged.gtf

```
