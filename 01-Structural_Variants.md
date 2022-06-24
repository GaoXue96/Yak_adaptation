# Structural Variants

## 1 SV detection

### 1-1 Reads alignment

#### 1-1-1 NGM-LR alignment

```shell
ngmlr -t 6 -r reference.fa -q sample.fq -o sample.sam -x ont
```

#### 1-1-2 Pacbio alignment

```shell
ngmlr -t 6 -r reference.fa -q sample.fq -o sample.sam
```

#### 1-1-3 Convert to BAM format

```
samtools view -bS sample.sam | samtools sort > sample.bam
samtools index sample.bam
```

### 1-2 SV calling

#### 1-2-1 NGM-LR SV calling

```shell
cuteSV sample.bam reference.fa sample.cuteSV.vcf./ --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads 6 --sample sample 
```

#### 1-2-2 Pacbio SV calling

```shell
cuteSV sample.bam reference.fa sample.cuteSV.vcf ./ --genotype --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 200 --diff_ratio_merging_DEL 0.5 --threads 6 --sample sample
```

## 2 SV genotyping

#### 2-1 Merge all SVs

```shell
ls *.cuteSV.vcf > vcf.txt
SURVIVOR merge vcf.txt 50 1 1 -1 -1 -1 merged_SURVIVOR_cuteSV_1.vcf
```

#### 2-2 NGM-LR SV genotyping

```shell
cuteSV sample.bam reference.fa sample.cuteSV.gt.vcf ./ --genotype -Ivcf merged_SURVIVOR_cuteSV_1.vcf --sample sample --threads 6
```

#### 2-3 Pacbio SV genotyping

```shell
cuteSV sample.bam reference.fa sample.cuteSV.gt.vcf ./ --genotype -Ivcf merged_SURVIVOR_cuteSV_1.vcf --sample sample --threads 6
```

#### 2-4 Merge all VCFs

```shell
ls *.cuteSV.gt.vcf > cuteSV.gt.vcf.txt
SURVIVOR merge cuteSV.gt.vcf.txt 50 1 1 -1 -1 -1 merged_SURVIVOR_cuteSV.gt.vcf
```

## 3 Fst calculation

```
vcftools --vcf merged_SURVIVOR_cuteSV.gt.filter.recode.vcf --weir-fst-pop cattle.txt --weir-fst-pop yak.txt --out cattle.vs.yak.fst
```