# Saccharome WGS: Processamento

[TOC]

## Pré-processamento: Diginorm

```shell
## Entrelaçando os pares
for fqpair1 in $(ls ./dados/preprocessed/*R1.fastq); do

	fqpair2=$(echo ${fqpair1} | sed 's/R1/R2/')
	fqbase=$(basename ${fqpair1} R1.fastq)
	
	interleave-reads.py \
		${fqpair1} \
		${fqpair2} \
		-o ./dados/interleaved/${fqbase}.pe.fq.gz
		# Fim do interleave-reads.py
done

## Diginorm: normalizando pra uma cobertura de 20
## Tentativa 1 - erro 01 abaixo
normalize-by-median.py \
	-p \
	-k 30 \
	-C 20 \
	-N 4 \
	-x 1e9 \
	--savetable ../analises/19-19-02_diginorm_k30/normC20k30.kh \
	../dados/interleaved/*pe.fq.gz
	
## Tentativa 2 - aumentei N (4->6) e K (20->30)
## Rodado de dentro da pasta analises/19-02_diginorm_k30
normalize-by-median.py \
	-p \
	-k 30 \
	-C 20 \
	-N 6 \
	-x 1e9 \
	--savetable ./normC20k30.kh \
	../../dados/interleaved/*pe.fq.gz
```

Erro 01 do diginorm:

![1550665306161](D:\Workspace\projects\001_saccharome_fer\wgs\assembly-megahit\assets\1550665306161.png)

## Pré-processamento: Diginorm (Parâmetros pipeline)

```shell
## Rodado de dentro da pasta analises/19-02_diginorm_k30
memlimitGB=100
memlimitMB=$((memlimitMB*1000))
khmer_k=23
khmer_byte_x="1e9"
khmer_byte_graph_x="$((memlimitGB/${khmer_N}))e9"
khmer_N="4"

## Paired-end
normalize-by-median.py \
	-p \
	-k ${khmer_k} \
	-C 5 \
	-N ${khmer_N} \
	-x ${khmer_byte_graph_x} \
	--savetable ./normC20k20.kh \
	../../dados/interleaved/*pe.fq
	
## Single-end
normalize-by-median.py \
	-k ${khmer_k} \
	-C 20 \
	-N ${khmer_N} \
	-x ${khmer_byte_x} \
	--savetable ./normC20k20.kh \
	--loadtable ./normC20k20.kh \
	../../dados/interleaved/all.se.fq > \
	./normC20k20.se.out.txt \
	2> normC20k20.se.err.txt
	
## Filtrar por abundância
### Poda leituras em k-mers que são pouco abundantes em leituras de alta cobertura. A opção -V é usada para datasets com cobertura variável.
filter-abund.py \
	-V ./normC20k20.kh \
	./*.keep > \
	./filter-abund.out.txt \
	2> ./filter-abund.err.txta
	
### Extração
for i in *.pe.fq.gz.keep.abundfilt; do
	extract-paired-reads.py ${i} > \
	${i}.extract-paired-reads.out.txt \
	2> ${i}.extract-paired-reads.err.txt
done

### Normalizar por 5 - paired
normalize-by-median.py \
	-C 5 \
	-k ${khmer_k} \
	-N ${khmer_N} \
	-x ${khmer_byte_graph_x} \
	-p \
	--savetable ./normC5k20.kh \
	*.pe.fq.keep.abundfilt.pe > \
	normC5k20.pe.out.txt \
	2> normC5k20.pe.err.txt
	
### Normalizar por 5 - single
normalize-by-median.py \
	-C 5 \
	--loadtable ./normC5k20.kh \
	--savetable ./normC5k20.kh \
	*.pe.fq.keep.abundfilt.se > \
	normC5k20.se.out.txt \
	2> normC5k20.se.err.txt
	
# Comprimir e combinar - paired
for i in `ls *.pe.fq.keep.abundfilt.pe.keep`
do 
	bn=`basename ${i} .pe.fq.keep.abundfilt.pe.keep`
	gzip -9c ${i} > ${bn}.pe.kak.fq.gz 
	
done

# Comprimir e combinar - single
for i in `ls *.pe.fq.keep.abundfilt.se.keep`
do 
	bn=`basename ${i} .pe.fq.keep.abundfilt.se.keep`
    gzip -9c ${i} > ${bn}.se.kak.fq.gz
done

for i in `ls *.se.fq.keep.abundfilt.keep`
do
	bn=`basename ${i} .se.fq.keep.abundfilt.keep`
    gzip -9c ${i} >> ${bn}.se.kak.fq.gz
done

```

## Montagem: Megahit

```shell
## Meta-large
megahit \
    --num-cpu-threads 20 \
    --out-dir ./analises/2019-03-24_Megahit_Sensitive \
    -1 ./dados/preprocessed/CRZ2_R1.fastq,./dados/preprocessed/CRZ3_R1.fastq,./dados/preprocessed/CRZ4_R1.fastq,./dados/preprocessed/ORZ2_R1.fastq,./dados/preprocessed/ORZ3_R1.fastq,./dados/preprocessed/ORZ4_R1.fastq \
    -2 ./dados/preprocessed/CRZ2_R2.fastq,./dados/preprocessed/CRZ3_R2.fastq,./dados/preprocessed/CRZ4_R2.fastq,./dados/preprocessed/ORZ2_R2.fastq,./dados/preprocessed/ORZ3_R2.fastq,./dados/preprocessed/ORZ4_R2.fastq \
    --read ./dados/preprocessed/CRZ2_R1_singletons.fastq,./dados/preprocessed/CRZ2_R2_singletons.fastq,./dados/preprocessed/CRZ3_R1_singletons.fastq,./dados/preprocessed/CRZ3_R2_singletons.fastq,./dados/preprocessed/CRZ4_R1_singletons.fastq,./dados/preprocessed/CRZ4_R2_singletons.fastq,./dados/preprocessed/ORZ2_R1_singletons.fastq,./dados/preprocessed/ORZ2_R2_singletons.fastq,./dados/preprocessed/ORZ3_R1_singletons.fastq,./dados/preprocessed/ORZ3_R2_singletons.fastq,./dados/preprocessed/ORZ4_R1_singletons.fastq,./dados/preprocessed/ORZ4_R2_singletons.fastq \
    --presets meta-sensitive \
    --verbose
    
## Default
megahit \
    --num-cpu-threads 10 \
    --out-dir ./analises/18-02-Megahit_NO_META \
    -1 ./dados/preprocessed/CRZ2_R1.fastq,./dados/preprocessed/CRZ3_R1.fastq,./dados/preprocessed/CRZ4_R1.fastq,./dados/preprocessed/ORZ2_R1.fastq,./dados/preprocessed/ORZ3_R1.fastq,./dados/preprocessed/ORZ4_R1.fastq \
    -2 ./dados/preprocessed/CRZ2_R2.fastq,./dados/preprocessed/CRZ3_R2.fastq,./dados/preprocessed/CRZ4_R2.fastq,./dados/preprocessed/ORZ2_R2.fastq,./dados/preprocessed/ORZ3_R2.fastq,./dados/preprocessed/ORZ4_R2.fastq \
    --read ./dados/preprocessed/CRZ2_R1_singletons.fastq,./dados/preprocessed/CRZ2_R2_singletons.fastq,./dados/preprocessed/CRZ3_R1_singletons.fastq,./dados/preprocessed/CRZ3_R2_singletons.fastq,./dados/preprocessed/CRZ4_R1_singletons.fastq,./dados/preprocessed/CRZ4_R2_singletons.fastq,./dados/preprocessed/ORZ2_R1_singletons.fastq,./dados/preprocessed/ORZ2_R2_singletons.fastq,./dados/preprocessed/ORZ3_R1_singletons.fastq,./dados/preprocessed/ORZ3_R2_singletons.fastq,./dados/preprocessed/ORZ4_R1_singletons.fastq,./dados/preprocessed/ORZ4_R2_singletons.fastq \
    --verbose
    
    
## Meta-large: convencional
megahit \
    --num-cpu-threads 30 \
    --out-dir ./analises/2019-04-12_megahit_conventional \
    -1 ./dados/preprocessed/CRZ2_R1.fastq,./dados/preprocessed/CRZ3_R1.fastq,./dados/preprocessed/CRZ4_R1.fastq \
    -2 ./dados/preprocessed/CRZ2_R2.fastq,./dados/preprocessed/CRZ3_R2.fastq,./dados/preprocessed/CRZ4_R2.fastq \
    --read ./dados/preprocessed/CRZ2_R1_singletons.fastq,./dados/preprocessed/CRZ2_R2_singletons.fastq,./dados/preprocessed/CRZ3_R1_singletons.fastq,./dados/preprocessed/CRZ3_R2_singletons.fastq,./dados/preprocessed/CRZ4_R1_singletons.fastq,./dados/preprocessed/CRZ4_R2_singletons.fastq \
    --verbose
```

Start: 15-02 8:57

## Montagem: SPADES comum

```shell
## Nova rodada: melhor uso de memória?? ILLUSION!
spades-311.py \
	--pe1-1 ./dados/preprocessed/CRZ2_R1.fastq \
	--pe2-1 ./dados/preprocessed/CRZ3_R1.fastq \
	--pe3-1 ./dados/preprocessed/CRZ4_R1.fastq \
	--pe1-2 ./dados/preprocessed/CRZ2_R2.fastq \
	--pe2-2 ./dados/preprocessed/CRZ3_R2.fastq \
	--pe3-2 ./dados/preprocessed/CRZ4_R2.fastq \
	--s1 ./dados/preprocessed/CRZ2_R1_singletons.fastq \
	--s1 ./dados/preprocessed/CRZ2_R2_singletons.fastq \
	--s2 ./dados/preprocessed/CRZ3_R1_singletons.fastq \
	--s2 ./dados/preprocessed/CRZ3_R2_singletons.fastq \
	--s3 ./dados/preprocessed/CRZ4_R1_singletons.fastq \
	--s3 ./dados/preprocessed/CRZ4_R2_singletons.fastq \
	--threads 25 \
	--memory 950 \
	--careful \
	--cov-cutoff 10 \
	-o ./analises/20-02_assembly_spades_CRZ
```

## Montagem: MIRA

**Adiada/cancelada**

Arquivo de manifesto:

```shell
## Manifeso: MIRA
project = saccharome
job = genome,denovo,accurate

## Sequências
readgroup = IlluminaPares
autopairing
data = ./dados/preprocessed/CRZ2_R1.fastq ./dados/preprocessed/CRZ2_R2.fastq ./dados/preprocessed/CRZ3_R1.fastq ./dados/preprocessed/CRZ3_R2.fastq ./dados/preprocessed/CRZ4_R1.fastq ./dados/preprocessed/CRZ4_R2.fastq ./dados/preprocessed/ORZ2_R1.fastq ./dados/preprocessed/ORZ2_R2.fastq ./dados/preprocessed/ORZ3_R1.fastq ./dados/preprocessed/ORZ3_R2.fastq ./dados/preprocessed/ORZ4_R1.fastq ./dados/preprocessed/ORZ4_R2.fastq
technology = solexa

## Singletons
readgroup = IlluminaUnpaired
data =./dados/preprocessed/CRZ2_R1_singletons.fastq ./dados/preprocessed/ORZ2_R1_singletons.fastq
./dados/preprocessed/CRZ2_R2_singletons.fastq ./dados/preprocessed/ORZ2_R2_singletons.fastq
./dados/preprocessed/CRZ3_R1_singletons.fastq ./dados/preprocessed/ORZ3_R1_singletons.fastq
./dados/preprocessed/CRZ3_R2_singletons.fastq ./dados/preprocessed/ORZ3_R2_singletons.fastq
./dados/preprocessed/CRZ4_R1_singletons.fastq ./dados/preprocessed/ORZ4_R1_singletons.fastq
./dados/preprocessed/CRZ4_R2_singletons.fastq ./dados/preprocessed/ORZ4_R2_singletons.fastq
technology = solexa

## Parâmetros
parameters = COMMON_SETTINGS -GENERAL number_of_threads=15 \
		-CLIPPING kmerjunk_completekill=off
```

Linha de comando:

```shell
mira \
	--cwd=/work/rcsilva/projects/saccharome-wgs/analises/18-02_Assembly_MIRA \
	--threads=20 \
	/work/rcsilva/projects/saccharome-wgs/analises/18-02_Assembly_MIRA/manifest.conf
```

## Montagem: MetaSPADES

```shell
spades-311.py \
	--meta \
	-1 ./dados/metaspades/Saccharome_R1.fastq \
	-2 ./dados/metaspades/Saccharome_R2.fastq \
	-o ./analises/18_02_Assembly_metaspades \
	--threads 15 \
	--memory 350
	
## Para o diginorm
spades-311.py \
	--meta \
	--12 ./analises/2019-02-21_kalamazoo_diginorm/saccharome-interleaved.kak.fq \
	-o ./analises/2019-03-27_metaspades-diginorm \
	--threads 15 \
	--memory 700
	
## Montando separado: CRZ
spades-311.py \
	--meta \
	-1 ./dados/preprocessed/concat_per_treatment/CRZ_R1.fastq \
	-2 ./dados/preprocessed/concat_per_treatment/CRZ_R2.fastq \
	-o ./analises/2019-04-01_metaspades_crz \
	--threads 15 \
	--memory 800
	
## Para o diginorm
spades-311.py \
	--meta \
	--12 ./analises/2019-02-21_kalamazoo_diginorm/CRZ.kak.int.fq \
	-o ./analises/2019-04-01_metaspades_crz \
	--threads 15 \
	--memory 800
	
## Montando separado: CRZ2
spades-311.py \
	--meta \
	-1 ./dados/preprocessed/CRZ2_1.fastq \
	-2 ./dados/preprocessed/CRZ2_2.fastq \
	-o ./analises/2019-04-28_metaspades-CRZ2 \
	--threads 25 \
	--memory 500
```

## Clusterização: CD-Hit

Objetivo de puxar centroides/representativos das proteínas do Prokka

```shell
# Montagem: primeira do Megahit
# Proteínas: Prokka (prodigal)
cd-hit-est \
	-i ./analises/18-02_prokka01/testasm.ffn \
	-o ./analises/2019-02-28_cdhit-megahit/megahit.cdhit.ffn \
	-c 0.95 \
	-n 3 \
	-l 10 \
	-aS 0.9 \
	-d 0 \
	-B 0 \
	-p 1 \
	-g 1 \
	-T 20 \
	-M 100000 \
	> ./analises/2019-02-28_cdhit-megahit/cd-hit-est.out.txt \
	2> ./analises/2019-02-28_cdhit-megahit/cd-hit-est.err.txt
```

## Predição: Prodigal

```shell
prodigal \
	-i ./analises/13-02_Megahit_Assembly01/final.contigs.fa \
	-o ./analises/18-02-PRODIGAL01/proteins.gbk
```

## Predição: MetaGeneMark

```shell
Usage:
        gmhmmp [options] [sequence file]
        -r Use RBS for gene start prediction
```

```shell
## Versão suja: não roda a limpa!
gmhmmp -r -o ./analises/2019-03-20_metagenemark/MGMout.lst -m /home/rcsilva/src/MetaGeneMark_linux_64/mgm/MetaGeneMark_v1.mod ./analises/2019-02-13_Megahit_Assembly01/final.contigs.fa
```

## Anotação: Prokka

```shell
## Anotação da primeira montagem do Megahit
prokka \
	./analises/13-02_Megahit_Assembly01/final.contigs.fa \
	--outdir ./analises/18-02_prokka01 \
	--prefix testasm \
	--metagenome \
	--cpus 15
```

## Anotação: USEARCH (q: megahit, s: M5NR-AAs) 

```shell
## 1. Genoma completo
infile="./analises/18-02_prokka01/testasm.faa"
udb_ref='./dados/databases/db-M5NR/udb/M5NR.proteins.usearch.udb'
out_dir='./analises/2019-02-24_annotation_megahit_X_M5NR'

usearch81 \
	-usearch_local ${infile} \
	-db ${udb_ref} \
	-evalue 1e-6 \
	-accel 0.8 \
	-id 0.6 \
	-query_cov 0.5 \
	-blast6out ${out_dir}/annotation.m8 \
	-top_hit_only \
	-maxrejects 2500000 \
	-threads 30
```

## Quantificação: METAKallisto

```shell
preproc_dir='/work/rcsilva/projects/saccharome-wgs/dados/preprocessed'

## Construção do índice
kallisto index \
	-i ./analises/2019-02-26_kallisto_count_megahit01/megahit \
	./analises/13-02_Megahit_Assembly01/final.contigs.fa
	
## Quantificação com as reads
kallisto quant \
	-i ./analises/2019-02-26_kallisto_count_megahit01/megahit \
	-o ./analises/2019-02-26_kallisto_count_megahit01/count \
	--bias \
	-t 25 \
	${preproc_dir}/CRZ2_R1.fastq ${preproc_dir}/CRZ2_R2.fastq \
	${preproc_dir}/CRZ3_R1.fastq ${preproc_dir}/CRZ3_R2.fastq \
	${preproc_dir}/CRZ4_R1.fastq ${preproc_dir}/CRZ4_R2.fastq \
	${preproc_dir}/ORZ2_R1.fastq ${preproc_dir}/ORZ2_R2.fastq \
	${preproc_dir}/ORZ3_R1.fastq ${preproc_dir}/ORZ3_R2.fastq \
	${preproc_dir}/ORZ4_R1.fastq ${preproc_dir}/ORZ4_R2.fastq
	
	
## Loop for
for r1file in $(ls ${preproc_dir}/*_R1.fastq); do

	r2file=$(echo ${r1file} | sed 's/R1/R2/')
	fqbase=$(basename ${r1file} .fastq)
	
	kallisto quant \
		-i ./analises/2019-02-26_kallisto_count_megahit01/megahit \
		-o ./analises/2019-02-26_kallisto_count_megahit01/count/${fqbase} \
		--bias \
		-t 25 \
		${r1file} ${r2file}
		
done
```

## Anotação/API: teste recuperando ID do M5NR

```shell
# ID de teste: 52788eb94121eb62274b2de1d9f5dd48 c907e004d5c5216486926cf97c5d1a7f
curl -X POST -d '{"data":["6587e425296e8c89d8d4c99bf0576ef2"]}' "http://api.metagenomics.anl.gov/m5nr/md5"

## Consultando a API:
http://api.mg-rast.org/m5nr/md5/6587e425296e8c89d8d4c99bf0576ef2
```

## Anotação: Puxando proteínas relevantes com pullseq

```shell
# Executado da pasta 2019-03-07_pull_prokka-megahit_proteins
pullseq \
	--input ./proteins.faa \
	--names ./headers.txt > ./pullseq-out/filtered-proteins.faa
	
## Remover velhas anotações
sed 's/\(^>PROKKA_[0-9]*\).*/\1/' filtered-proteins.faa > to-annotate.faa
```

## Anotação: InterProScan das proteínas DA

```shell
## Lembre-se de remover os asteriscos!
## Também é necessário que a Python do $PATH seja a Python3, ou o programa
## vai dar pau e pifar como uma grande bola de fogo!
/home/bioinfo_softwares/interproscan-5.21-60.0/interproscan.sh \
    --seqtype p \
    --pathways \
    --goterms \
    --iprlookup \
    --output-dir ./analises/2019-03-07_interproscan_megahit \
    --formats tsv,gff3,html,svg \
    --input ./analises/2019-02-18_prokka/proteome.150.clean.faa
    
# instalação: https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload
# deu pau na versão 5.21, mobidb error

## De todo o proteoma
## Lembre-se de remover os asteriscos!
/home/bioinfo_softwares/interproscan-5.21-60.0/interproscan.sh \
    --seqtype p \
    --pathways \
    --goterms \
    --iprlookup \
    --minsize 150 \
    --highmem \
    --output-dir ./analises/2019-03-19_iprscan_all \
    --formats tsv,gff3,html,svg \
    --tempdirname ./temp \
    --input ./analises/2019-02-18_prokka/proteome.clean.faa
    
# instalação: https://github.com/ebi-pf-team/interproscan/wiki/HowToDownload
# deu pau na versão 5.21, mobidb error
```

## Anotação baseada em ortólogos: eggNOG-mapper

```shell
outdir="/work/rcsilva/projects/saccharome-wgs/analises/2019-03-15_eggnog-mapper"

## Bacteria
python2 /data/bioinfo/eggnog-mapper/emapper.py \
	--database bact \
	-i ./analises/2019-02-18_prokka/testasm.faa \
	-m diamond \
	--usemem \
	--cpu 15 \
	--target_orthologs one2one \
	--output_dir ./analises/2019-03-17_eggnog_one2one \
	--output bact_orthologs
	
## Eucariotos
python2 /data/bioinfo/eggnog-mapper/emapper.py \
	-i ./analises/2019-02-18_prokka/testasm.faa \
	--database euk \
	--cpu 15 \
	--target_orthologs one2one \
	--output_dir ./analises/2019-03-17_eggnog_one2one \
	--output euk_hmm_orthologs_true

## Mais estringente?
python2 /data/bioinfo/eggnog-mapper/emapper.py \
	--database bact \
	-i ./analises/2019-02-18_prokka/proteome.150.clean.faa \
	-m diamond \
	--usemem \
	--cpu 15 \
	--id 30 \
	--seed_ortholog_evalue 0.00001 \
	--target_orthologs one2one \
	--output_dir ./analises/2019-03-25_eggnog_stringent \
	--output orthologs
	
## Evidência não-eletrônica
python2 /data/bioinfo/eggnog-mapper/emapper.py \
	--database bact \
	-i ./analises/2019-02-18_prokka/testasm.faa \
	-m diamond \
	--usemem \
	--cpu 15 \
	--go_evidence non-electronic \
	--report_orthologs \
	--target_orthologs many2one \
	--output_dir ./analises/2019-04-01_eggnog_not_electronic \
	--output ortho
	
## Eggnog para Metagenemark
python2 /data/bioinfo/eggnog-mapper/emapper.py \
	--database bact \
	-i ./dados/subsets/mgm.pep150.faa \
	-m diamond \
	--usemem \
	--cpu 15 \
	--go_evidence non-electronic \
	--report_orthologs \
	--target_orthologs many2one \
	--output_dir ./analises/2019-04-15_metagenemark_eggnog_annot \
	--output ortho
```

## Anotação: KOBAS (Hammer)

Anotando as proteínas FAA preditas pelo Prodigal:

```shell
## Rodado do Hammer, instalação do kobas toda deprecada/depende de python específica no Thor, que possui várias Python e Perl todas mexidas

run_kobas.py \
	-i testasm.faa \
	-E 1e-5 \
	-N 10 \
	-d K/B \
	-s ko
```

## Anotação taxonômica: Kraken

```shell
## Fastq
for fqfile1 in $(ls ./dados/preprocessed/*_R1.fastq)
do

	fqfile2=$(echo $fqfile1 | sed 's/R1/R2/')
	fqbase=$(basename $fqfile1 .fastq)
    
	kraken2 \
		--db /usr/local/bioinfo/kraken2/DB/  \
		--output ./analises/2019-03-15_kraken \
		--use-names \
		--use-mpa-style \
		--memory-mapping \
		--threads 15 \
		--confidence 0.8 \
		--paired \
		--report ./analises/2019-03-15_kraken \
		${fqfile1} ${fqfile2}	
done

## Proteoma
kraken2 \
	--db /usr/local/bioinfo/kraken2/DB/  \
	--output ./analises/2019-03-18_kraken_proteome \
	--use-names \
	--use-mpa-style \
	--memory-mapping \
	--threads 15 \
	--confidence 0.8 \
	--report ./analises/2019-03-15_kraken \
	./analises/2019-02-18_prokka/testasm.fna
	
## Kraken oficial: todas as amostras
kraken2 \
	--db /usr/local/bioinfo/kraken2/DB/  \
	--memory-mapping \
	--threads 30 \
	--use-names \
	./analises/2019-02-13_Megahit_Assembly01/final.contigs.fa >> ./analises/2019-04-30_kraken_contigs/out-names.txt
```

## Anotação taxonômica: Kaiju

```shell
## Etapa I: Atribuição do NR procariotos + microrganismos eucariotos
# Proteínas

kaiju \
	-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
	-f /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/kaiju_db_nr_euk.fmi \
	-i ./analises/2019-02-13_Megahit_Assembly01/final.contigs.fa \
	-o ./analises/2019-03-12_kaiju-megahit-prodigal/output_assembly/kaiju.out \
	-z 15 \
	-E 0.0001 \
	-p \
	-v

## Classificação das amostras direto do FASTq para o banco NR+Eucariotos
for fqfile1 in $(ls ./dados/preprocessed/*_R1.fastq)
do

	fqfile2=$(echo $fqfile1 | sed 's/R1/R2/')
	fqbase=$(basename $fqfile1 .fastq)
	
	if [[ ! -s ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.out ]]; then
		kaiju \
		-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
		-f /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/kaiju_db_nr_euk.fmi \
		-i ${fqfile1} \
		-j ${fqfile2} \
		-o ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.out \
		-z 20 \
		-E 0.0001 \
		-v
	fi
	
	if [[ ! -s ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.out.summary ]]; then
		kaijuReport \
		-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
		-n /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/names.dmp \
		-i ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.out \
		-r genus \
		-o ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.out.summary
	fi
	
done

## Classificação das amostras direto do FASTq para o banco Progenomes
for fqfile1 in $(ls ./dados/preprocessed/*_R1.fastq)
do

	fqfile2=$(echo $fqfile1 | sed 's/R1/R2/')
	fqbase=$(basename $fqfile1 .fastq)
	
	if [[ ! -s ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.progenomes.out ]]; then
		kaiju \
		-t /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/nodes.dmp \
		-f /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/kaiju_db_progenomes.fmi \
		-i ${fqfile1} \
		-j ${fqfile2} \
		-o ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.progenomes.out \
		-z 20 \
		-E 0.0001 \
		-v
	fi
	
	if [[ ! -s ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.progenomes.out.summary ]]; then
		kaijuReport \
		-t /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/nodes.dmp \
		-n /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/names.dmp \
		-i ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.progenomes.out \
		-r genus \
		-o ./analises/2019-03-12_kaiju-megahit-prodigal/${fqbase}.kaiju.progenomes.out.summary
	fi
	
done
	
## Etapa II: Gera relatório
kaijuReport \
	-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
	-n /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/names.dmp \
	-i ./analises/2019-03-12_kaiju-megahit-prodigal/kaiju.out \
	-r genus \
	-o ./analises/2019-03-12_kaiju-megahit-prodigal/kaiju.out.summary
	
## Adiciona os nomes (deprecated)
for outfile in $(ls ./*.kaiju.out)
do
	addTaxonNames \
		-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
		-n /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/names.dmp \
		-i ${outfile} \
		-o kaiju.names.out
done

## Classificação com o nome todo
for outfile in $(ls ./*.kaiju.out)
do
	baseout=$(basename ${outfile} .kaiju.out)
	
	kaijuReport \
		-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
		-n /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/names.dmp \
		-i ${outfile} \
		-c 4 \
		-r class \
		-l superkingdom,phylum,class,order,family,genus,species \
		-o ./taxa-levels/${baseout}.class.summary
done

## Gera um relatório por taxon
for outfile in $(ls ./*.kaiju.out)
do
	baseout=$(basename ${outfile} .kaiju.out)
	declare -a taxa=("phylum" "class" "order" "family" "genus" "species")
	
	for taxrank in "${taxa[@]}"
	do
	
		kaijuReport \
			-t /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/nodes.dmp \
			-n /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/names.dmp \
			-i ${outfile} \
			-r ${taxrank} \
			-o ./taxout/${baseout}.${taxrank}.summary
	done
done
	
## Parsear os relatórios finais
# Primeiro, extrai todos os taxons
for sumfile in $(ls ./*.summary)
do
	while IFS=$'\t' read -r -a input; do
	taxname=${input[2]}
	echo ${taxname}
	echo ${taxname} >> ./tmptax_kaiju.1c
	done < ${sumfile}
done

## Remove redundâncias
sort tmptax_kaiju.1c | uniq > tmptax_sorted_kaiju.1c

## Greppa os valores por amostra
for sum_file in $(ls ./*.summary)
do
	sum_base=$(basename ${sum_file} .kaiju.out.summary)
	for tax_line in $(cat ./tmptax_sorted_kaiju.1c); do
		tax_value=$(grep "${tax_line}" ${sum_file} | sed 's/.*\t\(.*\)\t.*/\1/')
		echo "${tax_line} ${tax_value}"
	done
done

## Gera arquivo Krona
kaiju2krona \
	-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
	-n /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/names.dmp \
	-i ./CRZ2.kaiju.out \
	-o CRZ2.kaiju.out.krona
	
## Gerando HTMLs do Krona
for kaijufile in $(ls *.kaiju.out); do

	base=$(basename ${kaijufile} .kaiju.out)
	
	kaiju2krona \
	-t /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/nodes.dmp \
	-n /data/db/metagenomics/db-Kaiju/db-NR_Euk_2019-02/names.dmp \
	-i ${kaijufile} \
	-o ./${base}.krona
	
	ktImportText -o ${base}.html ${base}.krona
done

## Teste todos os níveis
kaijuReport \
	-t /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/nodes.dmp \
	-n /data/db/metagenomics/db-Kaiju/db-Progenomes_2019-02/names.dmp \
	-i ./CRZ3.kaiju.out \
	-r phylum \
	-p \
	-o ./test/CRZ3.test.phy.summary
```

## Anotação taxonômica: taxMaps

```shell
## Molde
taxMaps \
	-1 ${fqfile1} \
	-2 ${fqfile2} \
	-c 15 \
	-d /data/db/metagenomics/db-taxmap/refseq_complete_bacarchvir.lcak300.gem \
	-t /data/db/metagenomics/db-taxmap/taxonomy.tbl \
	-o ./analises/2019-03-28_taxmaps
	
## Amostra-teste
taxMaps \
	-1 ./dados/preprocessed/CRZ2_R1.fastq \
	-2 ./dados/preprocessed/CRZ2_R2.fastq \
	-c 15 \
	-d /data/db/metagenomics/db-taxmap/refseq_complete_bacarchvir.lcak300.gem \
	-t /data/db/metagenomics/db-taxmap/taxonomy.tbl \
	-o ./analises/2019-03-28_taxmaps

```



## Anotação taxonômica: Metaphlan

Análise feita para o conjunto de *reads*, não para a montagem

```shell
## Juntando outputs
merge_metaphlan_tables.py \
	CRZ2.out CRZ3.out CRZ4.out ORZ2.out ORZ3.out ORZ4.out > test.txt

## Convertendo em heatmap: species
metaphlan_hclust_heatmap.py \
	--in test.txt \
	--out heatmap \
	-m average \
	-f jaccard
	
## Todos os níveis
python2 /usr/local/bioinfo/Metaphlan2/utils/metaphlan_hclust_heatmap.py \
	--in ./test.txt \
	--out ./heatmap_genus.png \
	--tax_lev 'a' \
	-m average
```

## Anotação e binning: Metameta

Primeiro é necessário configurar o arquivo exemplo `yaml`, e depois..

1.  Rodar o teste de sanidade

```shell
metameta --configfile ./analises/2019-03-25_metameta/myconfig.yaml -np
```

2.  Rodar a pipeline em si

```shell
## Tentatva 1: deu pau
## 3 parametros adicionados top-down
metameta \
	--configfile ./analises/2019-03-25_metameta/myconfig.yaml \
	--use-conda \
	--keep-going \
	--cores 15
```

## Teste de sanidade: busca da montagem VS cana

```shell
usearch81 \
	-ublast analises/2019-02-13_Megahit_Assembly01/final.contigs.fa \
	-db ./dados/refs/sugarcane-genome/cirad.udb \
	-evalue 1e-9 \
	-accel 0.8 \
	-alnout ./analises/2019-03-25_sanity_cane_alignment/out.aln \
	-maxaccepts 1 \
	-maxrejects 0 \
	-strand both \
	-blast6out ./analises/2019-03-25_sanity_cane_alignment/usearch.out.m8 \
    -qsegout ./analises/2019-03-25_sanity_cane_alignment/queries.fa \
	-matched ./analises/2019-03-25_sanity_cane_alignment/matched.fa
```

## Teste de sanidade: busca dos adaptadores

```shell
usearch81 \
	-search_oligodb ./analises/2019-02-15-Cleanse-ORZ4/prontos/ORZ4_R2.fastq \
	-db ./dados/adapters/nextera.fa \
	-strand both \
	-userout ./analises/2019-03-27_sanity_adapters/ORZ4_R2.txt \
	-userfields query+target+qstrand+diffs+tlo+thi+trowdots
	-threads 30
```

## Massagem de dados: EGGNOG para matriz KEGG

```shell
## Passos
# 1. Extrair do EGGNOG as colunas com PROKKA_ID e KEGG
# 2. Limpar a matriz do KEGG (somente vazios)
# 3. Parsear as vírgulas

# 4. Deixar apenas a coluna kegg e então montar o dicionário
## Só coluna 02
awk 'BEGIN { FS = "\t" } ; { OFS = "\t" } ; { print $2 }' kegg_headers_correct2.tsv >> kegg_only.1c.txt

## Busca o dicionário (shell script)
KEGG_getter.sh kegg_only.1c.txt keggdict.tsv

### Arquivos
# Feito no R!
```

Modelo de saída

| contig | proteína | KO   | origem | função | rota | nome | contagens |
| ------ | -------- | ---- | ------ | ------ | ---- | ---- | --------- |
|        |          |      |        |        |      |      |           |

## Visualização R - GOs

https://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2671-2

## Teste: busca de oligos (PASS)

```shell
usearch81 \
	-search_oligodb Saccharome_R1.fastq \
	-db /data/db/adapters/all_adapters.fa \
	-strand both \
	-userout ./R1.primers.tsv \
	-userfields query+target+qlo+qhi+qstrand
```

## Análises tentadas e abandonadas

-   Montagem com metaspades: sem memória
-   Montagem com SPADES: trava sem memória
-   Montagem com IDBA-UD: trava na tabela de k-mers
-   Montagem com kalamazoo -> metaspades: sem memória
-   Análise de GOs obtidos do eggnog: pouco sentido (tanto T, quanto P, quanto Q)
-   Kraken: pouca classificação ( < 20%)
-   TaxMaps: demoradíssmo (2 dias / amostra / 20 cpus)

# Citações

1. Atropos:

   1. hold

2. Fastqc:

   1. hold

3. Prinseq:

   1. Hold

4. Kalamazoo:

   1. MR Crusoe et al., 2014. http://dx.doi.org/10.6084/m9.figshare.979190
   2. CT Brown et al., arXiv:1203.4802 [q-bio.GN]

5. Megahit:

   1. hold

6. Spades:

   1. hold

7. USEARCH

    1. EDGAR, Robert C. Search and clustering orders of magnitude faster than BLAST. **Bioinformatics**, v. 26, n. 19, p. 2460-2461, 2010.

8. M5NR / MG-Rast

9. GNU Parallel

10. Kaiju

11. Kraken

12. eggnog-mapper

13. KOBAS

14. STAMP

15. Gene ontology

       
