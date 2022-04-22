srs_name=""
srr_id=""

mkdir ${srs_name}
cd ${srs_name}
prefetch -O ./ ${srr_id} -X 100000000000

#####convert sra to fastq
fastq-dump --split-3 --gzip *.sra

#####quality control & mapping to genome

if [ -e *_2.fastq.gz ]

then
java -jar ${trimmomatic}/trimmomatic-0.39.jar PE -phred33 *_1.fastq.gz *_2.fastq.gz ${srs_name}_1.clean.fq.gz ${srs_name}_1_unpaired.fastq.gz ${srs_name}_2.clean.fq.gz ${srs_name}_2_unpaired.fastq.gz -threads 5 \
ILLUMINACLIP:${trimmomatic}/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

STAR --runThreadN 5 \
--genomeDir ${Bovine_genome} \
--sjdbGTFfile ${Bovine_genome}/Bos_taurus.ARS-UCD1.2.96.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn ${srs_name}_1.clean.fq.gz ${srs_name}_2.clean.fq.gz \
--outFileNamePrefix ${srs_name}-STAR

else

java -jar ${trimmomatic}/trimmomatic-0.39.jar SE -phred33 *.fastq.gz ${srs_name}.clean.fq.gz -threads 5 \
ILLUMINACLIP:${trimmomatic}/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

STAR --runThreadN 5 \
--genomeDir  ${Bovine_genome} \
--sjdbGTFfile ${Bovine_genome}/Bos_taurus.ARS-UCD1.2.96.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn ${srs_name}.clean.fq.gz \
--outFileNamePrefix ${srs_name}-STAR

fi

######FPKM quantification
stringtie -p 5 -e -B -G ${Bovine_genome}/Bos_taurus.ARS-UCD1.2.96.gtf \
-o ./${srs_name}-STAR_STARgenome/${srs_name}.gtf -A ../Expression/${srs_name}.tsv ${srs_name}-STARAligned.sortedByCoord.out.bam


######Run featureCounts
featureCounts -T 10 -p -t exon -g gene_id -a ${Bovine_genome}/Bos_taurus.ARS-UCD1.2.96.gtf \
-o ../Expression/${srs_name}.featureCounts.txt ${srs_name}-STARAligned.sortedByCoord.out.bam
