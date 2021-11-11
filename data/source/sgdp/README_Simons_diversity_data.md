#Simons Diversity Project

This collection contains data from the [Simons Genome Diversity Project](https://www.simonsfoundation.org/life-sciences/simons-genome-diversity-project-dataset/).

**Please respect [Fort Lauderdale Principles](https://www.genome.gov/pages/research/wellcomereport0303.pdf) and see the [accompanying Data Reuse policy](https://github.com/igsr/1000Genomes_data_indexes/blob/master/data_collections/simons_diversity_data/README_Simons_diversity_datareuse_statement.md) with regard to this data.**

###About the Simons Genome Diversity Project

Information about the project is available on the [Simons Foundation website](https://www.simonsfoundation.org/life-sciences/simons-genome-diversity-project-dataset/):

"The sampling strategy differs from studies of human genome diversity that have aimed to maximize medical relevance by studying populations with large numbers of present-day people. This new study takes a different approach by sampling populations in a way that represents as much anthropological, linguistic and cultural diversity as possible, and thus includes many deeply divergent human populations that are not well represented in other datasets."

The genomes were sequenced to at least 30x coverage using Illumina technology at Harvard Medical School. Contacts there regarding the data are: Shop Mallick (shop@genetics.med.harvard.edu), Nick Patterson (nickp@broadinstitute.org) or David Reich (reich@genetics.med.harvard.edu).

###Samples in the Simons Diversity Project

Please note that some of the samples included in the Simons Diversity Project are also present in either the [1000 Genomes Project](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/) or the [HGDP](https://github.com/igsr/1000Genomes_data_indexes/blob/master/data_collections/HGDP/README_HGDP.md) collections.

###Questions

If you have any questions please email info@1000genomes.org.

####Programs and reference data
The Simons diversity data was aligned to the reference genome using the following programs and reference datasets:

1. [BWA-MEM](https://github.com/lh3/bwa/blob/master/bwakit/README.md) [bwakit-0.7.15](https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download)
2. [Samtools-1.9](http://www.htslib.org/doc/samtools-1.9.html)
3. [BioBamBam-0.0.191](https://github.com/gt1/biobambam/releases/tag/0.0.191-release-20150401083643)
4. [GATK-3.3-0](https://github.com/broadgsa/gatk-protected/tree/3.3)
5. [Cramtools-3.0](https://github.com/enasequence/cramtools/tree/cram3)
6. Indels for realignment 
   - Phase 3 biallelic indels from Shapeit2 lifted to GRCh38 using the NCBI remapper ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL.wgs.1000G_phase3.GRCh38.ncbi_remapper.20150424.shapeit2_indels.vcf.gz
   - High quality, experiment-validated indel set produced by Devine and Mills. The coordinates were lifted to GRCh38 by Alison Meynert from IGMM in Edinburgh using CrossMap and the UCSC chain files, filtered out ref == alt cases ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/Mills_and_1000G_gold_standard.indels.b38.primary_assembly.vcf.gz
7. SNPs for recalibration 
   - dbSNP 142 in GRCh38 coordinates ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/other_mapping_resources/ALL_20141222.dbSNP142_human_GRCh38.snps.vcf.gz

####Reference genome: GRCh38 with alternative sequences, plus decoys and HLA

The reference genome that the data was aligned to is in this file: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

This was sourced from: ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna

The HLA sequence in the file came from Heng Li's [bwakit distribution](http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download). 

This reference genome contains both a primary assembly and alternative sequences, produced by the [Genome Reference Consortium](http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/), to which decoy and HLA sequences have been added. Further information regarding this is available from [BWA](https://github.com/lh3/bwa/blob/master/README-alt.md).

####Command lines
1. Alignment at run level

          bwa mem  -t 1 -B 4 -O 6 -E 1 -M -R $rg_string $reference_fasta_file $fastq_file(1) $fastq_file(2) | k8 bwa-postalt.js -p  $prefix_hla_hit $reference_fasta_file.alt | samtools view -1 - > $bam_file

2. Local realignment around known indels by GATK

          java $jvm_args -jar GenomeAnalysisTK.jar -T IndelRealigner -R $reference_fasta -I $bam_file -o $realigned_bam_file -targetIntervals $intervals_file -known $known_indels_file(s) -LOD 0.4 -model KNOWNS_ONLY -compress 0 --disable_bam_indexing

3. Recalibrate base quality scores using known SNPs by GATK

          java $jvm_args -jar GenomeAnalysisTK.jar -T BaseRecalibrator -nt 1 -l INFO -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate -R $reference_fasta -o $recal_data.table -I $bam_file -knownSites $known_snps_from_dbSNP142
          java $jvm_args -jar GenomeAnalysisTK.jar -T PrintReads -l INFO -R $reference_fasta -o $recalibrated_bam -I $bam_file -BQSR $recal_data.table --disable_bam_indexing

4. Mark Duplicates using BioBamBam

          bammarkduplicates I=$input_bam O=$output_bam index=1 rmdup=0

5. Mergeing Library level bam files to Sample level bam files using BioBamBam

          bammerge I=$input_bam1 I=$input_bam2 ...... index=1 indexfilename=$output_index_file_name SO=coordinate > $merged_bam

6. Creating lossless CRAM files using Cramtools

          java cramtools-3.0.jar cram --input-bam-file $input_bam --output-cram-file $output_cram --capture-all-tags --ignore-tags OQ:CQ:BQ --preserve-read-names --lossy-quality-score-spec *8 --reference-fasta-file $reference_fasta



###Alignment availability

The alignment files created by this approach are listed in the index file: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/simons_diversity_data/simons_diversity_data.GRCh38DH.alignment.index

The alignments are located beneath this directory, under data. There, the files are origanised by population and then sample. In each sample directory, an alignment directory contains the files.

The alignments have been made available in CRAM format. Additional information on CRAM can be found in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README_using_1000genomes_cram.md. For each data set a .cram file, .crai index and .bas file is provided. Additional information on these file types can be found in ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/README_file_formats_and_descriptions.md.

###Further information

If you require further information, please contact info@1000genomes.org.


