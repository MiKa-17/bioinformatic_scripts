library(stats)
library(BiocGenerics)
library(IRanges)
library(GenomeInfoDb)
library(karyoploteR)


seqIds <- c("NC_000001.11",
            "NC_000002.12",
            "NC_000003.12",
            "NC_000004.12",
            "NC_000005.10",
            "NC_000006.12",
            "NC_000007.14",
            "NC_000008.11",
            "NC_000009.12",
            "NC_000010.11",
            "NC_000011.10",
            "NC_000012.12",
            "NC_000013.11",
            "NC_000014.9",
            "NC_000015.10",
            "NC_000016.10",
            "NC_000017.11",
            "NC_000018.10",
            "NC_000019.10",
            "NC_000020.11",
            "NC_000021.9",
            "NC_000022.11",
            "NC_000023.11",
            "NC_000024.10",
            "NC_012920.1")

seqLen <- c(248956422,
            242193529,
            198295559,
            190214555,
            181538259,
            170805979,
            159345973,
            145138636,
            138394717,
            133797422,
            135086622,
            133275309,
            114364328,
            107043718,
            101991189,
            90338345,
            83257441,
            80373285,
            58617616,
            64444167,
            46709983,
            50818468,
            156040895,
            57227415,
            16569
           )

isCircular <- c(FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE,
                FALSE
                )

SeqInfo <- Seqinfo(seqnames=seqIds, seqlengths=seqLen, isCircular=isCircular, genome="human_genome")

KpNoLabels <- plotKaryotype(SeqInfo, labels=NULL);

kpAxis(KpNoLabels, labels=c("chr1"," "," ",
                                "chr2"," "," ",
                                "chr3"," "," ",
                                "chr4"," "," ",
                                "chr5"," "," ",
                                "chr6"," "," ",
                                "chr7"," "," ",
                                "chr8"," "," ",
                                "chr9"," "," ",
                                "chr10"," "," ",
                                "chr11"," "," ",
                                "chr12"," "," ",
                                "chr13"," "," ",
                                "chr14"," "," ",
                                "chr15"," "," ",
                                "chr16"," "," ",
                                "chr17"," "," ",
                                "chr18"," "," ",
                                "chr19"," "," ",
                                "chr20"," "," ",
                                "chr21"," "," ",
                                "chr22"," "," ",
                                "chrX"," "," ",
                                "chrY"," "," ",
                                "mitochondrion"," "," ")); 
kpAddMainTitle(KpNoLabels, main="full_genome_coverage_human_9c10clone_trimmed_ltr_plasmids"); 
kpAddBaseNumbers(KpNoLabels); 


kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/5_mapping/all_reads_human_genome_ltr_plasmids.sorted.bam", col="darkorange");

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/5_mapping/reads_trimmed_with_plasmid/reads_trimmed_with_plasmid_human_genome_ltr_plasmids.sorted.bam", col="darkblue");

#### reads with index adapter 

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/5_mapping/reads_with_adapter_index_1/reads_with_adapter_index_1.sorted.bam", col="darkgreen");

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/5_mapping/reads_with_adapter_index_2/reads_with_adapter_index_2.sorted.bam", col="darkviolet");
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/5_mapping/reads_with_adapter_index_1_complement/reads_with_adapter_index_1_complement.sorted.bam", col="blue");
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/6_in_house_seq/2_vvgt_insert_sites/5_mapping/reads_with_adapter_index_2_complement/reads_with_adapter_index_2_complement.sorted.bam", col="violet");


#### new run with gag primer 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/4_9c10clone_gag_primer/3_mapping/all_reads_9c10clone_gag_human_genome.sorted.bam", col="violet");


### new run with ltr primer 
kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/2_9c10clone_ltr_primer/3_mapping/all_reads_9c10clone_ltr_human_genome.sorted.bam", col="darkviolet");

kpPlotBAMDensity(KpNoLabels, data="/home/michele/git/4_in_house_seq/2_vvgt_insert_sites/2_repetition_run_0823/2_9c10clone_ltr_primer/3_mapping/reads_trimmed_l1000_9c10clone_ltr_human_genome.sorted.bam", col="darkgreen");

