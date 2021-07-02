##################################################################################
###     Flye and Raven assemblies, polished with Racon x4 and Medaka x 1 
###     
###
###     ORF and Breseq QC metrics are produced 
###
###     
###         
###     Georgia Breckell  11.06.2021
###
###
###
###################################################################################

configfile:
    "configs/assemblies.yml"

rule all:
    input:
        expand("results/{strain}/{assembler}/polished_genome.fasta", strain=config["strain"], assembler=config["assembler"]),
        expand("results/medaka1/{strain}_{assembler}_pilon/pilon.fasta", strain=config["strain"], assembler=config["assembler"]),
        expand("results/medaka1/{strain}_{assembler}_ORF/ORF.pdf", strain=config["strain"], assembler=config["assembler"]),
        expand("results/medaka1/{strain}_{assembler}_breseq/output/output.gd", strain=config["strain"], assembler=config["assembler"]),

#ruleorder: copy_unicycler > Racon_illumina

#Randomly subsample 1000 Nanopore long reads. Generates a new fastq containing the subsampled reads.
rule subsample_ONT:
    input:
        original="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz"
    output:
        "/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_withheld.fastq.gz"
    shell:
        "seqtk sample -s 11 {input} 1000 > {output}"

#Identify which long reads were subsampled by comparing with main dataset. Generates a list of read IDs for the subsampled reads.
rule ID_subsampled_reads_ONT: 
    input:
        subsampled="/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_withheld.fastq.gz",
        original="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz"
    output:
        "/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_common.list"
    shell:
        "seqkit common {input.subsampled} {input.original} | grep '@' | cut -c 2-37 > {output} "
        
#Removes subsampled reads from the main dataset, according to the read IDs on the list above.
rule remove_subsampled_reads_ONT:
    input: 
        common_list="/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_common.list",
        original="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz"
    output:
        "/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_assembly.fastq.gz"
    shell:
        "seqkit grep -f {input.common_list} -v {input.original} -o {output}"

#Randomly subsample 5000 paired Illumina reads, Generates a new fastq containing the subsampled reads
rule Subsample_Illumina:
    input:
        R1="/data/Georgia/SC12A_reads/illumina/{strain}_R1.fastq.gz",
        R2="/data/Georgia/SC12A_reads/illumina/{strain}_R2.fastq.gz"
    output:
        subsampled_R1="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R1_WH.fastq.gz",
        subsampled_R2="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R2_WH.fastq.gz"
    run:
        shell("seqtk sample -s 11 {input.R1} 5000 > {output.subsampled_R1}")
        shell("seqtk sample -s 11 {input.R2} 5000 > {output.subsampled_R2}")

#Identify which long reads were subsampled by comparing with main dataset. Generates a list of read IDs for the subsampled reads.
rule ID_subsampled_reads_Illumina: 
    input:
        subsampled_R1="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R1_WH.fastq.gz",
        R1="/data/Georgia/SC12A_reads/illumina/{strain}_R1.fastq.gz",
        subsampled_R2="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R2_WH.fastq.gz",
        R2="/data/Georgia/SC12A_reads/illumina/{strain}_R2.fastq.gz"
    output:
        R1="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R1_common.list",
        R2="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R2_common.list"
    run:
        shell("seqkit common {input.subsampled_R1} {input.R1} | grep '@' | cut -c 2-37 > {output.R1}") 
        shell("seqkit common {input.subsampled_R2} {input.R2} | grep '@' | cut -c 2-37 > {output.R2}") 
 
##Removes subsampled reads from the main dataset, according to the read IDs on the list above.
#Generates an "Assembly" dataset from which all assemblies are created.    
rule remove_subsampled_reads_Illumina:
    input: 
        R1_common="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R1_common.list",
        R1="/data/Georgia/SC12A_reads/illumina/{strain}_R1.fastq.gz",
        R2_common="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R2_common.list",
        R2="/data/Georgia/SC12A_reads/illumina/{strain}_R2.fastq.gz"
    output:
        R1="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R1_assembly.fastq.gz",
        R2="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R2_assembly.fastq.gz"
    run:
        shell("seqkit grep -f {input.R1_common} -v {input.R1} -o {output.R1}")
        shell("seqkit grep -f {input.R2_common} -v {input.R2} -o {output.R2}")


#Raven long read only assembly with the "Assembly dataset of reads", no renmaing needed as output is in desired convention. 
#Raven output as single fasta file so no movement out of the "output_dir" is needed either 
rule raven_assembly:
    input:
        "/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        assembly="results/{strain}/raven/assembly.fasta",
        dir_removed=touch("results/{strain}/raven/output_dir_removed")
    log:
        "results/{strain}/logs/raven.log"    
    benchmark:
        "results/{strain}/benchmarks/raven.assembly.benchmark.txt"
    run:
       shell("raven {input} > {output.assembly} 2> {log}")

#flye long read only assembly with the "Assembly dataset of reads", no renaming is needed because flye output is in desired convention 
rule flye_assembly:
    input:
        "/data/Georgia/SC12A_reads/ONT_subsampled/{strain}_assembly.fastq.gz"
    output:
        fasta="results/{strain}/flye/flye_output/assembly.fasta",
        gfa="results/{strain}/flye/flye_output/assembly_graph.gfa"
    params:
        out_prefix="results/{strain}/flye/flye_output/"
    log:
        "results/{strain}/logs/flye.log"
    benchmark:
        "results/{strain}/benchmarks/flye.assembly.benchmark.txt"
    conda:
        "environments/assemblies_2_7.yml"
    shell:
        "flye --nano-raw {input} --out-dir {params.out_prefix} --plasmids 2> {log}" 
        

#Copies the flye assembly out of the flye output dir
rule copy_flye:
    input:
        assembly="results/{strain}/flye/flye_output/assembly.fasta",
        gfa="results/{strain}/flye/flye_output/assembly_graph.gfa"
    output:
        assembly="results/{strain}/flye/assembly.fasta",
        gfa="results/{strain}/flye/assembly_graph.gfa"
    run:
        shell("cp {input.assembly} {output.assembly}"),
        shell("cp {input.gfa} {output.gfa}")

#rename GFA file to standard convention
rule rename_flye:
    input:
        gfa="results/{strain}/flye/assembly_graph.gfa"
    output:
        "results/{strain}/flye/assembly.gfa"
    run:
        shell("mv {input.gfa} {output}")

#rule to remove flye assembly directory   
rule remove_flye_output_dir:
    input:
        assembly="results/{strain}/flye/assembly.fasta",
    output:
        touch("results/{strain}/flye/output_dir_removed")
    params:
        output_dir="results/{strain}/flye/flye_output/"
    shell:
        "rm -r {params.output_dir}"

###GENOME POLISHING

#First alignment for the first round of polishing. Assembly reads are aligned to the assembly using Minimap2
#This needs to be a seperate input for each assembler as they all have different 
#output formats, following this first round of polishing the different assemblers can be treated the same 

rule Minimap2_round1:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/assembly.fasta",
    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round1.sam"),
    conda:
        "environments/base.yml"
    shell:
        "minimap2 -a {input.assembly} {input.reads} > {output.alignment}"

#First round of Racon polishing using Minimap2 alignments, assembly reads and the assembly
rule Racon_round1:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz", 
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round1.sam",
        assembly="results/{strain}/{assembler}/assembly.fasta",
    output:
        Racon_round1="results/{strain}/{assembler}/polishing/Racon_round1.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round1}"),

#Second alignment for polishing, the output from the first round of polishing is used as the input assembly. 
rule Minimap2_round2:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/polishing/Racon_round1.fasta",

    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round2.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Second round of Racon polishing. The output from the first round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to first racon output.
rule Racon_round2:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round2.sam",
        assembly="results/{strain}/{assembler}/polishing/Racon_round1.fasta",
    output:
        Racon_round2="results/{strain}/{assembler}/polishing/Racon_round2.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round2}")

#Third alignment for polishing, the output from the second round of polishing is used as the input assembly.
rule Minimap2_round3:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/polishing/Racon_round2.fasta",
    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round3.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Third round of Racon polishing. The output from the second round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to second racon output.
rule Racon_round3:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round3.sam",
        assembly="results/{strain}/{assembler}/polishing/Racon_round2.fasta",
    output:
        Racon_round3="results/{strain}/{assembler}/polishing/Racon_round3.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round3}"),

#Fouth alignment for polishing, the output from the third round of polishing is used as the input assembly.
rule Minimap2_round4:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        assembly="results/{strain}/{assembler}/polishing/Racon_round3.fasta",
    output:
        alignment=temp("results/{strain}/{assembler}/polishing/Minimap2_round4.sam"),
    run:
        shell("minimap2 -a {input.assembly} {input.reads} > {output.alignment}"),

#Fouth round of Racon polishing. The output from the third round is used as the input assembly along with the original reads 
#and the alignment file of the raw reads to third racon output.
rule Racon_round4:
    input:
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz",
        alignment="results/{strain}/{assembler}/polishing/Minimap2_round4.sam",
        assembly="results/{strain}/{assembler}/polishing/Racon_round3.fasta",
    output:
        Racon_round4="results/{strain}/{assembler}/polishing/Racon_round4.fasta",
    run:
        shell("racon -m 8 -x -6 -g -8 -w 500 {input.reads} {input.alignment} {input.assembly} > {output.Racon_round4}"),

#Medaka polishing with long reads

rule Medaka: 
    input: 
        assembly="results/{strain}/{assembler}/polishing/Racon_round4.fasta",
        reads="/data/Georgia/SC12A_reads/ONT/{strain}.fastq.gz"
    params:
        model="r941_min_high_g360",
        outdir="results/medaka1/{strain}_{assembler}_medaka/"
    output:
        touch("results/medaka1/{strain}_{assembler}_medaka_done")
    shell:
        "medaka_consensus -i {input.reads} -d {input.assembly} -o {params.outdir} -m {params.model}"


#BWA indexing of medaka output.  
rule BWA_index_Racon_polished: 
    input:
        marker_file="results/medaka1/{strain}_{assembler}_medaka_done",
    params:
        assembly="results/medaka1/{strain}_{assembler}_medaka/consensus.fasta"
    output:
        touch("results/medaka1/{strain}_{assembler}_medaka_makeidx.done"),
    run:
        shell("bwa index {params.assembly}")

#BWA mem aligns Illumina reads to the medaka output (partially polished genome)
rule BWA_mem_illumina: 
    input:
        R1="/data/Georgia/SC12A_reads/illumina/{strain}_R1.fastq.gz", 
        R2="/data/Georgia/SC12A_reads/illumina/{strain}_R2.fastq.gz", 
        idxdone="results/medaka1/{strain}_{assembler}_medaka_makeidx.done",
    output:
        mapping=temp("results/medaka1/{strain}_{assembler}_Illumina_mapping.sam"),
    params:
        indexed_fasta="results/medaka1/{strain}_{assembler}_medaka/consensus.fasta",
    run:
        shell("bwa mem {params.indexed_fasta} {input.R1} {input.R2} > {output.mapping}"),

#BWA mem alignment is sorted for input into Pilon
rule sorted_BWA_mem_alignment: 
    input:
        mapping="results/medaka1/{strain}_{assembler}_Illumina_mapping.sam",   
    output:
        "results/medaka1/{strain}_{assembler}_Illumina_mapping_sorted.bam",   
    run:
        shell("samtools sort {input.mapping} -o {output}"),

#BWA mem sorted alignment is indexed 
rule indexed_BWA_mem_alignment:
    input:
        "results/medaka1/{strain}_{assembler}_Illumina_mapping_sorted.bam",
    output:
        touch("results/medaka1/{strain}_{assembler}_BWA_alignment_makeidx.done"),
    run:
        shell("samtools index {input}"),

#Pilon polish the genome using Illumina short reads
rule Pilon_polish: 
    input:
        assembly_flag="results/medaka1/{strain}_{assembler}_medaka_done",
        bam="results/medaka1/{strain}_{assembler}_Illumina_mapping_sorted.bam",
        idxdone="results/medaka1/{strain}_{assembler}_BWA_alignment_makeidx.done",
    output:
        pilon="results/medaka1/{strain}_{assembler}_pilon/pilon.fasta",
    params:
        outdir="results/medaka1/{strain}_{assembler}_pilon/",
        assembly="results/medaka1/{strain}_{assembler}_medaka/consensus.fasta"
    run:
        shell("pilon --genome {params.assembly} --bam {input.bam} --outdir {params.outdir}"),

#Making Final Genome stats file
rule genome_stats: 
        input:
            "results/{strain}/{assembler}/polished_genome.fasta"
        output:
            "results/{strain}/{assembler}/genome_stats.txt"
        shell:
            "seqkit stats -T {input} > {output}" 

rule pooled_stats:
        input:
            "results/{strain}/{assembler}/genome_stats.txt"
        output:
            touch("results/{strain}/{assembler}/genome_stats.pooled"),
        params:
            "results/{assembler}_genome_stats.txt"
        shell:
           "cat {input} >> {params}"

####Assembly QC 

###ORF assesment
#Run prodigal step of ORF assesment
rule ORF_prodigal:
        input: 
            "results/medaka1/{strain}_{assembler}_pilon/pilon.fasta"
        output: 
            "results/medaka1/{strain}_{assembler}_ORF/ORF.faa"
        shell: 
            "prodigal -a {output} -q -i {input}"

#Run Diamond step of ORF assesment
rule ORF_diamond:
        input: 
            "results/medaka1/{strain}_{assembler}_ORF/ORF.faa"
        output: 
            "results/medaka1/{strain}_{assembler}_ORF/ORF.data"
        threads: 
            32
        params:
                db="/data/databases/diamond/uniprot.dmnd",
                of="6 qlen slen"
        shell: 
            "diamond blastp --threads 32 --max-target-seqs 1 --db {params.db} --query {input} --outfmt {params.of} --out {output}"

#Run plotting step of ORF assesment
rule ORF_hist:
        input: 
            "results/medaka1/{strain}_{assembler}_ORF/ORF.data"
        output: 
            "results/medaka1/{strain}_{assembler}_ORF/ORF.pdf"
        shell: 
            "R --slave --no-restore --file=./scripts/hist_orf_lengths.R --args {input} {output}"
#Run breseq analysis
rule breseq:
    input:
        read1="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R1_assembly.fastq.gz",
        read2="/data/Georgia/SC12A_reads/illumina_subsampled/{strain}_R2_assembly.fastq.gz",
        assembly="results/medaka1/{strain}_{assembler}_medaka_done"
    output: "results/medaka1/{strain}_{assembler}_breseq/output/output.gd"
    params:
        outdir="results/medaka1/{strain}_{assembler}_breseq",
        assembly="results/medaka1/{strain}_{assembler}_medaka/consensus.fasta"
    shell:
        "breseq -r {params.assembly} -o {params.outdir} {input.read1} {input.read2}"


