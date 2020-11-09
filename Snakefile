configfile: "config.yaml"

rule all:
    input:
        "calls/all.vcf"

############## DOWNLOAD ##############
# Rule that downloads the genomes that will
# be aligned to the reference genome
rule download_genomes:
    params:
        samples = config["SRA_IDs"],
        output_dir = config["storage_dir"] + "data/",
        threads = config["threads"]
    output:
        config["storage_dir"]+"data/{sample}.fastq"
    threads: int(config["threads"])
    shell:
        "parallel-fastq-dump -t {params.threads} -s {params.samples} -O {params.output_dir} 2> {log}"

############## BWA ##############
# Rule that maps the genomes to the ref genome
rule bwa_map:
    input:
        config["refgenome"],
        config["storage_dir"]+"data/{sample}.fastq"
    output:
        config["storage_dir"]+"mapped_reads/{sample}.bam"
    threads: int(config["threads"])
    shell:
        "bwa mem -t 16 {input} | samtools view -Sb - > {output}"

############## SAM-TOOLS ##############
rule samtools_sort:
    input:
       config["storage_dir"]+"mapped_reads/{sample}.bam"
    output:
        config["storage_dir"]+"sorted_reads/{sample}.bam"
    threads: int(config["threads"])
    log:
        "logs/samtools_sort_{sample}.log"
    shell:
        "samtools sort -@ {threads} -T {output} -O bam {input} > {output} 2> {log}"
        
rule samtools_index:
    input:
        config["storage_dir"]+"sorted_reads/{sample}.bam"
    output:
        config["storage_dir"]+"sorted_reads/{sample}.bam.bai"
    log:
        "logs/samtools_index_{sample}.log"
    shell:
        "samtools index  {input} 2> {log}"

############## Variant Calling ##############
rule bcftools_call:
    input:
        fa=config["refgenome"],
        bam=expand(config["storage_dir"]+"sorted_reads/{sample}.bam", sample=config["SRA_IDs"]),
        bai=expand(config["storage_dir"]+"sorted_reads/{sample}.bam.bai", sample=config["SRA_IDs"]),
    output:
        "calls/all.vcf"
    threads: int(config["threads"])
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call --threads {threads} -mv - > {output} "

############## CLEAN ##############
# Remove the output directory and its contents
rule clean:
    shell:
        'rm data/* ; rm output/*'
