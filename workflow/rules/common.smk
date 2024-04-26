rule copy_version:
    output:
        f"results/{config['sample']}/workflow_version.txt"
    #run:
    #    with open(f"{output}", "w") as f:
    #        print(config['version'], file=f)
    shell:
        ("echo {config[version]} > {output}")

rule extract_meta:
    input:
        f"resources/datasets/{config['sample']}/NAT_raw"
    output:
        f"results/{config['sample']}/metadata.yaml",
        f"results/{config['sample']}/raw_dtype.txt"
    params:
        any_pod5 = lambda wildcards: find_any_file(wildcards, filetype=".pod5"),
        any_fast5 = lambda wildcards: find_any_file(wildcards, filetype=".fast5")
    conda:
        "../envs/meta.yaml"
    script:
        "../scripts/extract_meta.py"

rule choose_pipeline:
    input:
        f"results/{config['sample']}/metadata.yaml",
        "resources/guppy_workflows.csv"
    output:
        f"results/{config['sample']}/config.txt",
        f"results/{config['sample']}/pore.txt"
    priority: 100
    script:
        "../scripts/choose_pipeline.py"

rule find_reference:
    input:
        f"resources/datasets/{config['sample']}/reference/"
    output:
        ref = f"results/{config['sample']}/reference.fasta",
        ref_fp = f"results/{config['sample']}/reference_fp.txt"
    script:
        "../scripts/find_reference.py"

rule bam_to_fastq:
    input:
        f"results/{config['sample']}/{{fn}}.bam"
    output:
        temp(f"results/{config['sample']}/{{fn}}.fastq")
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools fastq -T mv,MM,ML {input} > {output}"

rule samtools_sort:
    input:
        "{fp}.sam"
    output:
        "{fp}.sorted.bam"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort {input} > {output}"

rule samtools_index:
    input:
        "{fp}.bam"
    output:
        "{fp}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input}"

rule starttime:
    input:
        "resources/guppy_workflows.csv",
        "resources/rerio/download_model.py",
    output:
        f"results/{config['sample']}/starttime.txt"
    shell:
        ("""
        mkdir -p results/{config[sample]}
        start=`date +%s`
        echo $start > {output}
        """)

rule starttime:
    input:
        "resources/guppy_workflows.csv",
        "resources/all_models_downloaded.txt",
    output:
        f"results/{config['sample']}/starttime.txt"
    shell:
        ("""
        mkdir -p results/{config[sample]}
        start=`date +%s`
        echo $start > {output}
        """)
