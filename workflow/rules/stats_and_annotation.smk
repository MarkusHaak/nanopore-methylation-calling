rule prokka:
    input:
        f"results/{config['sample']}/reference.fasta"
    output:
        outdir = directory(f"results/{config['sample']}/prokka")
    log:
        f"results/{config['sample']}/prokka/prokka.log"
    conda:
        "../envs/prokka.yaml"
    threads: config['threads']
    shell:
        ("""
        mkdir -p {output.outdir}
        prokka --outdir {output.outdir} --force --cpus {threads} {input} >{log} 2>&1
        """)

rule genomecov_from_bam:
    input:
        f"results/{config['sample']}/{{tool_set_model}}.mapped.sorted.bam"
    output:
        f"results/{config['sample']}/{{tool_set_model}}.genomecov"
    log:
        f"results/{config['sample']}/{{tool_set_model}}.log"
    wildcard_constraints:
        tool_set_model=r"\w+"
    params:
        "-d -split" #"-bg"  # optional parameters
    wrapper:
        "v3.7.0/bio/bedtools/genomecov"

rule contigcov:
    input:
        f"results/{config['sample']}/{{tool_set_model}}.genomecov"
    output:
        f"results/{config['sample']}/{{tool_set_model}}.contig_depth.csv"
    wildcard_constraints:
        tool_set_model=r"\w+"
    script:
        "../scripts/contigcov_from_genomecov.py"

rule nanoplot:
    input:
        f"results/{config['sample']}/{{tool_set_model}}.mapped.sorted.bam"
    output:
        outdir = directory(f"results/{config['sample']}/nanoplot/{{tool_set_model}}/")
    log:
        f"results/{config['sample']}/nanoplot/{{tool_set_model}}/nanoplot.log"
    wildcard_constraints:
        tool_set_model=r"\w+"
    conda:
        "../envs/nanoplot.yaml"
    threads: config['threads']
    shell:
        ("""
        mkdir -p {output.outdir}
        NanoPlot -t {threads} --bam {input} -o {output.outdir} > {log} 2>&1
        """)
