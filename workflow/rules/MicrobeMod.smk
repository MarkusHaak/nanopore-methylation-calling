rule MicrobeMod_annotate_rm:
    input:
        "resources/MicrobeMod/MicrobeMod/db/", # ensure installation of MicrobeMod
        bam = f"results/{config['sample']}/{{tool_set_model}}.mapped.sorted.bam", # ensure that dorado ran successfully (R10.4.1)
        reference = f"results/{config['sample']}/reference.fasta"
    output:
        f"results/{config['sample']}/MicrobeMod/{{tool_set_model}}.rm.genes.tsv"
    wildcard_constraints:
        tool_set_model=r"\w+"
    conda:
        "../envs/MicrobeMod.yaml"
    threads:
        config['threads']
    shell:
        ("""
        mkdir -p results/{config[sample]}/MicrobeMod
        MicrobeMod annotate_rm -f {input.reference} -o results/{config[sample]}/MicrobeMod/{wildcards.tool_set_model} -t {threads}
        """)

rule MicrobeMod_call_methylation:
    input:
        "resources/MicrobeMod/MicrobeMod/db/", # ensure installation of MicrobeMod
        f"results/{config['sample']}/{{tool_set_model}}.mapped.sorted.bam.bai",
        bam = f"results/{config['sample']}/{{tool_set_model}}.mapped.sorted.bam",
        reference = f"results/{config['sample']}/reference.fasta"
    output:
        f"results/{config['sample']}/MicrobeMod/{{tool_set_model}}.motifs.tsv",
        runtime = f"results/{config['sample']}/MicrobeMod/{{tool_set_model}}.runtime.txt"
    wildcard_constraints:
        tool_set_model=r"\w+"
    conda:
        "../envs/MicrobeMod.yaml"
    threads:
        config['threads']
    shell:
        ("""
        mkdir -p results/{config[sample]}/MicrobeMod
        start=`date +%s`
        MicrobeMod call_methylation -b {input.bam} -r {input.reference} \
            -o results/{config[sample]}/MicrobeMod/{wildcards.tool_set_model} -t {threads}
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)
