rule map_guppy_calls:
    input:
        fastq = f"results/{config['sample']}/{{raw_set}}_guppy_{{model}}/pass/",
        reference = f"results/{config['sample']}/reference.fasta"
    output:
        temp(f"results/{config['sample']}/{{raw_set}}_guppy_{{model}}.mapped.sam")
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+",
        model="[A-Za-z0-9]+"
    conda:
        "../envs/map_guppy.yaml"
    shell:
        ("""
        if compgen -G "{input.fastq}/*.bam" > /dev/null ; then
            for bamfile in {input.fastq}/*.bam ; do 
                samtools fastq -T mv,MM,ML "${{bamfile}}" > "${{bamfile}}.fastq"
            done
            minimap2 --secondary=no -ax map-ont -y {input.reference} {input.fastq}/*.bam.fastq > {output}
        else
            minimap2 --secondary=no -ax map-ont -y {input.reference} {input.fastq}/*.fastq > {output}
        fi
        """)

rule map_dorado_calls:
    input:
        fastq = f"results/{config['sample']}/{{raw_set}}_dorado.fastq",
        reference = f"results/{config['sample']}/reference.fasta"
    output:
        temp(f"results/{config['sample']}/{{raw_set}}_dorado.mapped.sam")
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --secondary=no -ax map-ont -y {input.reference} {input.fastq} > {output}"
