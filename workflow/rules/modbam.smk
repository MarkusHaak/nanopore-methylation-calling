rule modkit_auto_filter:
    input:
        reference = f"results/{config['sample']}/reference.fasta",
        bam = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam",
        ind = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam.bai"
    output:
        bed = f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.auto_filter.bed",
        runtime = f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.auto_filter.runtime.txt"
    log:
        f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.auto_filter.log"
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+"
    conda:
        "../envs/modkit.yaml"
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/modkit
        modkit pileup {input.bam} {output.bed} --log-filepath {log} \
            --ref {input.reference} --only-tabs
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)

rule modkit_fixed_filter:
    input:
        reference = f"results/{config['sample']}/reference.fasta",
        bam = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam",
        ind = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam.bai"
    output:
        bed = f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.fixed_filter.bed",
        runtime = f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.fixed_filter.runtime.txt",
    log:
        f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.fixed_filter.log"
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+"
    conda:
        "../envs/modkit.yaml"
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/modkit
        modkit pileup {input.bam} {output.bed} --log-filepath {log} --ref {input.reference} \
            --filter-threshold C:{config[modkit_filter_threshold]} \
            --filter-threshold A:{config[modkit_filter_threshold]} \
            --only-tabs
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)

rule modkit_no_filter:
    input:
        reference = f"results/{config['sample']}/reference.fasta",
        bam = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam",
        ind = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam.bai"
    output:
        bed = f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.no_filter.bed",
        runtime = f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.no_filter.runtime.txt",
    log:
        f"results/{config['sample']}/modkit/{{raw_set}}_{{basecaller}}.modkit_pileup.no_filter.log"
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+"
    conda:
        "../envs/modkit.yaml"
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/modkit
        modkit pileup {input.bam} {output.bed} --log-filepath {log} --ref {input.reference} \
            --no-filtering --force-allow-implicit \
            --only-tabs
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)

rule modbam2bed:
    input:
        reference = f"results/{config['sample']}/reference.fasta",
        bam = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam",
        ind = f"results/{config['sample']}/{{raw_set}}_{{basecaller}}.mapped.sorted.bam.bai"
    output:
        bed = f"results/{config['sample']}/modbam2bed/{{raw_set}}_{{basecaller}}.modbam2bed.bed",
        runtime = f"results/{config['sample']}/modbam2bed/{{raw_set}}_{{basecaller}}.runtime.txt",
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+"
    conda:
        "../envs/modbam2bed.yaml"
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/modbam2bed
        modbam2bed --combine {input.reference} {input.bam} > {output.bed}
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)
