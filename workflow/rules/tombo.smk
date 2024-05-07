rule resquiggle:
    input:
        reference = f"results/{config['sample']}/reference.fasta",
        single_fast5_dir = f"results/{config['sample']}/{{raw_set}}_guppy_canonical/single_fast5s/"
    output:
        f"results/{config['sample']}/{{raw_set}}_guppy_canonical/resquiggle.runtime.txt"
    log:
        failed = f"results/{config['sample']}/{{raw_set}}_guppy_canonical/resquiggle_failed_reads.txt",
        #stderr = f"results/{config['sample']}/{{raw_set}}_guppy_canonical/resquiggle.log",
    conda:
        "../envs/tombo.yaml"
    threads: config['threads']
    shell: # continue even if unexpected errors occure
        ("""
        start=`date +%s`
        tombo resquiggle {input.single_fast5_dir} {input.reference} --processes {threads} --failed-reads-filename {log.failed} --ignore-read-locks --overwrite || true
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output}
        """)

rule tombo_de_novo:
    input:
        f"results/{config['sample']}/NAT_guppy_canonical/resquiggle.runtime.txt",
        reference = f"results/{config['sample']}/reference.fasta",
        single_fast5_dir = f"results/{config['sample']}/NAT_guppy_canonical/single_fast5s/"
    output:
        f"results/{config['sample']}/tombo/de_novo.runtime.txt",
        f"results/{config['sample']}/tombo/de_novo.valid_coverage.plus.wig"
    log:
        f"results/{config['sample']}/tombo/de_novo.log"
    conda:
        "../envs/tombo.yaml"
    threads: config['threads']
    resources:
        tombo = 1
    params:
        per_read_stats = "--per-read-statistics-basename  results/{config[sample]}/tombo/de_novo" if config['tombo_per_read_stats'] else ""
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/tombo/
        tombo detect_modifications de_novo --fast5-basedirs {input.single_fast5_dir} \
            --statistics-file-basename results/{config[sample]}/tombo/de_novo {params.per_read_stats} \
            --processes {threads} 2>&1 | tee -a {log}
        tombo text_output browser_files \
            --fast5-basedirs {input.single_fast5_dir} \
            --statistics-filename results/{config[sample]}/tombo/de_novo.tombo.stats \
            --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd dwell \
            --browser-file-basename results/{config[sample]}/tombo/de_novo 2>&1 | tee -a {log}
        tombo text_output signif_sequence_context --statistics-filename results/{config[sample]}/tombo/de_novo.tombo.stats \
            --genome-fasta {input.reference} --sequences-filename results/{config[sample]}/tombo/de_novo.most_signif.fasta 2>&1 | tee -a {log}
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output}
        """)
        # --per-read-statistics-basename  results/{config[sample]}/tombo/de_novo \

rule tombo_alternative_model:
    input:
        f"results/{config['sample']}/NAT_guppy_canonical/resquiggle.runtime.txt",
        reference = f"results/{config['sample']}/reference.fasta",
        single_fast5_dir = f"results/{config['sample']}/NAT_guppy_canonical/single_fast5s/"
    output:
        f"results/{config['sample']}/tombo/alt.runtime.txt",
        f"results/{config['sample']}/tombo/alt.6mA.valid_coverage.plus.wig",
        f"results/{config['sample']}/tombo/alt.5mC.valid_coverage.plus.wig"
    log:
        f"results/{config['sample']}/tombo/alt.log"
    conda:
        "../envs/tombo.yaml"
    threads: config['threads']
    resources:
        tombo = 1
    params:
        per_read_stats = "--per-read-statistics-basename  results/{config[sample]}/tombo/alt" if config['tombo_per_read_stats'] else ""
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/tombo/
        for mod in 6mA 5mC ; do
            tombo detect_modifications alternative_model --alternate-bases $mod \
                --fast5-basedirs {input.single_fast5_dir} \
                --statistics-file-basename results/{config[sample]}/tombo/alt {params.per_read_stats} \
                --processes {threads} 2>&1 | tee -a {log}
            tombo text_output browser_files \
                --fast5-basedirs {input.single_fast5_dir} \
                --statistics-filename "results/{config[sample]}/tombo/alt.${{mod}}.tombo.stats" \
                --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd dwell \
                --browser-file-basename "results/{config[sample]}/tombo/alt.${{mod}}" 2>&1 | tee -a {log}
            tombo text_output signif_sequence_context --statistics-filename "results/{config[sample]}/tombo/alt.${{mod}}.tombo.stats" \
                --genome-fasta {input.reference} --sequences-filename "results/{config[sample]}/tombo/alt.${{mod}}.most_signif.fasta" 2>&1 | tee -a {log}
        done
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output}
        """)
        #--per-read-statistics-basename  results/{config[sample]}/tombo/alt \

rule tombo_model_sample_compare:
    input:
        f"results/{config['sample']}/NAT_guppy_canonical/resquiggle.runtime.txt",
        f"results/{config['sample']}/WGA_guppy_canonical/resquiggle.runtime.txt",
        reference = f"results/{config['sample']}/reference.fasta",
        NAT_single_fast5_dir = f"results/{config['sample']}/NAT_guppy_canonical/single_fast5s/",
        wga_single_fast5_dir = f"results/{config['sample']}/WGA_guppy_canonical/single_fast5s/"
    output:
        f"results/{config['sample']}/tombo/compare.runtime.txt",
        f"results/{config['sample']}/tombo/compare.valid_coverage.plus.wig"
    log:
        f"results/{config['sample']}/tombo/compare.log"
    conda:
        "../envs/tombo.yaml"
    threads: config['threads']
    resources:
        tombo = 1
    params:
        per_read_stats = "--per-read-statistics-basename  results/{config[sample]}/tombo/compare" if config['tombo_per_read_stats'] else ""
    shell:
        ("""
        start=`date +%s`
        mkdir -p results/{config[sample]}/tombo/
        tombo detect_modifications model_sample_compare \
            --fast5-basedirs {input.NAT_single_fast5_dir} \
            --control-fast5-basedirs {input.wga_single_fast5_dir} \
            --statistics-file-basename results/{config[sample]}/tombo/compare {params.per_read_stats} \
            --processes {threads} 2>&1 | tee -a {log}
        tombo text_output browser_files \
            --fast5-basedirs {input.NAT_single_fast5_dir} \
            --control-fast5-basedirs {input.wga_single_fast5_dir} \
            --statistics-filename "results/{config[sample]}/tombo/compare.tombo.stats" \
            --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd dwell difference \
            --browser-file-basename "results/{config[sample]}/tombo/compare" 2>&1 | tee -a {log}
        tombo text_output signif_sequence_context --statistics-filename "results/{config[sample]}/tombo/compare.tombo.stats" \
            --genome-fasta {input.reference} --sequences-filename "results/{config[sample]}/tombo/compare.most_signif.fasta" 2>&1 | tee -a {log}
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output}
        """)
