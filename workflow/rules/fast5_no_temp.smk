rule sanitize:
    input:
        f"{config['datasets_dir']}/{config['sample']}/{{raw_set}}_raw/",
    output:
        temp(directory(f"{config['datasets_dir']}/{config['sample']}/{{raw_set}}_sanitized/")),
        temp(directory(f"{config['datasets_dir']}/{config['sample']}/{{raw_set}}_merged/"))
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+"
    conda:
        "../envs/fast5api.yaml"
    threads: config['threads']
    params:
        any_fast5 = lambda wildcards, input: find_any_file(wildcards, filetype=".fast5", d=f"{input}")
    shell:
        ("""
        mkdir -p {output[0]}
        mkdir -p {output[1]}
        x=$(python workflow/scripts/match_single_fast5.py "{params.any_fast5}")
        if [ "$x" = "True" ]; then
            single_to_multi_fast5 --input_path {input} --save_path {output[1]} --threads {threads} --recursive
            compress_fast5 -t {threads} --recursive -i {output[1]} -s {output[0]} -c vbz --sanitize
        else
            compress_fast5 -t {threads} --recursive -i {input} -s {output[0]} -c vbz --sanitize
        fi
        """)

rule compress_to_gzip:
    input:
        f"results/{config['sample']}/{{raw_set}}_guppy_{{model}}/workspace/",
    output:
        directory(f"results/{config['sample']}/{{raw_set}}_guppy_{{model}}/gzip/")
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+",
        model="[A-Za-z0-9]+"
    conda:
        "../envs/fast5api.yaml"
    threads: config['threads']
    shell:
        ("""
        mkdir -p {output}
        if [ $(check_compression --recursive -i {input} | grep -c vbz) -ge 1 ] ; then
            compress_fast5 -t {threads} --recursive -i {input} -s {output} -c gzip
        else
            cp -r {input}/* {output}/
        fi
        """)

rule make_single_fast5s:
    input:
        f"results/{config['sample']}/{{raw_set}}_guppy_canonical/gzip/" # to ensure that guppy already ran successfully (and converted pod5 to fast5 if necessary)
    output:
        temp(directory(f"results/{config['sample']}/{{raw_set}}_guppy_canonical/single_fast5s/"))
    conda:
        "../envs/fast5api.yaml"
    params:
        any_fast5 = lambda wildcards, input: find_any_file(wildcards, filetype=".fast5", d=f"{input}")
    threads:
        config['threads']
    shell:
        ("""
        x=$(python workflow/scripts/match_single_fast5.py "{params.any_fast5}")
        if [ "$x" = "True" ]; then
            cp -r {input} {output}
        else
            multi_to_single_fast5 --input_path {input} --save_path {output} --threads {threads} --recursive
        fi
        """)

rule merge_single_fast5s:
    input:
        f"results/{config['sample']}/{{raw_set}}_guppy_canonical/single_fast5s/",
        f"results/{config['sample']}/{{raw_set}}_guppy_canonical/resquiggle.runtime.txt"
    output:
        directory(f"results/{config['sample']}/{{raw_set}}_guppy_canonical/merged_single_fast5s/")
    conda:
        "../envs/fast5api.yaml"
    threads:
        config['threads']
    shell:
        ("""
        single_to_multi_fast5 --input_path {input[0]} --save_path {output} --threads {threads} --recursive
        """)