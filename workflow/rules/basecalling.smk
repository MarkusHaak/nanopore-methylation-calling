rule dorado:
    input:
        f"resources/dorado-{config['dorado_version']}/bin/dorado",
        f"{config['datasets_dir']}/{config['sample']}/{{raw_set}}_raw/",
        f"results/{config['sample']}/metadata.yaml",
        "resources/all_models_downloaded.txt", #all_model_paths(),
        f"results/{config['sample']}/pore.txt",
        f"results/{config['sample']}/raw_dtype.txt"
    output:
        bam = temp(f"results/{config['sample']}/{{raw_set}}_dorado.bam"),
        runtime = f"results/{config['sample']}/{{raw_set}}_dorado.runtime.txt"
    #log:
    #    f"results/{config['sample']}/{{raw_set}}_dorado.log" #  2> >(tee -a {log} >&2)
    wildcard_constraints:
        raw_set="[A-Za-z0-9]+"
    conda:
        "../envs/pod5.yaml"
    params:
        raw_dtype = parse_raw_dtype,
        raw_dir = lambda wildcards, input: f"{input[1]}",#get_raw_dir,
        simplex = lambda wildcards: get_dorado_models_from_config(wildcards, "simplex"),
        mods = lambda wildcards: get_dorado_models_from_config(wildcards, "modifications"),
    resources:
        gpu = 1
    shell:
        ("""
        simplex={params.simplex}
        mods={params.mods}
        start=`date +%s`
        if [ ! -z "$simplex" ]; then
            if [ ! -z "$mods" ]; then
                if [[ {params.raw_dtype} == fast5 ]]; then
                    pod5 convert fast5 {params.raw_dir} --recursive --output {params.raw_dir} --one-to-one {params.raw_dir}
                fi
                resources/dorado-{config[dorado_version]}/bin/dorado basecaller {params.simplex} {params.raw_dir} \
                    --recursive \
                    --modified-bases-models {params.mods} \
                    --emit-moves -b {config[dorado_batchsize]} > {output.bam}
            fi
        fi
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)

rule guppy_canonical:
    input:
        "resources/ont-guppy/bin/guppy_basecaller",
        f"{config['datasets_dir']}/{config['sample']}/{{raw_set}}_sanitized/",
        f"results/{config['sample']}/metadata.yaml",
        "resources/all_models_downloaded.txt",#all_model_paths(),
        f"results/{config['sample']}/pore.txt",
        f"results/{config['sample']}/config.txt",
        f"results/{config['sample']}/raw_dtype.txt"
    output:
        temp(directory(f"results/{config['sample']}/{{raw_set}}_guppy_canonical/workspace/")),
        temp(directory(f"results/{config['sample']}/{{raw_set}}_guppy_canonical/pass/")),
        temp(directory(f"results/{config['sample']}/{{raw_set}}_guppy_canonical/fail/")),
        runtime = f"results/{config['sample']}/{{raw_set}}_guppy_canonical/runtime.txt"
    log:
        f"results/{config['sample']}/{{raw_set}}_guppy_canonical/guppy.log"
    conda:
        "../envs/pod5.yaml"
    resources:
        gpu = 1
    priority: 50 # running tombo takes longer than the modkit jobs
    params: 
        raw_dtype = parse_raw_dtype,
        raw_dir = lambda wildcards, input: f"{input[1]}",#get_raw_dir,
        config = parse_config,
        mods = lambda wildcards: get_guppy_models_from_config(wildcards, "modifications"),
    shell:
        ("""
        config={params.config}
        mods={params.mods}
        start=`date +%s`
        if [ ! -z "$config" ]; then
            if [ ! -z "$mods" ]; then
                mkdir -p results/{config[sample]}/{wildcards.raw_set}_guppy_canonical
                if [[ {params.raw_dtype} == pod5 ]]; then
                    if [ ! -d "./{params.raw_dir}/fast5s/" ]; then
                        mkdir -p ./{params.raw_dir}/fast5s/
                        pod5 convert to_fast5 ./{params.raw_dir}/*.fast5 --output ./{params.raw_dir}/fast5s/ 2>&1 | tee -a {log}
                    fi
                    resources/ont-guppy/bin/guppy_basecaller -i {params.raw_dir}/fast5s/ \
                        -s results/{config[sample]}/{wildcards.raw_set}_guppy_canonical/ \
                        -c {params.config} --recursive -x 'cuda:0' --fast5_out 2>&1 | tee -a {log}
                else
                    resources/ont-guppy/bin/guppy_basecaller -i {params.raw_dir} \
                        -s results/{config[sample]}/{wildcards.raw_set}_guppy_canonical/ \
                        -c {params.config} --recursive -x 'cuda:0' --fast5_out 2>&1 | tee -a {log}
                fi
            fi
        fi
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)

rule guppy_modified:
    input:
        "resources/ont-guppy/bin/guppy_basecaller",
        f"{config['datasets_dir']}/{config['sample']}/NAT_sanitized/",
        f"results/{config['sample']}/metadata.yaml",
        "resources/all_models_downloaded.txt", #all_model_paths(),
        f"results/{config['sample']}/pore.txt",
        f"results/{config['sample']}/config.txt",
        f"results/{config['sample']}/raw_dtype.txt"
    output:
        temp(directory(f"results/{config['sample']}/NAT_guppy_modified/pass/")),
        temp(directory(f"results/{config['sample']}/NAT_guppy_modified/fail/")),
        runtime = f"results/{config['sample']}/NAT_guppy_modified/runtime.txt"
    log:
        f"results/{config['sample']}/NAT_guppy_modified/guppy.log"
    conda:
        "../envs/pod5.yaml"
    resources:
        gpu = 1
    params: 
        raw_dtype = parse_raw_dtype,
        raw_dir = lambda wildcards, input: f"{input[1]}",#get_raw_dir,
        config = parse_config,
        mods = lambda wildcards: get_guppy_models_from_config(wildcards, "modifications"),
    shell:
        ("""
        config={params.config}
        mods={params.mods}
        start=`date +%s`
        if [ ! -z "$config" ]; then
            if [ ! -z "$mods" ]; then
                mkdir -p results/{config[sample]}/NAT_guppy_modified
                if [[ {params.raw_dtype} == pod5 ]]; then
                    if [ ! -d "./{params.raw_dir}/fast5s/" ]; then
                        mkdir -p ./{params.raw_dir}/fast5s/
                        pod5 convert to_fast5 ./{params.raw_dir}/*.fast5 --output ./{params.raw_dir}/fast5s/ 2>&1 | tee -a {log}
                    fi
                    resources/ont-guppy/bin/guppy_basecaller -i {params.raw_dir}/fast5s/ \
                        -s results/{config[sample]}/NAT_guppy_modified/ \
                        -c {params.mods} --recursive -x 'cuda:0' --bam_out 2>&1 | tee -a {log}
                else
                    resources/ont-guppy/bin/guppy_basecaller -i {params.raw_dir} \
                        -s results/{config[sample]}/NAT_guppy_modified/ \
                        -c {params.mods} --recursive -x 'cuda:0' --bam_out 2>&1 | tee -a {log}
                fi
            fi
        fi
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)