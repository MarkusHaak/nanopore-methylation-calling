rule snapper:
    input:
        nat_fast5 = f"results/{config['sample']}/NAT_guppy_canonical/merged_single_fast5s/",
        wga_fast5 = f"results/{config['sample']}/WGA_guppy_canonical/merged_single_fast5s/",
        reference = f"results/{config['sample']}/reference.fasta"
    output:
        out_dir = directory(f"results/{config['sample']}/snapper/"),
        runtime = f"results/{config['sample']}/snapper/runtime.txt"
    conda:
        "../envs/snapper.yaml"
    threads:
        config['threads']
    shell: # requires that output dir does not exist
        ("""
        if [ -d "{output.out_dir}" ]; then
            rm -r {output.out_dir}
        fi
        start=`date +%s`
        snapper \
            -sample_fast5dir {input.nat_fast5} \
            -control_fast5dir {input.wga_fast5} \
            -reference {input.reference} \
            -threads {threads} \
            -outdir {output.out_dir} \
            -threads {threads}
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)
