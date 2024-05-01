rule nanodisco_preprocess:
    input:
        nat_fast5 = f"results/{config['sample']}/NAT_guppy_canonical/gzip/",
        wga_fast5 = f"results/{config['sample']}/WGA_guppy_canonical/gzip/",
        reference = f"results/{config['sample']}/reference.fasta"
    output:
        rds = f"results/{config['sample']}/nanodisco/analysis/ND_difference.RDS",
        runtime = f"results/{config['sample']}/nanodisco/preprocess_runtime.txt",
        preprocessed = temp(directory(f"results/{config['sample']}/nanodisco/preprocessed/")),
        difference_dir = temp(directory(f"results/{config['sample']}/nanodisco/difference/")),
        analysis_dir = directory(f"results/{config['sample']}/nanodisco/analysis/"),
        tmp_dir = temp(directory(f"results/{config['sample']}/nanodisco/tmp/"))
    log:
        f"results/{config['sample']}/nanodisco/nanodisco_preprocess.log"
    threads:
        config['threads_nanodisco']
    singularity:
        "library://fanglab/default/nanodisco"
    shell:
        ("""
        start=`date +%s`
        
        if [ -d "{output.preprocessed}" ]; then
            rm -r {output.preprocessed}
        fi
        nanodisco preprocess -p {threads} -f {input.wga_fast5} -s WGA -o {output.preprocessed} -r {input.reference} 2>>{log}
        nanodisco preprocess -p {threads} -f {input.nat_fast5} -s NAT -o {output.preprocessed} -r {input.reference} 2>>{log}
        
        mkdir -p {output.tmp_dir}
        export TMPDIR={output.tmp_dir}
        if [ -d "{output.difference_dir}" ]; then
            rm -r {output.difference_dir}
        fi
        echo "running nanodisco difference -nj {threads} -nc 1 -p 2 -i {output.preprocessed} -o {output.difference_dir} -w WGA -n NAT -r {input.reference}"
        nanodisco difference -nj {threads} -nc 1 -p 2 -i {output.preprocessed} -o {output.difference_dir} -w WGA -n NAT -r {input.reference} 2>>{log}

        if [ -d "{output.analysis_dir}" ]; then
            rm -r {output.analysis_dir}
        fi
        echo "running nanodisco merge -d {output.difference_dir} -o {output.analysis_dir} -b ND"
        nanodisco merge -d {output.difference_dir} -o {output.analysis_dir} -b ND 2>>{log}
        end=`date +%s`
        runtime=$((end-start))
        echo $runtime > {output.runtime}
        """)

#rule nanodisco_discovery:
#    input:
#        rds = f"results/{config['sample']}/nanodisco/analysis/ND_difference.RDS",
#        reference = f"results/{config['sample']}/reference.fasta",
#        tmp_dir = f"results/{config['sample']}/nanodisco/tmp/"
#    output:
#        discovery_dir = directory(f"results/{config['sample']}/nanodisco/discovery/"),
#        runtime = f"results/{config['sample']}/nanodisco/discovery_runtime.txt"
#    threads:
#        config['threads_nanodisco']
#    singularity:
#        "library://fanglab/default/nanodisco"
#    shell:
#        ("""
#        pwd
#        ls
#        echo tada
#        ls /mnt/
#        echo tadum
#        ls
#        export TMPDIR={input.tmp_dir}
#        start=`date +%s`
#        nanodisco motif -p 1 -b ND -d {input.rds} -o {output.discovery_dir} -r {input.reference} -a
#        end=`date +%s`
#        runtime=$((end-start))
#        echo $runtime > {output.runtime}
#        """)
#
#rule nanodisco_typing_fine_mapping:
#    input:
#        analysis_dir = f"results/{config['sample']}/nanodisco/analysis/",
#        rds = f"results/{config['sample']}/nanodisco/analysis/ND_difference.RDS",
#        reference = f"results/{config['sample']}/reference.fasta",
#    output:
#        motif_dir = directory(f"results/{config['sample']}/nanodisco/motif/"),
#        runtime = f"results/{config['sample']}/nanodisco/typing_fine_mapping_runtime.txt"
#    threads:
#        config['threads_nanodisco']
#    singularity:
#        "library://fanglab/default/nanodisco"
#    shell:
#        ("""
#        start=`date +%s`
#        nanodisco characterize -p {threads} -b ND -d {input.rds} -o {output.motif_dir} -m GATC,CCWGG,GCACNNNNNNGTT,AACNNNNNNGTGC -t nn,rf,knn -r {input.reference}
#        end=`date +%s`
#        runtime=$((end-start))
#        echo $runtime > {output.runtime}
#        """)
