rule download_dorado:
    output:
        temp(f"resources/dorado-{config['dorado_version']}.tar.gz")
    shell:
        "wget -P resources/ https://cdn.oxfordnanoportal.com/software/analysis/dorado-{config[dorado_version]}.tar.gz"

rule extract_dorado:
    input:
        f"resources/dorado-{config['dorado_version']}.tar.gz"
    output:
        f"resources/dorado-{config['dorado_version']}/bin/dorado"
    #cache:
    #    True
    shell:
        "tar -xzf {input} -C resources/"

rule download_guppy:
    output:
        temp(f"resources/ont-guppy_{config['guppy_version']}.tar.gz")
    shell:
        "wget -P resources/ https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_{config[guppy_version]}.tar.gz"

rule extract_guppy:
    input:
        f"resources/ont-guppy_{config['guppy_version']}.tar.gz"
    output:
        f"resources/ont-guppy/bin/guppy_basecaller"
    #cache:
    #    True
    shell:
        "tar -xzf {input} -C resources/"

rule export_guppy_workflows:
    input:
        "resources/ont-guppy/bin/guppy_basecaller"
    output:
        "resources/guppy_workflows.csv"
    shell:
        ("""
        {input} --print_workflows > {output}
        """)

rule clone_rerio:
    output:
        "resources/rerio/download_model.py"
    shell:
        "git clone https://github.com/nanoporetech/rerio.git resources/rerio"

rule install_MicrobeMod:
    output:
        directory("resources/MicrobeMod/MicrobeMod/db/")
    conda:
        "../envs/MicrobeMod.yaml"
    shell:
        ("""
        cd resources/
        rm -rf MicrobeMod
        git clone https://github.com/cultivarium/MicrobeMod.git
        cd MicrobeMod/MicrobeMod/
        python download_db.py
        cd ../
        pip install .
        """)
