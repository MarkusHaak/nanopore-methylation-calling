rule download_dorado_models:
    input:
        f"resources/dorado-{config['dorado_version']}/bin/dorado"
    output:
        directory("resources/dorado_models/{model}/")
    shell:
        "mkdir -p resources/dorado_models && resources/dorado-{config[dorado_version]}/bin/dorado download --model {wildcards.model} --directory resources/dorado_models/"

rule download_rerio_dorado_models:
    input:
        "resources/rerio/download_model.py"
    output:
        directory("resources/rerio/dorado_models/{model}/")
    shell:
        "python {input} {output}_url"

rule download_rerio_guppy_models:
    input:
        "resources/rerio/download_model.py"
    output:
        "resources/ont-guppy/data/{model}.cfg",
        "resources/ont-guppy/data/{model}.jsn"
    shell:
        ("""
        if [ ! -f resources/ont-guppy/data/{wildcards.model}.cfg ]
        then
            python {input} resources/rerio/basecall_models/{wildcards.model}
            cp resources/rerio/basecall_models/{wildcards.model}.* resources/ont-guppy/data/
        fi
        """)

rule download_required_models:
    input:
        all_model_paths()
    output:
        "resources/all_models_downloaded.txt"
    shell:
        "touch {output}"