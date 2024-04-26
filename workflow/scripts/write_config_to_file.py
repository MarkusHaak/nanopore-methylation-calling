import os, yaml

os.makedirs(os.path.dirname(f"{snakemake.output}"), exist_ok=True)
with open(f"{snakemake.output}", 'w+') as ff:
    yaml.dump(snakemake.config, ff)