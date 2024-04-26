import re

from get_workflow import get_guppy_config_name, config_name_to_pore

config_name = get_guppy_config_name(snakemake.input[0], snakemake.input[1])
print("config name:", config_name)
with open(f"{snakemake.output[0]}", 'w') as f:
    print(config_name, file=f)

pore = config_name_to_pore(config_name)
print("pore:", pore)
with open(f"{snakemake.output[1]}", 'w') as f:
    print(pore, file=f)