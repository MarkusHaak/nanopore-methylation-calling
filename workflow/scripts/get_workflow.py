import sys
import yaml
import pandas as pd
import re

def get_guppy_config_name(metadata_fn, guppy_workflows_fn):
    columns = ["flowcell", "kit", "barcoding", "config_name", "model_version"]
    data = []
    with open(guppy_workflows_fn, 'r') as f:
        for _ in range(4):
            f.readline()
        for line in f.readlines():
            line = line.strip().split()
            if len(line) == 5:
                data.append(line)
            elif len(line) == 4:
                data.append(line[:2] + [None] + line[2:])
            elif len(line) == 0:
                pass
            else:
                raise ValueError("Unexcepted number of inputs")
    df = pd.DataFrame(data, columns=columns)
    
    with open(metadata_fn, 'r') as f:
        meta = yaml.safe_load(f)
    flowcell = meta.get('flow_cell_product_code', None)
    if flowcell:
        flowcell = flowcell.upper() # guppy worklow entries are always uppercase
    sequencing_kit = meta.get('sequencing_kit', None)
    if sequencing_kit:
        sequencing_kit = sequencing_kit.upper() # guppy worklow entries are always uppercase
    if flowcell is None:
        print("WARNING: flow_cell_product_code field missing in fast5 (--> v2.0), assuming FLO-MIN106")
        flowcell = "FLO-MIN106"
    elif (df.flowcell == flowcell).sum() == 0:
        raise ValueError(f"WARNING: unknown flowcell {flowcell}")
        #flowcell = "FLO-MIN106"
    kit_known = True
    if sequencing_kit is None:
        print("WARNING: kit not found in meta data")
    elif (df.kit == sequencing_kit.upper()).sum() == 0:
        print(f"WARNING: unknown kit {sequencing_kit.upper()}")
        kit_known = False

    print(f"Deduced: Flowcell: {flowcell}, kit: {sequencing_kit}")

    if sequencing_kit and kit_known:
        d = df.loc[(df.flowcell == flowcell.upper()) & (df.kit == sequencing_kit.upper()) & (df.config_name.str.startswith('dna'))]
    else:
        d = df.loc[(df.flowcell == flowcell.upper()) & (df.config_name.str.startswith('dna'))]
        print(f"WARNING: sequencing kit missing or unknown, slect first config for flowcell {flowcell}")
    if len(d) > 0:
        # return first entry if multiple match
        return d.iloc[0]['config_name']
    elif len(d) == 0:
        raise ValueError(f"No configs for flowcell {flowcell} (and kit {sequencing_kit}) listed in Guppy workflows")
    return d['config_name']

def config_name_to_pore(config_name):
    return re.search('dna_([^_]+)', config_name).group(1)

if __name__ == "__main__":
    config_name = get_guppy_config_name(sys.argv[1], sys.argv[2])
    print("config name:", config_name)
    pore = config_name_to_pore(config_name)
    print("pore:", pore)