import glob
import shutil
import os
# find any fasta file anywhere in the sample
fns = (
    glob.glob(f"{snakemake.config['datasets_dir']}/{snakemake.config['sample']}/reference/**/*.fasta", recursive=True) +
    glob.glob(f"{snakemake.config['datasets_dir']}/{snakemake.config['sample']}/reference/**/*.fas", recursive=True) + 
    glob.glob(f"{snakemake.config['datasets_dir']}/{snakemake.config['sample']}/reference/**/*.fna", recursive=True) + 
    glob.glob(f"{snakemake.config['datasets_dir']}/{snakemake.config['sample']}/reference/**/*.fa", recursive=True))
if len(fns) > 1:
    raise IOError(f"Found more than one potential reference fasta file: {fns}.")
elif len(fns) == 0:
    raise IOError("Found no potential reference fasta file.")
if not os.path.exists(f"results/{snakemake.config['sample']}/"):
    os.makedirs(f"results/{snakemake.config['sample']}/")
with open(f"{snakemake.output.ref_fp}", "w") as f:
    print(fns[0], file=f)
shutil.copyfile(fns[0], f"{snakemake.output.ref}")