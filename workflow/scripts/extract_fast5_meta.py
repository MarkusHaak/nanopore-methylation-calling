import os
import glob
import yaml

from ont_fast5_api.fast5_interface import get_fast5_file

#fast5_fns = glob.glob(snakemake.input, '')
#assert len(fast5_fns) >= 1, "No fast"

#fast5_fn = str(snakemake.input)
fast5_fn = str(snakemake.params.any_fast5)
out_fn = str(snakemake.output[0])
out_dtype_fn = str(snakemake.output[1])
os.makedirs(os.path.dirname(out_fn), exist_ok=True)

with get_fast5_file(fast5_fn, mode="r") as f5:
    for r in f5.get_reads():
        read = r
        break
    with open(out_fn, "w+") as ff:
        # extract metadata
        tracking_id = read.get_tracking_id()
        yaml.dump(tracking_id, ff)
        context_tags = read.get_context_tags()
        yaml.dump(context_tags, ff)

with open(out_dtype_fn, "w") as f:
    print("fast5", file=f)