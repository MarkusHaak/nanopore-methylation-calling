import os
import glob
import yaml
import re
import sys

from ont_fast5_api.fast5_interface import get_fast5_file

def main(fast5_fn, pod5_fn, out_fn, out_dtype_fn):
    if os.path.dirname(out_fn):
        os.makedirs(os.path.dirname(out_fn), exist_ok=True)
    if pod5_fn != "missing":
        tmp_meta_fn = out_fn + ".tmp"
        os.system(f"pod5 inspect debug {pod5_fn} > {tmp_meta_fn}")
        os.system(f"echo pod5 > {out_dtype_fn}")
        # convert to yaml
        parsed = {}
        with open(tmp_meta_fn, 'r') as f:
            for line in f.readlines():
                m = re.fullmatch(r"(\w+)[ ]?:[ ]?(.*)", line.strip())
                if m is not None:
                    key = m.group(1)
                    val = m.group(2).strip()
                    if val.startswith('{') and val.endswith('}'):
                        for s in val.split(','):
                            m = re.fullmatch(r"'(\w+)'[ ]?:[ ]?'(.*)'", s.strip())
                            if m is not None:
                                key = m.group(1)
                                val = m.group(2).strip()
                                parsed[key] = val
                    else:
                        parsed[key] = val
        os.remove(tmp_meta_fn)
        with open(out_fn, 'w') as f:
            for key,val in parsed.items():
                print(f"{key}: {val}", file=f)
    elif fast5_fn != "missing":
        with get_fast5_file(fast5_fn, mode="r") as f5:
            for r in f5.get_reads():
                read = r
                break
            with open(out_fn, "w+") as ff:
                # extract metadata
                tracking_id = read.get_tracking_id()
                yaml.dump(tracking_id, ff)
                context_tags = read.get_context_tags()
                if context_tags:
                    yaml.dump(context_tags, ff)

        with open(out_dtype_fn, "w") as f:
            print("fast5", file=f)

if __name__ == '__main__':
    try:
        fast5_fn = str(snakemake.params.any_fast5)
        pod5_fn = str(snakemake.params.any_pod5)
        out_fn = str(snakemake.output[0])
        out_dtype_fn = str(snakemake.output[1])
    except:
        in_fn = sys.argv[1]
        fast5_fn = in_fn if in_fn.endswith('.fast5') else 'missing'
        pod5_fn = in_fn if in_fn.endswith('.pod5') else 'missing'
        out_fn = sys.argv[2]
        out_dtype_fn = sys.argv[3]

    main(fast5_fn, pod5_fn, out_fn, out_dtype_fn)