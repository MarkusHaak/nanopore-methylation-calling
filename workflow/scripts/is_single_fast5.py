#import os, sys, re
#m = re.fullmatch("[^-]+-[^-]+-[^-]+-[^-]+-[^-]+.fast5", os.path.basename(sys.argv[1]))
#if m is not None:
#    print("True")
#else:
#    print("False")

import os, sys
from ont_fast5_api.fast5_interface import get_fast5_file

def main(fast5_fn):
    if fast5_fn != "missing":
        with get_fast5_file(fast5_fn, mode="r") as f5:
            rids = f5.get_read_ids()
        if len(rids) == 1:
            print("True")
        else:
            print("False")
    else:
        raise ValueError("ERROR: No fast5 input file provided.")


if __name__ == '__main__':
    try:
        fast5_fn = str(snakemake.params.any_fast5)
    except:
        in_fn = sys.argv[1]
        if os.path.isdir(in_fn):
            fast5_fn = glob.glob(in_fn.rstrip('/') + "/**/*.fast5", recursive=True)
            if len(fast5_fn) == 0:
                fast5_fn = 'missing'
            else:
                fast5_fn = fast5_fn[0]
        elif in_fn.endswith('.fast5'):
            fast5_fn = sys.argv[1]
        else:
            fast5_fn = 'missing'

    main(fast5_fn)