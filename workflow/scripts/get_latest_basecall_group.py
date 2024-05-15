import os, sys
import glob
from ont_fast5_api.fast5_interface import get_fast5_file

def main(fast5_fn, out_fn):
    if fast5_fn != "missing":
        with get_fast5_file(fast5_fn, mode="r") as f5:
            for r in f5.get_reads():
                read = r
                break
            group = read.get_latest_analysis('Basecall_1D')
            attr = read.get_analysis_attributes(group)
            if out_fn:
                with open(out_fn, "w") as f:
                    print(group, file=f)
            else:
                print(group)
                for key,val in attr.items():
                    print(key,":",val)
                #print('model_version_id :', attr.get('model_version_id'),',', 'time_stamp :', attr.get('time_stamp'))
    else:
        raise ValueError("ERROR: No fast5 input file provided.")


if __name__ == '__main__':
    try:
        fast5_fn = str(snakemake.params.any_fast5)
        out_fn = str(snakemake.out_fn)
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
        if len(sys.argv) > 2:
            out_fn = sys.argv[2]
        else:
            out_fn = None

    main(fast5_fn, out_fn)