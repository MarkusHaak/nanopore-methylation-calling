import os, sys, re
m = re.fullmatch("[^-]+-[^-]+-[^-]+-[^-]+-[^-]+.fast5", os.path.basename(sys.argv[1]))
if m is not None:
    print("True")
else:
    print("False")