import sys, re
m = re.fullmatch("[^-]+-[^-]+-[^-]+-[^-]+-[^-]+.fast5", sys.argv[1])
if m is not None:
    print("True")
else:
    print("False")