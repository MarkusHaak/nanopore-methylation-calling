import pandas as pd

df = pd.read_csv(f"{snakemake.input}", sep="\t", header=None, names=['contig', 'pos', 'cov'])
d = df.groupby('contig').agg(['sum','max'])
con = pd.DataFrame(d[('cov', 'sum')] / d[('pos', 'max')], columns=['depth'])
con['size'] = d[('pos', 'max')]
for i in snakemake.config['contig_coverages']:
    con[f'min_{i}'] = df.groupby('contig')['cov'].apply(lambda c: (c >= i).sum() / len(c) * 100.)
con.to_csv(f"{snakemake.output}", float_format='%.2f', sep='\t')