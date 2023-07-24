import pandas as pd
import numpy as np

data_path = '/Users/milesvollmer/Panama_Tank_Field/intermediate_files/normalized_field_asv_counts.csv'

df = pd.read_csv(data_path)

df = df.drop(['cpm_norm', 'cpm', 'n_reads', 'lib.size', 'norm.factors'], axis=1)
df["health"] = np.where(df["health"] == "D", 0, 1)

print(df.columns)
# sample as X axis, health and ASVs as Y

sampledf = df.drop(['year', 'season', 'site', 'dataset', 'domain', 'phylum', 'class', 'order', 'family', 'genus'], axis=1)
healthdf = sampledf.drop(['asv_id', 'log2_cpm_norm'], axis=1)
healthdf = healthdf.groupby('sample_id').max()
print(healthdf.head())
sampledf = sampledf.pivot(index= 'sample_id', columns= 'asv_id', values='log2_cpm_norm')
print(sampledf.head())
fullSampledf = pd.concat([sampledf, healthdf], axis=1)
print(fullSampledf.head())





