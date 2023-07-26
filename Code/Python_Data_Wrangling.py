import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

data_path = '/Users/milesvollmer/Panama_Tank_Field/intermediate_files/normalized_field_asv_counts.csv'

df = pd.read_csv(data_path)

df = df.drop(['cpm_norm', 'cpm', 'n_reads', 'lib.size', 'norm.factors'], axis=1)
df["health"] = np.where(df["health"] == "D", 0, 1)
df['asv_id'] = df['asv_id'].astype(str)
df['family'] = df['family'].astype(str)
df['genus'] = df['genus'].astype(str)
df['asv_id'] = df[['asv_id', 'family', 'genus']].agg('-'.join, axis=1)


# sample as X axis, health and ASVs as Y
sampledf = df.drop(['year', 'season', 'site', 'dataset', 'domain', 'phylum', 'class', 'order', 'family', 'genus'], axis=1)
healthdf = sampledf.drop(['asv_id', 'log2_cpm_norm'], axis=1)
healthdf = healthdf.groupby('sample_id').max()

sampledf = sampledf.pivot(index= 'sample_id', columns= 'asv_id', values='log2_cpm_norm')
fullSampledf = pd.concat([sampledf, healthdf], axis=1)

dataset = fullSampledf

features = dataset.drop(['health'], axis=1).columns
target = ['health']
y=dataset.loc[:, target]
X=dataset.loc[:, features]

X_train, X_test, y_train, y_test = train_test_split(X,              #the input features
                                                    y,              #the label
                                                    test_size=0.3,  #set aside 30% of the data as the test set
                                                    random_state=7, #reproduce the results
                                                   )



