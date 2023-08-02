import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

data_path = '/Users/milesvollmer/Panama_Tank_Field/intermediate_files/normalized_tank_asv_counts.csv'

tdf = pd.read_csv(data_path)

tdf = tdf.drop(['cpm_norm', 'cpm', 'n_reads', 'lib.size', 'norm.factors', 'tank', 'anti_health', 'geno', 'plate', 'species'], axis=1)
tdf = tdf.where(tdf['anti'] == 'N')
tdf = tdf.where(tdf['health'] == 'D')
tdf["resist"] = np.where(tdf["resist"] == "S", 1, 0)
tdf['asv_id'] = tdf['asv_id'].astype(str)
tdf['family'] = tdf['family'].astype(str)
tdf['genus'] = tdf['genus'].astype(str)
tdf['asv_id'] = tdf[['asv_id', 'family', 'genus']].agg('-'.join, axis=1)
tdf['health'] = tdf['resist']
tdf = tdf.dropna()
ex8tdf = tdf.where(tdf['time'] == '8_exp').dropna()
ex2tdf = tdf.where(tdf['time'] == '2_exp').dropna()

# sample as X axis, health and ASVs as Y
sampletdf = ex2tdf.drop(['year', 'season', 'site', 'dataset', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'resist', 'anti', 'time'], axis=1)
healthtdf = sampletdf.drop(['asv_id', 'log2_cpm_norm'], axis=1)
healthtdf = healthtdf.groupby('sample_id').max()
sampletdf = sampletdf.pivot(index= 'sample_id', columns= 'asv_id', values='log2_cpm_norm')
fullSampletdf = pd.concat([sampletdf, healthtdf], axis=1)

tdatFeatures = fullSampletdf.drop(['health'], axis=1).columns
tdatTarget = ['health']
tdat_y=fullSampletdf.loc[:, tdatTarget]
tdat_X=fullSampletdf.loc[:, tdatFeatures]

'''
LR Model
full tdf:
Accuracy: 0.6176470588235294
Classification Report:
              precision    recall  f1-score   support

           0       0.71      0.53      0.61        19
           1       0.55      0.73      0.63        15

    accuracy                           0.62        34
   macro avg       0.63      0.63      0.62        34
weighted avg       0.64      0.62      0.62        34

[[10  9]
 [ 4 11]]

ex8tdf:
Accuracy: 0.5625
Classification Report:
              precision    recall  f1-score   support

         0.0       0.71      0.50      0.59        10
         1.0       0.44      0.67      0.53         6

    accuracy                           0.56        16
   macro avg       0.58      0.58      0.56        16
weighted avg       0.61      0.56      0.57        16

[[5 5]
 [2 4]]

ex2tdf:
Accuracy: 0.6666666666666666
Classification Report:
              precision    recall  f1-score   support

         0.0       0.71      0.56      0.63         9
         1.0       0.64      0.78      0.70         9

    accuracy                           0.67        18
   macro avg       0.68      0.67      0.66        18
weighted avg       0.68      0.67      0.66        18

[[5 4]
 [2 7]]
 
 
 RF Model
 full tdf:
The accuracy of the model is: 0.5882352941176471
Classification Report:
              precision    recall  f1-score   support

           0       0.67      0.53      0.59        19
           1       0.53      0.67      0.59        15

    accuracy                           0.59        34
   macro avg       0.60      0.60      0.59        34
weighted avg       0.60      0.59      0.59        34

[[10  9]
 [ 5 10]]

ex2tdf:
The accuracy of the model is: 0.6666666666666666
Classification Report:
              precision    recall  f1-score   support

         0.0       0.71      0.56      0.63         9
         1.0       0.64      0.78      0.70         9

    accuracy                           0.67        18
   macro avg       0.68      0.67      0.66        18
weighted avg       0.68      0.67      0.66        18

[[5 4]
 [2 7]]

ex8tdf:
The accuracy of the model is: 0.5
Classification Report:
              precision    recall  f1-score   support

         0.0       0.62      0.50      0.56        10
         1.0       0.38      0.50      0.43         6

    accuracy                           0.50        16
   macro avg       0.50      0.50      0.49        16
weighted avg       0.53      0.50      0.51        16

[[5 5]
 [3 3]]


'''