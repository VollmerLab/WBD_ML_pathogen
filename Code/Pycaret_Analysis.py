import pycaret
import pandas as pd

data_path = ''

df = pd.read_csv(data_path)

dataset = df
dataset = dataset.drop(["time", "tank", "anti_health","anti","geno","resist","plate","dataset"], axis=1)

data = dataset.sample(frac=0.7, random_state=786)
data_unseen = dataset.drop(data.index)
data.reset_index(inplace=True, drop=True)
data_unseen.reset_index(inplace=True, drop=True)
print('Data for Modeling: ' + str(data.shape))
print('Unseen Data For Predictions: ' + str(data_unseen.shape))

numericalDataset = dataset.drop(["ID", "year", "season", "site", "health", "yrSeason"], axis=1)
NFI = numericalDataset.columns
NF = NFI.tolist()

from pycaret.classification import *
# feature_selection when set to True, a subset of features are selected using a combination of various permutation importance techniques
# feature_selection_threshold is the threshold for feature selection
exp_clf101 = setup(data = data, target = 'health', session_id=123, numeric_features=NF, feature_selection = True, feature_selection_threshold = 0.8)

best_model = compare_models()
print(best_model.head())