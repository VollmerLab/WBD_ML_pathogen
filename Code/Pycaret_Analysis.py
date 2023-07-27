import pycaret
import pandas as pd
import mlflow
from Python_Data_Wrangling import *

dataset = fullSampledf

data = dataset.sample(frac=0.6, random_state=786)
data_unseen = dataset.drop(data.index)
data.reset_index(inplace=True, drop=True)
data_unseen.reset_index(inplace=True, drop=True)
print('Data for Modeling: ' + str(data.shape))
print('Unseen Data For Predictions: ' + str(data_unseen.shape))


from pycaret.classification import *
# feature_selection when set to True, a subset of features are selected using a combination of various permutation importance techniques
# feature_selection_threshold is the threshold for feature selection
# Initialize PyCaret
clf1 = setup(dataset, target='health')


best_model = compare_models()
print(best_model.head())