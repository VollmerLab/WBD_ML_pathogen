import numpy as np
from sklearn.model_selection import train_test_split
from Pycaret_Analysis import *

dataset = dataset.drop(["time", "tank", "anti_health","anti","geno","resist","plate","dataset"], axis=1)
numericalDataset = dataset.drop(["ID", "year", "season", "site", "health", "yrSeason"], axis=1)
NFI = numericalDataset.columns
NF = NFI.tolist()

y=dataset.health
x=dataset.drop(["ID", "year", "season", "site", "health", "yrSeason"], axis=1)

training, test =train_test_split(dataset,test_size=0.40)
train_x = [x for x in training]
train_y = [y for y in training]

from sklearn.ensemble import RandomForestRegressor
rf = RandomForestRegressor(random_state = 42)

from sklearn.model_selection import RandomizedSearchCV
# Number of trees in random forest
n_estimators = [int(x) for x in np.linspace(start = 200, stop = 2000, num = 10)]
# Number of features to consider at every split
max_features = ['auto', 'sqrt']
# Maximum number of levels in tree
max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
max_depth.append(None)
# Minimum number of samples required to split a node
min_samples_split = [2, 5, 10]
# Minimum number of samples required at each leaf node
min_samples_leaf = [1, 2, 4]
# Method of selecting samples for training each tree
bootstrap = [True, False]
# Create the random grid
random_grid = {'n_estimators': n_estimators,
               'max_features': max_features,
               'max_depth': max_depth,
               'min_samples_split': min_samples_split,
               'min_samples_leaf': min_samples_leaf,
               'bootstrap': bootstrap}

# Use the random grid to search for best hyperparameters
# First create the base model to tune
rf = RandomForestRegressor()
# Random search of parameters, using 3 fold cross validation,
# search across 100 different combinations, and use all available cores
rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)
# Fit the random search model
rf_random.fit(train_x, train_y)