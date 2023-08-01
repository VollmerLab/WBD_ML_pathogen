import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from Python_Data_Wrangling import *

"""
The best estimator across ALL searched params:
RandomForestClassifier(max_depth=20, min_samples_leaf=2, min_samples_split=5, n_estimators=1000)

The best score across ALL searched params:
0.9758734249713631

The best parameters across ALL searched params:
{'n_estimators': 1000, 'min_samples_split': 5, 'min_samples_leaf': 2, 'max_features': 'sqrt', 'max_depth': 20, 'bootstrap': True}
The accuracy of the model is: 0.967741935483871

[[41  1]
 [ 2 80]]
"""

from sklearn.ensemble import RandomForestClassifier
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
rf = RandomForestClassifier()
# Random search of parameters, using 3 fold cross validation,
# search across 100 different combinations, and use all available cores
rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)

# Fit the random search model
rf_random.fit(X_train, y_train.values.ravel())

print("\n The best estimator across ALL searched params:\n", rf_random.best_estimator_)
print("\n The best score across ALL searched params:\n", rf_random.best_score_)
print("\n The best parameters across ALL searched params:\n", rf_random.best_params_)

y_pred = rf_random.predict(X_test)

# Prediction accuracy
print('The accuracy of the model is: {}'.format(rf_random.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))

pcr_rf = RandomForestClassifier(**{'bootstrap': True, 'ccp_alpha': 0.0, 'class_weight': None, 'criterion': 'gini', 'max_depth': None, 'max_features': 'sqrt', 'max_leaf_nodes': None, 'max_samples': None, 'min_impurity_decrease': 0.0, 'min_samples_leaf': 1, 'min_samples_split': 2, 'min_weight_fraction_leaf': 0.0, 'n_estimators': 100, 'n_jobs': -1, 'oob_score': False, 'random_state': 3802, 'verbose': 0, 'warm_start': False})

pcr_rf.fit(X_train, y_train)

y_pred = pcr_rf.predict(X_test)

# Prediction accuracy
print('The accuracy of the model is: {}'.format(pcr_rf.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))