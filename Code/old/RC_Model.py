from Python_Data_Wrangling import *
from sklearn.linear_model import RidgeClassifier
from sklearn.model_selection import RandomizedSearchCV

# ridge, elastic net, lasso
'''
 The best estimator across ALL searched params:
 RidgeClassifier(alpha=0.2)

 The best score across ALL searched params:
 0.9342425544100802

 The best parameters across ALL searched params:
 {'alpha': 0.2}
The accuracy of the model is: 0.8548387096774194
[[35  7]
 [11 71]]
'''

# Create the random grid
random_grid = {'alpha': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
               }

# Use the random grid to search for best hyperparameters
# First create the base model to tune
rc = RidgeClassifier()
# Random search of parameters, using 3 fold cross validation,
# search across 100 different combinations, and use all available cores
rc_random = RandomizedSearchCV(estimator = rc, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)

# Fit the random search model
rc_random.fit(X_train, y_train.values.ravel())

print("\n The best estimator across ALL searched params:\n", rc_random.best_estimator_)
print("\n The best score across ALL searched params:\n", rc_random.best_score_)
print("\n The best parameters across ALL searched params:\n", rc_random.best_params_)

y_pred = rc_random.predict(X_test)

# Prediction accuracy
print('The accuracy of the model is: {}'.format(rc_random.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))