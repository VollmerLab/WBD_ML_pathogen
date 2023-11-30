from Python_Data_Wrangling import *
from sklearn.linear_model import Lasso
from sklearn.model_selection import RandomizedSearchCV

# lasso
'''
 The best estimator across ALL searched params:
 Lasso(alpha=0.1, precompute=True, tol=0.01)

 The best score across ALL searched params:
 0.804173003139232

 The best parameters across ALL searched params:
 {'tol': 0.01, 'precompute': True, 'alpha': 0.1}
The accuracy of the model is: 0.8333948631880885
'''

# Create the random grid
random_grid = {'alpha': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
               'precompute': [False, True],
               'tol': [.0001, .0005, .001, .005, .01, .05]
               }

# Use the random grid to search for best hyperparameters
# First create the base model to tune
rc = Lasso()
print(rc.get_params())
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