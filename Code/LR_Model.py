from Python_Data_Wrangling import *
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import RandomizedSearchCV

'''
 The best estimator across ALL searched params:
 LinearRegression(fit_intercept=False, n_jobs=17)

 The best score across ALL searched params:
 0.5851343778714985

 The best parameters across ALL searched params:
 {'positive': False, 'n_jobs': 17, 'fit_intercept': False, 'copy_X': True}
The accuracy of the model is: 0.10039821182346598
'''

# Create the random grid
random_grid = {'copy_X': [True, False],
               'fit_intercept': [True, False],
               'n_jobs': range(1, 31),
               'positive': [True, False]
               }

# Use the random grid to search for best hyperparameters
# First create the base model to tune
lr = LinearRegression()
print(lr.get_params())
# Random search of parameters, using 3 fold cross validation,
# search across 100 different combinations, and use all available cores
lr_random = RandomizedSearchCV(estimator = lr, param_distributions = random_grid, n_iter = 100, cv = 3, verbose=2, random_state=42, n_jobs = -1)

# Fit the random search model
lr_random.fit(X_train, y_train.values.ravel())

print("\n The best estimator across ALL searched params:\n", lr_random.best_estimator_)
print("\n The best score across ALL searched params:\n", lr_random.best_score_)
print("\n The best parameters across ALL searched params:\n", lr_random.best_params_)

y_pred = lr_random.predict(X_test)

# Prediction accuracy
print('The accuracy of the model is: {}'.format(lr_random.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))