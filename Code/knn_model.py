from Python_Data_Wrangling import *

'''
The best estimator across ALL searched params:
 KNeighborsClassifier(n_neighbors=1)

 The best score across ALL searched params:
 0.9517241379310345

 The best parameters across ALL searched params:
 {'weights': 'uniform', 'n_neighbors': 1}
The accuracy of the model is: 0.9838709677419355
[[41  1]
 [ 1 81]]
 '''

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import RandomizedSearchCV

# Create the random grid
random_grid = dict(n_neighbors=range(1, 31), weights=['uniform', 'distance'])

# Use the random grid to search for best hyperparameters
# First create the base model to tune
knn = KNeighborsClassifier()
# Random search of parameters, using 3 fold cross validation,
# search across 100 different combinations, and use all available cores
knn_random = RandomizedSearchCV(knn, random_grid, n_iter = 100, cv = 10, random_state=42)

# Fit the random search model
knn_random.fit(X_train, y_train)

print("\n The best estimator across ALL searched params:\n", knn_random.best_estimator_)
print("\n The best score across ALL searched params:\n", knn_random.best_score_)
print("\n The best parameters across ALL searched params:\n", knn_random.best_params_)

y_pred = knn_random.predict(X_test)

# Prediction accuracy
print('The accuracy of the model is: {}'.format(knn_random.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))