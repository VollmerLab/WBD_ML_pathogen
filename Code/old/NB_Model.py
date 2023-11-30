from Python_Data_Wrangling import *

from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import RandomizedSearchCV
'''
The accuracy of the model is: 0.9274193548387096
[[36  6]
 [ 3 79]]
'''


# First create the base model to tune
nb = GaussianNB()

# Fit the random search model
nb.fit(X_train, y_train.values.ravel())


y_pred = nb.predict(X_test)

# Prediction accuracy
print('The accuracy of the model is: {}'.format(nb.score(X_test, y_test)))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))