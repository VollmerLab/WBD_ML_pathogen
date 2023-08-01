from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt

from Code.Python_Data_Wrangling import *

best_lr_params = {
    'C': 1.0,
    'class_weight': None,
    'dual': False,
    'fit_intercept': True,
    'intercept_scaling': 1,
    'l1_ratio': None,
    'max_iter': 1000,
    'multi_class': 'auto',
    'n_jobs': None,
    'penalty': 'l2',
    'random_state': 2807,
    'solver': 'lbfgs',
    'tol': 0.0001,
    'verbose': 0,
    'warm_start': False
}

# Create the LR model with the best hyperparameters
lr_model = LogisticRegression(**best_lr_params)

# Fit the LR model to the training data
lr_model.fit(X_train, y_train)

# Make predictions on the test set
y_pred = lr_model.predict(X_test)

# Calculate accuracy
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Get classification report
print("Classification Report:")
print(classification_report(y_test, y_pred))

from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))

from sklearn.inspection import permutation_importance

# Compute permutation feature importance
result = permutation_importance(lr_model, X_test, y_test, n_repeats=10, random_state=42)
perm_importance = result.importances_mean



# Sort permutation feature importance values in descending order
sorted_indices = perm_importance.argsort()[::-1]
top_10_perm_importance = perm_importance[sorted_indices][:10]

# Get the corresponding feature names for the top 10 importance values
top_10_feature_names = X_test.columns[sorted_indices][:10]

# Plot the top 10 permutation feature importance
plt.figure()
plt.barh(top_10_feature_names, top_10_perm_importance, color='skyblue')
plt.xlabel('Permutation Feature Importance')
plt.ylabel('Feature Names')
plt.title('Top 10 Permutation Feature Importance')
plt.gca().invert_yaxis()
plt.subplots_adjust(left=.5)
plt.show()

gini_importance = abs(lr_model.coef_[0])
# Sort Gini feature importance values in descending order
gini_sorted_indices = gini_importance.argsort()[::-1]
top_10_gini_importance = gini_importance[gini_sorted_indices][:10]

# Get the corresponding feature names for the top 10 importance values
top_10_feature_names = X_train.columns[gini_sorted_indices][:10]

# Plot the top 10 Gini feature importance
plt.figure()
plt.barh(top_10_feature_names, top_10_gini_importance, color='lightcoral')
plt.xlabel('Gini Feature Importance')
plt.ylabel('Feature Names')
plt.title('Top 10 Gini Feature Importance')
plt.gca().invert_yaxis()
plt.subplots_adjust(left=.5)
plt.show()


import shap

explainer = shap.Explainer(lr_model, X_train)
shap_values = explainer(X_test)

shap.plots.beeswarm(shap_values, show=False)
fig, ax = plt.gcf(), plt.gca()
plt.subplots_adjust(left=.5)
plt.show()