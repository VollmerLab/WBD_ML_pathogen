from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, ConfusionMatrixDisplay
import matplotlib.pyplot as plt

from Code.Python_Data_Wrangling import *
from Code.TankDataWrangling import *

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

from sklearn import metrics
#define metrics
y_pred_proba = lr_model.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)

#create ROC curve
plt.plot(fpr,tpr,label="AUC="+str(auc))
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc=4)
plt.show()

from sklearn.metrics import confusion_matrix
disp = ConfusionMatrixDisplay(confusion_matrix=confusion_matrix(y_test, y_pred), display_labels=lr_model.classes_)
disp.plot(cmap=plt.cm.terrain)
plt.show()


from sklearn.inspection import permutation_importance

# Compute permutation feature importance
result = permutation_importance(lr_model, X_test, y_test, n_repeats=10, random_state=42)
perm_importance = result.importances_mean


nfeatures = 20
# Sort permutation feature importance values in descending order
sorted_indices = perm_importance.argsort()[::-1]
top_10_perm_importance = perm_importance[sorted_indices][:nfeatures]

# Get the corresponding feature names for the top 10 importance values
top_10_feature_names = X_test.columns[sorted_indices][:nfeatures]

# Plot the top 10 permutation feature importance
plt.figure()
plt.barh(top_10_feature_names, top_10_perm_importance, color='skyblue')
plt.xlabel('Permutation Feature Importance')
plt.ylabel('Feature Names')
plt.title('Top ' + str(nfeatures) + ' Permutation Feature Importance')
plt.gca().invert_yaxis()
plt.subplots_adjust(left=.5)
plt.show()

gini_importance = abs(lr_model.coef_[0])
# Sort Gini feature importance values in descending order
gini_sorted_indices = gini_importance.argsort()[::-1]
top_10_gini_importance = gini_importance[gini_sorted_indices][:nfeatures]

# Get the corresponding feature names for the top 10 importance values
top_10_feature_names = X_train.columns[gini_sorted_indices][:nfeatures]

# Plot the top 10 Gini feature importance
plt.figure()
plt.barh(top_10_feature_names, top_10_gini_importance, color='lightcoral')
plt.xlabel('Gini Feature Importance')
plt.ylabel('Feature Names')
plt.title('Top ' + str(nfeatures) + ' Gini Feature Importance')
plt.gca().invert_yaxis()
plt.subplots_adjust(left=.5)
plt.show()


import shap

explainer = shap.Explainer(lr_model, X_train)
shap_values = explainer(X_test)

shap.plots.beeswarm(shap_values, show=False, max_display=21)
fig, ax = plt.gcf(), plt.gca()
plt.subplots_adjust(left=.5)
plt.show()


# Generate SHAP force plots for a specific prediction (you can change the index)

indexes = {0, 10}
for i in indexes:
    print(y_test.iloc[i])
    plot = shap.plots.force(shap_values[i], feature_names=X_test.columns, show=False)
    shap.save_html('healthySample' + str(i) + '.html', plot)

indexes = {2, 122}
for i in indexes:
    print(y_test.iloc[i])
    plot = shap.plots.force(shap_values[i], feature_names=X_test.columns, show=False)
    shap.save_html('diseasedSample' + str(i) + '.html', plot)