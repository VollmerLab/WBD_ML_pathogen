import shap
import imgkit
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier


from Python_Data_Wrangling import *

rf = RandomForestClassifier(max_depth=20, min_samples_leaf=2, min_samples_split=5, n_estimators=1000)

rf.fit(X_train, y_train)

from sklearn.inspection import permutation_importance

# Compute permutation feature importance
result = permutation_importance(rf, X_test, y_test, n_repeats=30, random_state=42)
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

gini_importance = rf.feature_importances_
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


# Create a TreeExplainer
explainer = shap.TreeExplainer(rf, X_train)

# Calculate SHAP values for the test set
shap_values = explainer(X_test)
shap.summary_plot(shap_values)
print(shap_values)