import shap
import imgkit
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier


from Python_Data_Wrangling import *

rf = RandomForestClassifier(max_depth=20, min_samples_leaf=2, min_samples_split=5, n_estimators=1000)

rf.fit(X_train, y_train)

# Initialize the SHAP explainer
explainer = shap.TreeExplainer(rf)

# Calculate SHAP values for the test data
shap_values = explainer.shap_values(X_test)

from sklearn.inspection import permutation_importance

# Compute permutation feature importance
result = permutation_importance(rf, X_test, y_test, n_repeats=10, random_state=42)
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

gini_importance = abs(rf.coef_[0])
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

explainer = shap.Explainer(rf, X_train)
shap_values = explainer(X_test)

shap.plots.beeswarm(shap_values, show=False)
fig, ax = plt.gcf(), plt.gca()
plt.subplots_adjust(left=.5)
plt.show()

# Choose a specific data point from the test set (e.g., the first data point)
data_point_index = 0

# Visualize the SHAP values for that specific data point
shap.initjs()
shap.force_plot(explainer.expected_value[1], shap_values[1][data_point_index], X_test.iloc[data_point_index])

# Summary plot for feature importance
shap.summary_plot(shap_values[1], X_test)

# Calculate the average SHAP values for each feature
mean_shap_values = np.abs(shap_values[1]).mean(axis=0)

# Get the indices of the top 10 features
top_10_indices = np.argsort(mean_shap_values)[-10:]

# Get the names of the top 10 features
top_10_features = X_test.columns[top_10_indices]

# Get the corresponding mean SHAP values for the top 10 features
top_10_shap_values = mean_shap_values[top_10_indices]

# Create a bar plot to visualize feature importance for the top 10 features
plt.figure()
plt.barh(top_10_features, top_10_shap_values)
plt.xlabel("Mean |SHAP Value|")
plt.ylabel("Features", labelpad=20)
plt.title("Top 10 Feature Importance - Mean SHAP Values")
plt.subplots_adjust(left=.5)
plt.show()


shap.plots.beeswarm(shap_values, show=False)
fig, ax = plt.gcf(), plt.gca()
plt.subplots_adjust(left=.5)
plt.show()


# Choose 10-15 instances from the test set for which you want to create the force plot
selected_instances = X_test.iloc[:15]  # Change the slice to select desired instances

# Calculate SHAP values for the selected instances
shap_values = explainer.shap_values(selected_instances)

# Create the force plot for the selected instances
shap.initjs()
force_plot = shap.force_plot(explainer.expected_value[1], shap_values[1], selected_instances)

# Save the force plot as an HTML file
shap.save_html("force_plot.html", force_plot)

# Convert the HTML file to a PNG image
imgkit.from_file("force_plot.html", "force_plot.png")
