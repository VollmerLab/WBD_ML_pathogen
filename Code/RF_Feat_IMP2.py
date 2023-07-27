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
