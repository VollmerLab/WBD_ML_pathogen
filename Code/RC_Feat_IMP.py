import pandas as pd
import numpy as np
from IPython.core.display_functions import display
from sklearn.model_selection import train_test_split
from sklearn.linear_model import RidgeClassifier
import shap
import matplotlib.pyplot as plt
from Python_Data_Wrangling import *

ridge_classifier = RidgeClassifier(alpha=0.2)

# Train the Ridge Classifier on the training data
ridge_classifier.fit(X_train, y_train)

# Initialize the Kernel Explainer
explainer = shap.KernelExplainer(ridge_classifier.decision_function, X_train)

# Calculate SHAP values for the test data
shap_values = explainer.shap_values(X_test)

# Calculate the average SHAP values for each feature
mean_shap_values = np.abs(shap_values).mean(axis=0)

# Get the indices of the top 10 features
top_10_indices = np.argsort(mean_shap_values)[-10:]

# Get the names of the top 10 features
top_10_features = X_test.columns[top_10_indices]

# Get the corresponding mean SHAP values for the top 10 features
top_10_shap_values = mean_shap_values[top_10_indices]

# Create a bar plot to visualize feature importance for the top 10 features
plt.figure(figsize=(10, 6))
plt.barh(top_10_features, top_10_shap_values)
plt.xlabel("Mean |SHAP Value|")
plt.ylabel("Features")
plt.title("Top 10 Feature Importance - Mean SHAP Values (Ridge Classifier)")
plt.show()

# Create a summary plot for feature importance
plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values, X_test)
plt.title("Summary Plot of Feature Importance (Ridge Classifier)")
plt.show()

# Choose 10-15 instances from the test set for which you want to create the force plot
selected_instances = X_test.iloc[:15]  # Change the slice to select desired instances

# Calculate SHAP values for the selected instances
shap_values = explainer.shap_values(selected_instances)

force_plot = shap.force_plot(explainer.expected_value, shap_values, selected_instances)

# Save the force plot as an HTML file
shap.save_html("RCforce_plot.html", force_plot)

# Convert the HTML file to a PNG image
import imgkit

imgkit.from_file("RCforce_plot.html", "RCforce_plot.png")


