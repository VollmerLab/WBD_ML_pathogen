import shap
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance

from Python_Data_Wrangling import *

rf = RandomForestClassifier(max_depth=20, min_samples_leaf=2, min_samples_split=5, n_estimators=1000)
# not sorting values correctly

rf.fit(X_train, y_train.values.ravel())

result = permutation_importance(
    rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)

rf_importances = pd.Series(result.importances_mean, index=features).nlargest(10)

fig, ax = plt.subplots()
rf_importances.plot.bar(ax=ax)
ax.set_title("Feature importances using permutation on Random Forest model")
ax.set_ylabel("Mean accuracy decrease")
fig.tight_layout()
plt.show()


explainer = shap.TreeExplainer(rf)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values, X_test, plot_type="bar")
