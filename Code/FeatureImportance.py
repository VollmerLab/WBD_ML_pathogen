import shap
from matplotlib import pyplot as plt
from sklearn.inspection import permutation_importance

from RF_Model import *

# not sorting values correctly

result = permutation_importance(
    rf_random, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)

rf_importances = pd.Series(result.importances_mean, index=features).nlargest(10)

fig, ax = plt.subplots()
rf_importances.plot.bar(ax=ax)
ax.set_title("Feature importances using permutation on full model")
ax.set_ylabel("Mean accuracy decrease")
fig.tight_layout()
plt.show()


explainer = shap.TreeExplainer(rf_random)
shap_values = explainer.shap_values(X_test).nlargest(10)

shap.summary_plot(shap_values, X_test, plot_type="bar")
