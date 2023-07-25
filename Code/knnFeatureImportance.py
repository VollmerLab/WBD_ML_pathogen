import shap
from matplotlib import pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
from sklearn.inspection import permutation_importance

from Python_Data_Wrangling import *

knn = KNeighborsClassifier(n_neighbors=1, weights='uniform')

knn.fit(X_train, y_train.values.ravel())

result = permutation_importance(
    knn, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)

knn_importances = pd.Series(result.importances_mean, index=features).nlargest(10)

fig, ax = plt.subplots()
knn_importances.plot.bar(ax=ax)
ax.set_title("Feature importances using permutation on KNN model")
ax.set_ylabel("Mean accuracy decrease")
fig.tight_layout()
plt.show()


explainer = shap.TreeExplainer(knn)
shap_values = explainer.shap_values(X_test)
shap.summary_plot(shap_values, X_test, plot_type="bar")
