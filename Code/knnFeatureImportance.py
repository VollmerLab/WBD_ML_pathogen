import shap
from matplotlib import pyplot as plt
from sklearn.neighbors import KNeighborsClassifier
from sklearn.inspection import permutation_importance

from Python_Data_Wrangling import *

knn = KNeighborsClassifier(metric='manhattan', n_neighbors=4)

knn.fit(X_train, y_train.values.ravel())

result = permutation_importance(
    knn, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)

knn_importances = pd.Series(result.importances_mean, index=features).nlargest(10)

fig, ax = plt.subplots()
knn_importances.plot.bar(ax=ax)
ax.set_title("Feature importances using permutation on KNN model")
ax.set_ylabel("Mean accuracy decrease")
plt.subplots_adjust(bottom=0.5)
plt.show()


explainer = shap.KernelExplainer(knn.predict_proba, X_train)
shap_values = explainer.shap_values(X_test)
shap.plots.waterfall(shap_values)
