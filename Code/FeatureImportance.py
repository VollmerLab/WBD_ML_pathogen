import shap
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn import metrics

from Python_Data_Wrangling import *

rf = RandomForestClassifier(max_depth=20, min_samples_leaf=2, min_samples_split=5, n_estimators=1000)


rf.fit(X_train, y_train)

result = permutation_importance(
    rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)

rf_importances = pd.Series(result.importances_mean, index=features).nlargest(10)

fig, ax = plt.subplots()
rf_importances.plot.bar(ax=ax)
ax.set_title("Feature importances using permutation on Random Forest model")
ax.set_ylabel("Mean accuracy decrease")
plt.subplots_adjust(bottom=0.5)
plt.show()


#define metrics
y_pred_proba = rf.predict_proba(X_test)[::,1]
fpr, tpr, _ = metrics.roc_curve(y_test,  y_pred_proba)
auc = metrics.roc_auc_score(y_test, y_pred_proba)

#create ROC curve
plt.plot(fpr,tpr,label="AUC="+str(auc))
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.legend(loc=4)
plt.show()

explainer = shap.TreeExplainer(rf)
shap_values = explainer(X_train)
shap.summary_plot(shap_values, X_train.ravel())
