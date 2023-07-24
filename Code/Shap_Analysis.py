from RF_Model import *
import shap

gg = df.iloc[1]
explainer = shap.TreeExplainer(tuned_rf)

shap_values = explainer.shap_values(numericalDataset)
print(shap_values)