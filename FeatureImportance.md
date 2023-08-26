Feature Importance Analysis:
To determine which bacteria have the greatest impact on wether a Coral is healthy or diseased we preformed a feature importance study.

I began by selecting the best machine learning algorithm using pycaret.
![plot](/NewData/PycaretBM.png)
After testing each ML Algorithm over 10 folds Logistic Regression Preformed the best in all categories besides recall.
Next I used pycarets hyperparameter tuning to find the best preforming parameters for our model.
![plot](/NewData/LR_Tune.png)
Now that we have our best model and best paramaters for it I created an LR model in SciKit which preformed as shown:
Accuracy: 0.9919354838709677
Classification Report:
              precision    recall  f1-score   support

           0       1.00      0.99      0.99        82
           1       0.98      1.00      0.99        42

    accuracy                           0.99       124
   macro avg       0.99      0.99      0.99       124
weighted avg       0.99      0.99      0.99       124
![plot](/NewData/LR_AUC.png)
![plot](/NewData/LRConMat.png)
As seen by the accuracy, auc, and confusion matrix our model predicts healthy vs diseased corals very well.
Now that we know the model works well we can begin the feature importance analysis.
I decided to use 3 different metrics of feature importance and to then compare the results to find the most prevalent features accross all 3 metrics.
I started with permutation importance, plotting the 20 most important features.
![plot](/NewData/LR_T20P_IMP.png)
I then looked at the Gini Importance.
![plot](/NewData/LR_T20G_IMP.png)
And finally the Shap Values
![plot](/NewData/LR_T20S_IMP.png)
