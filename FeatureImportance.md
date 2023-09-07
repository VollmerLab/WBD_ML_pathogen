Feature Importance Analysis:
To determine which bacteria have the greatest impact on whether a Coral is healthy or diseased we performed a feature importance study.

The Dataset has bacteria prevalence data from field corals sampled in 2016 and 2017.

This is a summary of the data with each bacteria as an index and the quantity(log2 of the counts per million) of bacteria in the sample, the health, the year, the family, and the genus as the features.
![plot](/NewData/BacteriaDFSummary.png)

This is a summary of the data with the samples as an index and bacteria and health as the features. This is the dataset we will give to the AI.
![plot](/NewData/SampleDFSummary.png)


To prepare this data for the AI to use we perform a training test split on the data and set 30% of the data aside as the test set.
![plot](/NewData/TrainTestSplit.png)


I began by selecting the best machine learning algorithm using Pycaret.
![plot](/NewData/PycaretBM.png)
After testing each ML Algorithm over 10 folds Logistic Regression Performed the best in all categories besides recall.
Next, I used Pycaret's hyperparameter tuning to find the best-performing parameters for our model.
![plot](/NewData/LR_Tune.png)
Now that we have our best model and best parameters for it I created an LR model in SciKit which performed as shown:
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
As seen by the accuracy, auc, and confusion matrix our model predicts healthy vs. diseased corals very well.
Now that we know the model works well we can begin the feature importance analysis.
I decided to use 3 different metrics of feature importance and then compare the results to find the most prevalent features across all 3 metrics.
I started with permutation importance, plotting the 20 most important features.
![plot](/NewData/LR_T20P_IMP.png)
I then looked at the Gini Importance.
![plot](/NewData/LR_T20G_IMP.png)
And finally the Shap Values
![plot](/NewData/LR_T20S_IMP.png)
