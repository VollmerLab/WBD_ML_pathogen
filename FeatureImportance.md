# Feature Importance Analysis of Disease Outcomes and Bacterial Associations on Staghorn Coral with White Band Disease

Aim: 
1. Use ML algorythims to classify healthy and diseased corals in the field using only their bacterial types (i.e. ASVs) and abundances, and
2. Identify which bacterial ASVs best explain the disease outcomes using feature importance values and SHAP plots.

## The Data
To determine which bacteria have the greatest impact on whether a Coral is healthy or diseased we performed a feature importance study.
The Dataset has bacteria prevalence data from field corals sampled in 2016 and 2017.
This is a summary of the data with each bacteria as an index and the quantity(log2 of the counts per million) of bacteria in the sample, the health, the year, the family, and the genus as the features.
![plot](/NewData/BacteriaDFSummary.png)
This is a summary of the data with the samples as an index and bacteria and health as the features. This is the dataset we will give to the AI.
![plot](/NewData/SampleDFSummary.png)
We sampled 270 healthy corals 143 diseased for a total of 413 individual corals. For each coral sample, we have 342 distinct bacteria(ASVs) with the prevalence for  
To prepare this data for the AI to use we perform a training test split on the data and set 30% of the data aside as the test set.
'''
features = dataset.drop(['health'], axis=1).columns  # Separate the health from the other features.
target = ['health']  # Health is the target which the AI is predicting.
y = dataset.loc[:, target]  # Dataset with the sample IDs as the index and the Health as the Feature
X = dataset.loc[:, features]  # Dataset with the sample IDs as the index and the Bacterial Prevalence as the Features


X_train, X_test, y_train, y_test = train_test_split(X,               # the input features.
                                                    y,               # the target(health).
                                                    test_size=0.3,   # set aside 30% of the data as the test set.
                                                    random_state=7)  # reproduce the results.
'''

## PyCaret ML Model Testing

PyCaret fits and compares the top ML models with our data to help select which models or sets of models perform best as a classifier.

I then fed this data into Pycaret to find the best machine-learning algorithm to use on this data.
![plot](/NewData/PycaretBM.png)

After testing each ML Algorithm over 10 folds Logistic Regression Performed the best in all categories besides recall.

Explain. figure Models are ranked based on 
Explain values are.
Based these results, LR, RF etc. all performed we and can be fined turned downstream. 

## Logistic Regression Model and Feature Importance.

Next, I used Pycaret's hyperparameter tuning to find the best-performing parameters for our model.
![plot](/NewData/LR_Tune.png)

## Model Setup
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

Explain matrix. 

## Feature Importance
I decided to use 3 different metrics of feature importance and then compare the results to find the most prevalent features across all 3 metrics.
I started with permutation importance, plotting the 20 most important features.
![plot](/NewData/LR_T20P_IMP.png)
I then looked at the Gini Importance.
![plot](/NewData/LR_T20G_IMP.png)
And finally the Shap Values
![plot](/NewData/LR_T20S_IMP.png)
