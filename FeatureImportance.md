# Feature Importance Analysis of Disease Outcomes and Bacterial Associations on Staghorn Coral with White Band Disease

Aim: 
1. Use ML algorythims to classify healthy and diseased corals in the field using only their bacterial types (i.e. ASVs) and abundances, and
2. Identify which bacterial ASVs best explain the disease outcomes using feature importance values and SHAP plots.

## The Data

16s rDNA amplicon sequencing of bacterial from 413 Staghorn corals (270 healthy; 143 diseased)
Total number of bacterial ASVs (i.e. species) is 342 ASVs - these represent bacteria found in 5% of samples.
Bacterial abundances are normalized as log2 counts per million to account for different sequence library sizes and make them more normaly distributed.

### Data Structure
 
Bacteria prevalence data from field corals sampled in 2016 and 2017.
Each bacteria as an index and the quantity(log2 of the counts per million) of bacteria in the sample, the health, the year, the family, and the genus as the features.

![plot](/NewData/BacteriaDFSummary.png)

### Data Summary

Data with the samples as an index and bacteria and health as the features. This is the dataset passed to the AI/ML models.

![plot](/NewData/SampleDFSummary.png)

## Test/Train Split

70/30 training test split

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

Each ML Algorithm tested with 10 fold cross validation

### Model Test Results

Table of model accuracy parameters ranked by accuracy. 

ML accuracies show how well the ML models predict the coral disease outcome (healthy or disease) based solely on their bacteria communities (which includes potential but currently unknown pathogens)

Logistic regression has the highest accuracy at 99.3%!!! with the top 5 modesl including different tree based classifiers including Random Forest also performing well.

![plot](/NewData/PycaretBM.png)

Accuracy is the prediction accuracy
AUC stands for Area under the Curve and is a standard performance metric for ML modeling.

## Logistic Regression Model and Feature Importance.

### LR Model
Hyperparameter tuning the best fit LR models with additional CV validation produces similar model results with 98.6% accuracy.

![plot](/NewData/LR_Tune.png)

![plot](/NewData/LR_AUC.png)

## LR Model Performance on Test Data

Accuracy: 0.992 # Can't get much better

Classification Report:
                         precision    recall  f1-score   support
           0 (Healthy)      1.00      0.99      0.99        82
           1 (Disease)      0.98      1.00      0.99        42

    accuracy                           0.99       124
   macro avg       0.99      0.99      0.99       124
weighted avg       0.99      0.99      0.99       124


### Confusion Matrix

Health state - Healthy = 0 and Disease = 1
X-axis is Predicted vs Y-axis is Observed Health State

The LR model correctly predicted health state in all but one case where a Healthy corals was predicted as disease, which biologically makes sense since some asymptomatic corals may have contrated disease but are not showing symptoms.

Implication! = we can use 16s bacterial sequencingn alone to accurately identify disease corals without any additional data and without knowing the actual pathogen. This will allow us to screen nursery raised corals for disease before being outplanted in wild.

![plot](/NewData/LRConMat.png)

## Using Feature Importance to Identify Key Bacteria (including Pathogens) Driving Disease Outcomes

Traditional statistics typically IDs hundreds of associations and potential pathogens hamstringing coral pathogen research.

We are using feature importance to identify the most important bacteria driving the models - starting with our best model. In the future, I plan to look at overlap between different ML models and look for consensus within the top 20 most important features.

### Top 20 - Gini Importance.

Gini importance is a standard importance plots and shows our top bacterial candidates.

ASV25 is a Francisellaceae bacterial strain that we think is a top pathogen from other work. 

![plot](/NewData/LR_T20G_IMP.png)

### Shap values

Shap plots are a more modern way of looking at the top features that we are incorporating in our work because they show the direction of importants (more or less on disease) and how they impact individual coral samples.

There is overlap in our top20 bacterial candidates, but SHAP plots like this show the direction. For example, ASV25 is strongly associated with disease as are the other top 3 ASVs each of which is a target pathogen for cultivation in the lab.

![plot](/NewData/LR_T20S_IMP.png)
