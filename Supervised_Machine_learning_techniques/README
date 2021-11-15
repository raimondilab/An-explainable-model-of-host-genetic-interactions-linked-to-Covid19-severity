This folder contain files and scripts used to achieved task nr. 2. 

1. In sub-task 2a, we used the "Feature_selection.ipynb" Jupyter notebook script to perform feature selection (removal of multicollinearity) to reduce 
the number of considered features initially screened through the Log-Odds-Ratio statistics (see task 1, sub-task 1b). We removed  multicollinearity
from features by considering variant allele counts with correlation coefficients (corr.≤|0.8|). 

2. In sub-task 2b, we performed Supervised Machine Learning with the Jupyter notebook script
"Supervised Machine Learning Approach_Focus on Stratified 5-Folds CV dataset 2000 cohort.ipynb". In this task we trained supervised learning models
for binary classification task by employing several algorithms, i.e. Support Vector Machine, Logistic Regression, Random Forest, and Extreme 
Gradient Boosting (XGBoost) classifiers.

3. sub-task 2c, we utilized same script (see sub-task 2b) and saved the Feature importance scores (weights) for each trained model (SVC, Logistic Regression, Random Forest, 
XGBoost) classifiers. The feature importance scores assigned from these models were aggregated across the 5-folds to generate
a non-zero panel of variants for further downstream analysis.

3a. We identified only 16 variants that consistently received non-zero coefficients in all decision tree-based models. We focused our attention on these 16 variants.
We re-evaluate the predictive power of these variants using the Jupyter notebook script "Supervised_Machine_Learning_focus_on_16_Full_supported_variants_fold_1_5.ipynb" 

3b. We saved ML models (SVC, Logsitic Regression, Random Forest and XGBoost) for the 16 full supported variants (see 3a.) "Supervised_ML_16_full_support_variants_saved_models.rar"

4. In sub-task 2d, we utilized the Jupyter notebook "Validation of Deployed 16 full supported variants ML models on new batch (3000 cohort).ipynb" to validate
best performing models trained using most supported variants with and without covariates on a followup cohort of sequenced, Italian patients. 

4a. the saved models' prediction probabilities aggregated across stratified 5-folds (see 4.) aggregated_median_prediction_probabilities_5_folds_full_support_covariates.csv

5. We saved feature generation matrix (see sub-task 3) for each stratified 5-folds as "dataset_feature_matrix.ra" used for the ML analyses (see 1, 2, 3a).

6. We saved Supervised Machine Learning models' (see sub-task 2) performance metrics across the stratified 5-folds as "fold_12345_performance_metrics_2000_cases.csv".
