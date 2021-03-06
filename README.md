# An Explainable Model of Host Genetic Interactions Linked to Covid-19 Severity
## Introduction
This project employed a multifaceted computational strategy to identify the genetic factors contributing to increased risk of severe COVID-19 infection from a Whole Exome Sequencing (WES) dataset of a cohort of 2000 Italian patients. We used a stratified k-fold screening, to rank variants more associated with severity, with training of multiple supervised classifiers, to predict severity on the basis of screened features. More details are found on (DOI:10.21203/rs.3.rs-1062190/v1).

### Project Tasks were organized as follows:

1. Variant screening;

   1a. Stratified K-fold split of sample cohort into train and test sets;

   1b. Variant screening;
   
   1c. Feature matrix generation; 

2. Supervised Machine Learning Techniques;

   2a. Feature selection: Removal of Multicollinearity;
   
   2b. Supervised Binary Classification;
   
   2c. Feature importance scores;
   
   2d. Final testing on a follow-up cohort;

3. Unsupervised Machine Learning Technqiues;

   3a. Principal Component Analysis (PCA) and clustering

4. Variant functional analysis;

    4a.  Pathway Enrichment Analysis;

    4b.  Retrieving associations between variants and disease traits or phenotypes.

**Authors**: Anthony Onoja, Nicola Picchiotti, Chiara Fallerini, Margherita Baldassarri, Francesca Fava, Francesca Colombo, Francesca Chiaromonte, Alessandra Renieri, Simone Furini, Francesco Raimondi.

***Software and Packages:*** All the analyses were performed using a customized Python script,Jupyter notebook, with the following libraries: scipy 1.2.0, numpy 1.19.4, scikit-learn 0.23.2., statsmodels 0.11.0 and matplotlib 3.2.1

**Note**: The WES dataset (variant stat., and feature count matrices) that support the findings of this study are available from GEN-COVID Multicenter Study group (https://clinicaltrials.gov/ct2/show/NCT04549831) and Francesco Raimondi (Ph.D.) head of Bioinformatics group at Bio@SNS lab but restrictions apply to the availability of these data, which were used under license for the current study, and so are not publicly made available. Data are however available from the authors upon reasonable request and with permission of GEN-COVID Multicenter Study group.


**Date**: 10/11/2021
