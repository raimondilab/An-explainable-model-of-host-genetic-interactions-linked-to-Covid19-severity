The scripts in this folder focused on task nr. 1

1. In sub-task 1a, we utilized the python script "Stratified_K-fold split_of_sample_cohort_into_train_and_test sets.py". We embedded a strategy for variant screening in, 
stratified k-fold cross-validation (using the Strati¦edKFold function from the scikit-learn library) scheme to generate 5 random splits, into a training and testing test,
of the original dataset. Each fold was constituted by an 80 % training set which was also employed for variant screening and a 20 % testing set. 

2. In sub-task nr. 1b, we first run the python script "Variant_screening_training_sets.py" to generate the Log-Odds Ratio (LOR) statistics. This performed 
a case-control association and screened the variants associated with either severe or asymptomatic patients in each of the training sets for each of the 
five folds generated (task nr. 1). We used the python script "FDR_adjustment_of_LOR_p_value.py" for Fisher's exact test, to make adjustment to the original 
LOR p-value column of the contingency table. The python script "filter_variant_stat_grading_improved_trainset_5_fold_CV.py" was used to filter significant 
variants (i.e., variants with p-values < 0.05 and LOR > 1 are enriched among severe, while those with LOR < -1 are enriched among asymptomatics). we recycled 
task nr. 2 four-times to account for stratified 2 - 5 folds. 

3. In sub-task nr. 1c, for each split, we generated a feature matrix using the python script "Feature_matrix_generation_5_fold_CV.py" for the training set by 
assigning the allele counts of each screened variant for each sample of the training: i.e. 0 for genotype 0/0, 1 for genotypes 1/0 or 0/1, 2 for genotype 1/1. 
The feature matrix for the test set was defined by considering only variants identified as significant after screening the training set of the corresponding 
split and by assigning the allele count of each sample of the test set. We also included as additional features age, which was normalized, and gender, which 
was binarized by setting males to 0 and females to 1. Severe  patients from group “3+4+5” were given the classification label “1”, the asymptomatic patients 
from group 0 were given the label “0”.

     3a. Output files (see 3) were saved as "filtered_variants_phenotype_stratified_5_foldCVs.rar"

4. The python script "genotype_extraction_and_merging_follow_up_cohort_data.py" was used to extract VCF files and merge with relevant CSV files of 24 
Chromosomes in follow-up cohort study of 3000 patients. This unified the WES data into desirable format for variants match extractions identified from 2000 
cohort analysis. 

     4a. The python script "Feature_matrix_generation_full_supported_variants_follow_up_cohort.py" was used to generate the feature matrix from (4). 
