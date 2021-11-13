The scripts in this folder focused on task nr. 2 and nr. 3.

In task nr. 2 we employed a Log-Odds Ratio (LOR) statistics to perform case-control association and to screen variants 
associated with either severe or asymptomatic patients in each of the training sets for each of the five folds generated.

In task nr. 3 for each split, we generated a feature matrix for the training set by assigning the allele counts of each 
screened variant for each sample of the training: i.e. 0 for genotype 0/0, 1 for genotypes 1/0 or 0/1, 2 for genotype 1/1. 
The feature matrix for the test set was defined by considering only variants identified as significant after screening the 
training set of the corresponding split and by assigning the allele count of each sample of the test set. We also included as
additional features age, which was normalized, and gender, which was binarized by setting males to 0 and females to 1. Severe 
patients from group “3+4+5” were given the classification label “1”, the asymptomatic patients from group 0 were given the label “0”.
