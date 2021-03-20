## Contents
- [Overview](#overview)
- [Files](#files)
- [Usage](#usage)
- [Publication](#publication)


## Overview
Improving catalytic ability of enzymes is critical to the success of many metabolic engineering projects, but the search space of possible protein mutants is too large to explore exhaustively through experiments. Here, we demonstrate that an optimization algorithm based on a regression model we had developed can effectively design short peptide tags to improve solubility of a few model enzymes. Based on the protein sequence information, a support vector regression model we recently developed was used to evaluate protein solubility after small peptide tags were introduced to a target protein. The optimization algorithm guided the sequences of the tags to evolve towards variants that had higher solubility. The optimization results were validated successfully by measuring solubility and activity of the model enzyme with and without the identified tags. The solubility of one protein (tyrosine ammonia lyase) was more than doubled and its activity was improved by 250%. This strategy successfully increased solubility of another two enzymes (aldehyde dehydrogenase and 1-deoxy-D-xylulose-5-phosphate synthase) we tested. 


## Files
* SVM6.mat: SVM model with improved performance
* cleandata_1.csv: the whole dataset we use to train SVM and includes 3148 proteins, with amino acid composition and protein solubility
* fitness5.m and GA5.m: Optimization algorithm GA
* matlab_5.csv: amino acid compostion of 6 proteins from our lab
* new_sample_10.csv: amino acid composition of 10 proteins selected from the 3148 proteins for optimization. 

## Usage
* Run GA5.m to optimize protein sequence with higher protein solubility. The output is the adding number of each amino acid and their sum is 20 since we want to add peptide sequences with 20 amino acids in total. To convert it to sequence, we design their sequence manually, which means there are multiple types of solutions since our inputs only consider amino acid composition rather than sequence. To generate the sequence of the two tags to be added to a protein from the number of amino acids we minimized the occurrence of amino acid repeats, which reduced the difficulty in synthesizing the related DNA oligonucleotides.

## Publication
Han, X., Ning, W., Ma, X., Wang, X., & Zhou, K. (2020). Improving protein solubility and activity by introducing small peptide tags designed with machine learning models. Metabolic engineering communications, 11, e00138.




