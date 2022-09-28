# CTD-squared Pancancer Chemosensitivity: Singular Value Decomposition and Regression Trees Approach
Rumiah Kessel
 
## Introduction
The goal of this study was to develop a model that predicts drug sensitivity across 515 distinct cell lines using drug response data from eleven cell lines and gene viability data from 515 cell lines. Two models were tested: both used principal component analysis (PCA) to reduce dimensionality of the training data, but one utilized sparse linear regression (lasso) while the other used bagged regression trees. The decision was made to move forward with the regression tree algorithm for reasons cited below.

## Methods

### Pre-Processing
Data was pre-processed in two ways before being fed into the bagged regression tree: first, a single drug's response data was averaged across the eleven cell lines and then multiplied by the Achilles data. This was both to apply the data across 515 cell lines as well as a way of incorporating the Achilles data into the training to aid the later PCA. However, due to limited Cancer Therapeutics Response Portal (“CTRP”) response data, the training set was limited to 442 of those cell lines. Second, this new data set was fed to a singular value decomposition (SVD) PCA algorithm, and the top 1,000 principal components were selected, accounting for 99.9% of the variance in the training set. Finally, the CTRP AUC sensitivity data was linearly mapped to be between 0 and 1.

### Algorithms
Two algorithms were tested: lasso linear regression and bagged regression trees. After testing, it became apparent that the linear regression was significantly worse at identifying the outlier sensitive combinations even though it was accurate in its estimates for the less sensitive combinations. Additionally, the linear regression was far more prone to giving false positives which ultimately led to the decision to use the bagged regression trees model fed with the PCA data. This resulted in better performance for middle-of-the-road sensitivities and significantly better performance for the more sensitive combinations. 

## Conclusion
The best performance was observed when a combination of a nonlinear approach and PCA was applied. One challenge was a lack of computing power to perform a larger PCA on the final prediction data set which necessitated reliance on the PCA from training. However, the approach of both PCA and the regression trees provided sound predictions. Based on these results, nonlinear regression algorithms warrant further research especially in their effective use alongside dimensionality reduction algorithms such as PCA.

## References
* CTD-squared Pancancer Chemosensitivity DREAM Challenge (syn21763589)
* "Correlating chemical sensitivity and basal gene expression reveals mechanism of action" Rees et al., Nat Chem Biol, 12, 109-116 (2016)
* "Harnessing Connectivity in a Large-Scale Small-Molecule Sensitivity Dataset" Seashore-Ludlow et al., Cancer Discovery, 5, 1210-1223 (2015)
* "An Interactive Resource to Identify Cancer Genetic and Lineage Dependencies Targeted by Small Molecules" Basu, Bodycombe, Cheah, et al., Cell, 154, 1151-1161 (2013)
* "A community effort to assess and improve drug sensitivity prediction algorithms" Costello JC, Heiser LM, Georgii E, et al. Nat Biotechnol. 2014;32(12):1202-1212. doi:10.1038/nbt.2877

