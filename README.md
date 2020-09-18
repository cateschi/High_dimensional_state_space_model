# A dynamic factor model approach to incorporate Big Data in state space models for official statistics

This set of files contains the R codes and dataset for the estimation of the high-dimensional state space model proposed in the paper "A dynamic factor model approach to incorporate Big Data in state space models for official statistics" by Caterina Schiavoni, Franz Palm, Stephan Smeekes and Jan van den Brakel. Each R script can be employed to estimate the different models discussed in the paper. Namely,
+ The script _LFS_model.R_ contains the function that estimates the Labour Force Model by using the Labour Force series as unique observed series.
+ The script _LFS_model_with_CC.R_ contains the function that estimates the Labour Force Model augmented with the auxiliary series of claimant counts.
+ The script _LFS_model_with_GT.R_ contains the function that estimates the Labour Force Model augmented with the auxiliary series of Google Trends.
+ The script _LFS_model_with_CC_and_GT.R_ contains the function that estimates the Labour Force Model augmented with both auxiliary series of claimant counts and Google Trends.

The file _GoogleTrends_CC.xlsx_ contains the data for the claimant counts and the Google Trends that are used in the analysis.

In case of queries you can contact the corresponding author at c[dot]schiavoni[at]maastrichtuniversity[dot]nl.
