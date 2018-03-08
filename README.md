# kdeborder
Kernel density estimation correcting for border bias (R package)

Arthur Charpentier and Ewen Gallic

This package proposes some R codes to compute the kernel density estimates of two-dimensional data points, using an extension of Ripley's circumference method to correct for border bias.
The method is described in our article: Charpentier, A. & Gallic, E. (2015). Kernel density estimation based on Ripley’s correction. GeoInformatica, 1-22. [Springer](https://link.springer.com/article/10.1007/s10707-015-0232-z).


To install the package in R, from Github:
```R
# install.packages("devtools")
devtools::install_github("ripleyCorr/kdeborder")
```
