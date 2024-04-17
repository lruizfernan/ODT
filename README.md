# Introduction to ODT

We present *ODT*, an R package to assign optimal treatment to each patient. 

## Installation
ODT can be installed from CRAN repository:

`install.packages("ODT")`


## Introduction

The package library has two main parts: 

* Training of a decision tree based ong genomic or butational data.
* Prediction of the optimal treatment for new patients whose genomic or mutational data is known. 


 *ODT* offers different functions depending on the desired action:

* **`trainTree`** if the user wants to train the model.
* **`predictTree`** if the user wants to predict optimal tretment with new data
* **`niceTree`** will be used when the user wants to obtain agraphical display of the tree. Unlike the `plot` function, in this case the users can save it as an image in a selected directory.

