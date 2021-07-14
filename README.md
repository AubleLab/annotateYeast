# annotateYeast
 
 ## Installation: 
 `devtools::install_github("AubleLab/annotateYeast")`
 
 ## Example use: 
 
 ``` 
 library(annotateYeast)
 data("testData")
 data("yeastGenes")
 
 distances = calcFeatureDist_aY(testData, yeastGenes)
```
