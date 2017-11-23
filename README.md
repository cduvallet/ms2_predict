# ms2_predict

Predict chemical class label from MS2 data using classifiers trained on
HMDB and MassBankd databases.

# Directory Structure

```
|-- License? Citation?
|-- README.md
|-- requirements.txt
|
|-- data
|    |-- raw data downloaded from HDMB and MassBank (don't commit files though!)
|    |-- cleaned up raw data, and maybe feature engineering outputs
|    |-- analysis results files
|         |-- things like classifier results, etc that will be used to make figures
|
|-- src
|    |-- data
|    |	  |-- scripts to download data
|    |	  |-- scripts to clean/process data
|    |    |-- maybe scripts to do feature engineering? or maybe in a separate directory
|    |
|    |-- analysis
|    |	  |-- scripts to build and evaluate various classifiers
|    |
|    |-- figures
|         |-- scripts to make figures from the outputs of things in ../analysis
|
|-- final
|    |-- figures
|    |    |-- save all of the figures here
|    |
|    |-- paper
|         |-- we can also move the paper to a separate git repo
|
|-- docs (host on readthedocs)
```

# AShvin test line