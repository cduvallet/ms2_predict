# ms2_predict

Currently there doesn't exist a method to map from Mass Spectrometer data to the chemical class of the substance that generated that data. This would be very useful for diagnosing diseases from, for example, blood samples. For our 6.867 final project we wrote a data pipeline and experimented with a few methods for this reverse mapping. We used data from the HMDB database.


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
