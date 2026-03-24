# Assessing the potential of bee-collected pollen sequence data to train machine learning models for geolocation of sample origin

Repository for manuscript Assessing the potential of bee-collected pollen sequence data to train machine learning models for geolocation of sample origin by Rebecca A. Hayes, Andrew D. Kern, and Lauren C. Ponisio. 

This repository includes code and data products to reproduce analyses conducted for this manuscript. Note that for the PNW Forests project, latitude and longitude coordinates of samples are removed due to landowner privacy policies. These may be available upon written request and permission of landowners. 

The environment.yaml stores the python environment used to run this code.



The repository is organized by the training data used for each analysis: either **raw** representing the models trained using raw, unclassified RBCL sequence data, or **taxonomic** representing the models trained using taxonomically classified RBCL sequence data. 

Within each folder (**raw** or **taxonomic**), there are several files. 

First are the cleaned data for analysis (excluding PNW Survey due to privacy agreement)
 1. X_train.npy: The NumPy array containing the scaled and standardized training features (the raw or taxonomically classified RBCL sequences) with values representing the relative abundance of each feature in each sample
 2. y_train.npy: The NumPy array containing the scaled and standardized training targets (the latitude and longitude of sample origin) with values representing the latitude and longitude of each sample
 3. X_test.npy: The NumPy array containing the scaled and standardized test features (the raw or taxonomically classified RBCL sequences) with values representing the relative abundance of each feature in each sample
 4. y_test.npy: The NumPy array containing the scaled and standardized testing targets (the latitude and longitude of sample origin) with values representing the latitude and longitude of each sample
 5. X_cols_raw.npy: The column names for the raw training dataset, to be used in feature importance calculations
 6. X_cols_tax.npy: The column names for the taxonomically classified training dataset, to be used in feature importance calculations
 7. project_test.npy and project_train.npy: The project labels for the training and test data, to be used in visualizations and comparisons

These NumPy arrays will be loaded into the following notebooks for analysis
 1. runTaxonomicClassifiedModels.ipynb and runRawModels.ipynb: these notebooks run hyperparameter grid search and run all models using those hyperparameters, outputting model comparison summary tables
 2. visualizeTaxonomicClassifiedModels.ipynb and visualizeRawModels.ipynb: these notebooks create the main figures for the manuscript, including the sample-by-sample error comparison, real versus predicted scatterplot maps, and the feature importance bar graphs

Additionally, this repo includes the scripts to generate the supplementary figures for the manuscript. These are titled supplementRaw.ipynb and supplementTax.ipynb. Also, the radar chart figure was created using R in the script radar_chart.R . 

The src folder includes a script that houses any internal functions written for the project. 



