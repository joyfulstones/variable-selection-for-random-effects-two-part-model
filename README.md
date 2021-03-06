# variable-selection-for-random-effects-two-part-model

We provide some SAS and Matlab code when the sample size n is 100 in this folder. Please change the current paths before you run the program.  

In the folder "SAS code", we provide the code for the full model method, the oracle method, the MIC method and the MIC method for bootstrap samples (bootstrap). For the full model method, the oracle method, and the MIC method, we only provide the code based on one replication (the data in the file "two_part21.xls"). For the MIC method for bootstrap samples, we provide the code based on 100 bootstrap samples. These bootstrap samples are generated using the data in the file "two_part21.xls". 
  
Since the calculation by using bootstrap method takes too long, we conduct our simulation on the server.

In the subdirectory \code\matlab code\generate data, we provide the code used to generate the data. You only need to run the program "maincsv.m", then the data and bootstrap samples will be generated in the subdirectory\code\matlab code\generate data. "data1.m" is a function called by "maincsv.m". Here we only provide the code that can generate one dataset and 100 bootstrap samples using this dataset.

In addition, we provide the sas code that can be used to obtain the estimation of parameters using the oracle method. You only need to run the program "oracle.sas" in the subdirectory \code\sas code\oracle, then the estimation of parameters are saved in a file named "estoracleall.csv" in the subdirectory \code\sas code\oracle\Results_save1.

Similarly, we provide sas code used to obtain the estimation of parameters for the full model method. The operation is under the subdirectory \code\sas code\full.

Finally, we provide the code of the MIC method. We first run SAS code "MIC.sas"to obtain the estimation of parameters, second derivative and log likelihood in the subdirectory \code\sas code\MIC. These result are saved in the files "estL0all.csv", "hesL0all.csv", "liL0all.csv"......."estL4all.csv", "hesL4all.csv", "liL4all.csv" (a total of five tuning parameters) in the subdirectory \code\sas code\MIC\Results_save1.

Then we run SAS code "bootstrap.sas" in the subdirectory \code\sas code\bootstrap to obtain the estimation of parameters, hessian matrix and log likelihood using the bootstrap samples. These result are saved in the files "estL0allr1.csv", "hesL0allr1.csv", "liL0allr1.csv"......."estL4allr1.csv", "hesL4allr1.csv", "liL4allr1.csv" in the subdirectory \code\sas code\bootstrap\Results_save1.

Next you should copy the files in the subdirectory \code\sas code\MIC\Results_save1 and the subdirectory \code\sas code\bootstrap\Results_save1 to the subdirectory \code\matlab code\MIC, then run the program "gcvsimulation.m" in the subdirectory \code\matlab code\MIC to obtain the estimate and confidence intervals of parameters. 

   


