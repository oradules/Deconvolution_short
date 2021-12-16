%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% Copyright : This is published under 3-clause BSD
%%%%% Last change March 2021 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
How to use the code
1) create a folder code and put all the programs there
1) create a folder data and put all experimental files in that folder. these are excel files 
with some format constraints: i) data is on columns starting from E3 ending 
before CE700 (if these constraints are not true change Read_data.m) ii) pay attention
to the Windows language when using . or , for decimals in the excel file
2) run Read_data.m ; this will create two folders  images containing plots of data
and matdatafiles containing saved Matlab variables to be used for the next step
3) run Deconv.m ; this computes polymerase positions by using the GA. it generates
a folder matresultfiles containing Matlab variables to be used for the next step
this step is time consuming and should be better run in batch on a multicore machine
(at least 50 cores); on our meso@lr machines it takes 400 - 600 s per dataset
(for those used as example here) with 50 cores. execution time will be further
reduced in future versions by precomputing a matrix of polymerase signals for all
possible positions. 
4) run Fit2.m ; this estimates the parameters of a two exponential survival function.
Use the pooled option if you want to pool all the data files, otherwise
files are analysed separately. For this version pooling is all or none.  
It generates  a folder Results2 containing the results of the fit. It is medium in terms
of execution time  
5) run Fit3.m ; this estimates the parameters of a two exponential survival function.
Use the pooled option if you want to pool all the data files, otherwise
files are analysed separately. For this version pooling is all or none.  
It generates  a folder Results3 containing the results of the fit.

How to use the python version
1) create a folder 'code' and put all the scripts there
2) modify the path in the code line: sys.path.append('/home/mdouaihy/temporary/deconv_python/') (-- mention which file has this line.) to the path of the folder containing the folder 'code'
3) create a folder for data and put all experimental files in that folder.
These are excel files with some format constraints:
	i) data is on columns starting from E3 ending  before CE700 (if these constraints are not true change readDrosoData.py)
	ii) the experimental data file should be in an .xlsx format.
	iii) the intensity data corresponding to the same gene should be stored by either putting all of data files in the same folder or one excel sheet for the gene.
4) When you run the second part of the notebook it will output the numbers of threads on your PC. Be careful to use not more than half the total number of threads when running parrallel computing for deconvolution.
5) Assign the path corresponding to the folder containing your data to the variable 'DataFilePath'
6) Assign the path corresponding to the folder that you want your results to be stored in to the variable 'outputFolder'
7) change combined to correspond to the way you stored your excel files.
	If you have stored all of the excel files corresponding to the same gene in the folder then set combined equal to 1.
	else use combined=0.
  Make sure you change the parameters corresponding the gene (MS2 structure, time resolution...) in the file (...)
