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