Prediction of metabolic fluxes from gene expression data with Huber penalty convex optimization function
Shao-Wu Zhang*, Wang-Long Gou, Yan Li

Key Laboratory of Information Fusion Technology of Ministry of Education, School of Automation, Northwestern Polytechnical University, Xi¡¯an,710072, China 
The Supplementary Datasets
--------------------------

1. Matlab scripts (Dataset S1):

   - main.m (the main function that calls HPCOF algorithm that convert  gene expression data into flux distribution.

   - task1_setProblem.m (to set datafiles to read and map expression data onto reactions) and the following sub-function files:
     task1a_AandB.m (to take minimum for 'A AND B' rule)
     task1b_AorB_max.m (to take maximum for 'A OR B' rule)
     task1c_AorB_sum.m (to take summation for 'A OR B' rule)

   - task2_HPCOF.m  (to perform simulations for flux prediction) 

   - task3_HPCOF.m" (to display the results)
     
2. Network model and data files used for simulations of E. coli metabolism :
   - Ecoli_iAF1260.mat (Genome-scale network model for E. coli metabolism) 
   - data_reference (gene expression data for the reference condition of E. coli) 
   - data_NOX (gene expression data for NOX E. coli) 
   - data_ATP(gene expression data for ATPase E. coli) 

3. Network model and data files used for simulations of S. cerevisiae metabolism (Dataset S3):
   - yeast_5.21_MCISB.mat (Genome-scale network model for S. cerevisiae metabolism)
   - Lee_xx_genedata (gene expression data for S. cerevisiae)
   - Lee_xx_fluxdata.txt (flux data for S. cerevisiae)
   - Lee_fluxname_model.txt (names of experimentally measured fluxes by Lee et al.)


Usage notes
-----------

Implentation of the MATLAB scripts requires the installation of the cplex package.

1. Unzip and locate all files under the same folder

2. Run MATLAB 

3. Open "main.m"

4. Select the problem number as guided therein 

5. Run "main.m" to perform the following three tasks:
	- Task 1: set problem 
	- Task 2: perform simulations
	- Task 3: Report results 


