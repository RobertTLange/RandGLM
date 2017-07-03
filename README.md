# Randomized Numerical Linear Algebra for GLMs with Big Datasets
## Author: Robert Tjarko Lange
## Institution: Barcelona Graduate School of Economics
## Supervisors: Prof. Ioannis Kosmidis (UCL) and Prof. Omiros Paspapiliopoulos (UPF)

This thesis explores the potential application of randomized numerical linear algebra (RandNLA) methods for the fast approximate estimation of generalized linear models (GLMs). The following files are submitted with the main thesis. In order to replicate the simulation results and the figures displayed within the thesis please execute the "pipeline_thesis.R" script in R. 

It contains the following main scripts:

* "F1\_visualization\_transform.R": Reproduces figure 1 in the thesis. 
	* Running time (on "c3.8xlarge" AWS machine): ca. 2 secs

* "F2\_lev\_score\_approx.R": Reproduces figure 2 in the thesis.
	* Running time (on "c3.8xlarge" AWS machine): ca. 12 secs

* "F3\_sim\_ls.R": Reproduces figure 3 in the thesis. Simulates the randomized least squares estimator for the three different data-generating processes. All simulations are parallelized.
	* Running time (on "c3.8xlarge" AWS machine): ca. 2 mins

* "F4\_sim\_glm.R": Reproduces figure 4 in the thesis. Simulates the randomized GLM estimator for the three different data-generating processes.
	* Running time (on "c3.8xlarge" AWS machine): ca. 2 hours

* "F5\_glm\_trace.R" : Reproduces figure 5 in the thesis. Displays how the trace of a GLM estimator evolves over the iterations and for different data-generating processes. 
	* Running time (on "c3.8xlarge" AWS machine): ca. 3 mins
	
--------------------------------------------

Furthermore, the following scripts contain the functions needed:
	
* "multiplot.R" : Function that aligns subplots nicely.

* "A\_Basic.R": Defines basic functions used throughout the code (Cholesky decomposition, etc.). Also includes functions for Hadamard matrix generation

* "B\_RandLS.R": Defines functions used to compute the Random Sampling LS estimator.

* "C\_RandGLM.R": Defines functions used to compute the Random Sampling GLM estimator. 

--------------------------------------------

The simulations are parallelized and the running time heavily depends on the specs of your machine. Especially, the GLM simulations might take some time since the estimator does not converge often (at least not for very strict stopping criterions).
I recommend running the simulations on a virtual machine. This can be easily done in the following steps:

1) Launch an AWS EC2 instance. The bigger the machine, the faster the termination of the simulations. Make sure that it has a running R distribution installed. For my simulations I used the largest (non-GPU) available machine on AWS ("c3.8xlarge" - 32 CPU).

2) Go to your terminal and ssh into the instance.

	Shell: ssh -i "keyname.pem" ubuntu@Public_DNS

3) Open a new terminal window and scp the above scripts into the instance

	Shell: scp -i "keyname.pem" -r CODE_THESIS ubuntu@Public_DNS:Rjob
	
4) Execute the pipeline R script from the shell and obtain log of execution

	Shell: Rscript pipeline_thesis.R > sim_log.txt

5) Open new terminal window and scp all files back to your local machine in order to obtain all figures and simulation results.

	Shell: scp -i "keyname.pem" -r ubuntu@Public_DNS:Rjob .
	
June 2017,
	
Robert Lange
