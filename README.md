Deconvolution, multiexponential regression and model inference package for MS2 data

version : April 2021

Written by Ovidiu Radulescu, University of Montpellier
Copyright : This is published under 3-clause BSD
 
The package consists of the following main programs, to be executed in order:  

1) read_raw_long_movie_data.m reads long movie data from raw files. Needs the data files

2) read_short_movie_data.m reads short movie data from files. Needs the data files

3) Deconvol.m applies the deconvolution algorithm. It uses the variable DataExp and generates the Matlab variables 
DataPred, cPosPred and Fit, corresponding to predicted intensity data, predicted Pol II positions, and best fit 
values of the objective function, respectively. The genetic algorithm is implemented using the Matlab function ga. 
By default it uses 50 workers (50 cores). 
This step is computationally expensive. It generates the file result_name_cPosPred.mat. The next step can be tested
independently if the result file is already available. 

4) fit_unconstrained_2exp.m or fit_unconstrained_3exp.m perform the multi-exponential unconstrained regression fitting and compute 
the model parameters for two and 
three exponential and states, respectively. They take as input the Matlab variables DataExp, DataPred, cPosPred 
resulting from one or several movies corresponding to the same phenotype.  The nonparametric estimate of the survival 
function is performed using the Matlab function ecdf (Meyer-Kaplan method). The multiexponential regression of the survival 
function is performed using the Matlab function lsqnonlin. 
The output of fit_unconstrained_2exp.m is a set of parameters for the ON/OFF telegraph model. 
The output of fit_unconstrained_3exp.m is a set of parameters for two types of 3 states models, referenced here as M1 and M2.   

5) fit_constrained_3exp.m performs the 3-exponential regression fitting constrained by the relation 
A1 lambda1 + A2 lambda2 + A3 lambda3 = 0 corresponding to model M3

 

Copyright (c) 2020, Ovidiu Radulescu, University of Montpellier
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the University of Montpellier nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.