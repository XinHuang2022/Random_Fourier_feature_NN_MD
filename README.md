# Implementation of random Fourier feature neural network approximation of potential field in molecular dynamics

## Purpose

This repository contains the **MATLAB** code for the manuscript titled [*Convergence Rates for Random Feature Neural Network Approximation in Molecular Dynamics*](https://arxiv.org/abs/2406.14791) with an application of Fourier feature network representation $\bar{v}_r(x)$ to approximate the target potential function $v(x)$ in molecular dynamics simulation. 

The training of the network is based on a sampled training data set $`\big\{ \big(x_j,V(x_j), \nabla V(x_j)\big)\ |\ j=1,\ldots,J \big\}`$, with the aim to approximate the correlation observables by using the reconstructed potential $\bar{v}_r(x)$ instead of the computationaly costly true potential $V(x)$ to obtain the force fields for molecular dynamics.

The **complete workflow** is decomposed into **four pipelines**, namely: 
* **Data sampling pipeline**: For sampling of the training data set and the corresponding testing data set utilizing the overdamped Langevin dynamics, under the true potential function $V(x)$ with fixed inverse temperature $\beta$.
* **Training pipeline**: For the optimization of the Fourier feature network with updates on the trainable frequency parameters $`\{\omega_k\}_{k=1}^K`$ and the corresponding amplitude coefficients $`\{\eta_k\}_{k=1}^K,`$ in order to minimize the regularized loss function

  $$\mathcal{L}_R(\omega,\eta)=\frac{1}{J} \sum _{j=1}^J \big(\alpha_1|v(x_j)-\bar v(x_j;\omega,\eta)|^2 + \alpha_2|\nabla v(x_j)-\nabla \bar v(x_j;\omega,\eta)|^2\big) +\lambda_1\sum _{k=1}^K|\eta_k|^2  +\lambda_2\big(\sum _{k=1}^K|\eta_k|^2\big)^2+\lambda_3\max\big(\sum _{k=1}^K|\eta_k|-C'',0\big). $$

* **Inference pipeline**: For the numerical evaluation of correlation observables at varying correlation time $\tau_i=i \Delta \tau,$ $i=0,1,2,\dots,$ with $\Delta \tau=0.1$. To approximate the integral in the phase space, we apply the Monte Carlo integration based on sample points generated by overdamped Langevin dynamics under the real part of the trained potential $\bar{v}_r(x):=\textup{Real}(\bar{v}(x))$.
*  **Visualization pipeline**: For plotting the generalization error of the Fourier feature network approximation, and visualizing the corresponding results for approximation of correlation observables.
  
## Getting Started

The main function can be accessed through the file named `main.m`, which contains four sections corresponding to the four pipelines introduced above. 
Some numerical tests of interest by varying the parameters could be:
* Changing the largest integer in the array `K_values`, which presents different performances of the Fourier feature network in reconstructing the potential function and in approximating the correlation function.
* Varying the training data set size parameter `J`, e.g., by contrasting the training results obtained with $J=10^4$ and $J=10^5$, we observe a decreased generalization error (testing loss) for the larger data set size $J$.
* Switching the paremeter `Mix_beta_true`. Setting `Mix_beta_true = 1` will implement the sampling under one fixed inverse temperature $\beta=1$, while setting `Mix_beta_true = 2` will implement a hybrid sampling technique under two inverse temperatures $\beta_1=1$ and $\beta_2=0.3$, each contributing half of the more comprehensive training data set.
* Increasing (or decreasing) the sample size parameter `M_MC` for the Monte Carlo integral in the phase space to approximately evaluate the correlation function will give smaller (or larger) statistical uncertainty in the numerical results.


## Usage Reference
The specification of the parameters in each pipeline are summarized as following:

**Data sampling pipeline**
*  `dim`:                                   the dimensionality of the space. For the convenience of visualization, our code assumes dimensionality 2. That is, $(x,p)\in\mathbb{R}^2\times \mathbb{R}^2$.
*  `alpha` and `gamma`:                     the parameters in the target potential function $V(x)$.
*  `R_a`, `R_b`, and `R_c`:                 the radius parameters of the smoothed cutoff functions.
*  `sigma_a`, `sigma_b`, and `sigma_c`:     the acuteness parameters of the smoothed cutoff functions.
*  `h_x_sample`:                            the time step size for the overdamped Langevin dynamics sampling.
*  `J`:                                     the size of the training data set and the testing data set.
*  `Num_replica`:                           the number of independent replicas of the training procedure, in order to provide an evaluation of the statistical uncertainty of the results.
*  `Num_parallel_worker`:                   the number of available CPU cores for the training process and the further computation of correlation observables.
*  `Mix_beta_true`:                         the parameter for flagging two different sampling strategies to generate the training data set. 

**Training pipeline**
*  `lambda_1`, `lambda_2`, and `lambda_3`:  the weight parameters of the penalization terms in the regularized loss function $\mathcal{L}_R(\omega, \eta)$.
*  `C_const_bound`:                         the constant parameter $C''$ in the third penalization term of the regularized loss function $\mathcal{L}_R(\omega, \eta)$.
*  `K_values`:                              an array containing a series of integer numbers, preferrably in the form of $2^M$, where $M=4,5,\dots,10$, denoting the number of nodes in the Fourier feature network.

**Inference pipeline**
* `M_MC`:                                   the sample size for the Monte Carlo integral method to evaluate the correlation function in the phase space.
* `tau`:                                    the largest correlation time in the numerical approximations of the correlation function.
* `Delta_tau_plot`:                         the $\Delta \tau$ parameter for the step length in the plot of the correlation function curve, with correlation time points $\tau_i=i \Delta \tau$, $i=0,1,\dots,\tau/\Delta\tau$.

**Visualization pipeline**

The visualization pipeline will not introduce any new parameters, but will require to load the data from an pre-stored array, which contains the *reference values* of the correlation function at each time point $\tau_i$ with high accuracy, so as to compute the $L^1$-error of the approximated correlation function curve.

The code for generating the reference correlation function curve is provided in the subdirectory `e_Reference_curve_generating`.


 
