# Dynamic Functional Connectivity Analysis
## About the thesis
This is the GitHub repository of my Bachelor Thesis in Physics. 
In this work, we studied whether and how a particular generative causal model (called MINDy) can reproduce a scalar feature of the time series from functional magnetic resonance imaging, the global dynamic functional connectivity speed. Understanding the complexity of brain activity, even during the resting state, is an open challenge in neuroscience.

## About the repository
The repository is organized into different files:
* <code>functions.py</code>: definitions of all the main functions used in the analysis.
* <code>notebook.ipynb</code>: Jupyter Notebook (10.3 MB) with the experiments and the analyses. It's organized itself into sections and subsections to be more readable. The variables whose calculation is particularly long are saved into the <code>files</code> directory.
* <code>training.m</code>: MATLAB simple pipeline to call the fitting function and save the obtained parameters into the <code>MINDy_parameters</code> directory.
* <code>sigma_plots.jl</code> does few more plots made in Julia.
















