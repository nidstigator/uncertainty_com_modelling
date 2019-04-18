# Neural Circuit Mechanism of Uncertainty and Change-of-Mind

Decision-making is often accompanied by a degree of confidence on whether a choice is correct. Decision uncertainty, or lack in confidence, may lead to change-of-mind. Studies have identified the behavioural characteristics associated with decision confidence or change-of-mind, and their neural correlates. Although several theoretical accounts have been proposed, there is no neural model that can compute decision uncertainty and explain its effects on change-of-mind. We propose a neuronal circuit model that computes decision uncertainty while accounting for a variety of behavioural and neural data of decision confidence and change-of-mind, including testable model predictions. Our theoretical analysis suggests that change-of-mind occurs due to the presence of a transient uncertainty-induced choice-neutral stable-steady state and noisy fluctuation within the neuronal network. Our distributed network model indicates that the neural basis of change-of-mind is more distinctively identified in motor-based neurons. Overall, our model provides a framework that unifies decision confidence and change-of-mind.

## Getting Started

These instructions will get you a copy of the source code up and running on your local machine for the purpose of reproducing the results outlined in our main manuscript and supplementary information.

### Prerequisites


```
MATLAB 2018a+ (not tested on versions prior to this)
XPPAUT 8.0
Optional for figure exporting in MATLAB:
	export_fig (https://github.com/altmany/export_fig)
	Ghostscripts (you can install using homebrew on Mac OS/Unix)
	Xpdf (you can install using homebrew on Mac OS/Unix)
```

### Description of script files

```
atiya_etal_simulate.m: the main script file that contains source code for our neural circuit model and a simple response time perceptual decision making experiment (see Simulation and Analysis section in Methods)

atiya_etal_vandenberg2016.m: source code of neural circuit model with a paradigm motivated by van den Berg et al., (2016) to test effect of confidence on multi-stage decisions.

pubgraph.m: used to make graphs look nicer- enhances colours, fonts and background. Used to maintain a consistent style for all figues. (see copyright notice below)

plot_figure_x: used to reproduce figure X in manuscript/supplementary information, given that a simulation struct is loaded in memory prior to running the script file.

sensorimotor.ode: xpp source code corresponding for sensorimotor module.

convert_y_to_s: helper function to convert firing rates of sensorimotor populations (y_1, y_2) to NMDA-averaged synaptic gating variables (S1, S2) trajectories on the phase-plane.

contourfcmap-pkg-master: A better library to plot contour maps than the built-in MATLAB algorithm. See copyright notice below
```

## Basic steps to reproduce simulation results

1. Please add all files and folders under the root directory of this source code to path.
2. Run atiya_etal_simulate.m. Enter 100 for simulating under all evidence quality levels, and run for 8000 trials.
3. Please note that under a 3.6GHZ Core i7 machine with 32GB of RAM, the simulation takes about 1hr and requires 10GB of free RAM. 
4. Run plot_figure_x files. Please note that plot_figure_3 is designed to plot results from the other simulation file (atiya_etal_vandenberg2016.m). See Notes section below for special cases of plotting.
5. Please note that in cases of script files that plot multiple figures (e.g. plot_supp_7_8.m), you would need to change certain values/boolean flags (e.g.: is_motor_correct, is_motor_com and coh)
6. a) For the multi-stage paradigm, run the atiya_etal_vandenberg2016.m script. Enter 100 for all evidence quality levels, and run for 8000 trials. 
   b) Please note that in the case of the van den Berg et al. (2016) paradigm code, 2 trials are simulated for each trial.
   c) This means that it takes twice the amount of time (2hrs) on the same hardware specifications specified in point (3)
7. Finally, the atiya_etal_parameter_analysis can simulate the model under a given number of trials (prompted for) for an evidence quality (default is epsilon=6.4) under various parameter values. 
    a) Please note that this simulation is designed to test all permutations between given values for:
        i) The uncertainty timing parameters (top-down (dis)inhibition).
        ii) uncertainty-encoding population parameters, mu and j_mc0.
    b) this simulation takes about 24 hours to complete for 8000 trials per permutation (default sets of parameters found in script file).
    c) the resulting struct can be plotted using the plot_supp_4.m or plot_supp_11.m.

## Notes: 

- plot_supp_7_and_8 is useful for plotting single trials. You can change the flags in the that file to get any type of trial from the simulation struct. 

- To reproduce results Figure 1 in the paper, set t_inh (line 201) to 400ms. Then follow steps above.
- To reproduce results in Figures 2,4,5, run the atiya_etal_simulate.m with default values (t_inh = 500).
- To reproduce results in Figure 3, run atiya_etal_vandenberg2016 and follow steps above.
- To reproduce parameter analysis results (Supplementary Figure 4 and 11), run atiya_etal_parameter_analysis and see point 7 above. 
- To reproduce Supplementary Figure 5, atiya_etal_simulate.m with default values, except j_mc0: set j_mc0=0 in order to effectively kill the feedback loop. After that, use the plot_figure_2 file to plot the '<' pattern.

- Plotting Figure 5 and Supplementary Figure 9-10 requires the relevant nullclines data (exported from XPP) to be loaded into MATLAB. 

The analysis (plotting) code is tested with 8000 trials per condition with the above parameters. 

## XPP analysis
XPP code files can be found under folder xpp_code. For more info on running XPP code, see: http://www.math.pitt.edu/~bard/bardware/tut/xpptut.html

Once the sensorimotor.m code is simulated in XPP, the phase plane configurations presented in the paper can be verified.

- For evidence quality e=3.2 with no uncertainty input, set the coh parameter to 3.2 and mc to 0.
- For evidence quality e=25.6 with no uncertainty input, set the coh parameter to 25.6 and mc to 0.
- For evidence quality e=3.2 with uncertainty input, set the coh parameter to 3.2 and mc to any value between 0.003-0.007.

## Contributing

The code is not open for public contribution yet. This may change in the future.

## Authors

* **Nadim Atiya** - (https://github.com/nidstigator)

Main contributors to this research project: Inaki Rano, Girijesh Prasad and KongFatt Wong-Lin


## License

Copyright 2019 Nadim Atiya

uncertainty_com_modelling is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

uncertainty_com_modelling is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with uncertainty_com_modelling.  If not, see <https://www.gnu.org/licenses/>.


## External packages used
Contourfcmap- Plotting a better contour in Matlab. Copyright (c) 2015 Kelly Kearney
pubgraph - Make Pretty (Publishable) Graphs. Copyright (c) 2012 Ruth Livingstone


