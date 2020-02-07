# Sensitivity study to the neutrino mass hierarchy in the JUNO experiment

The aim of this study is to determine the sensitivity to the determination of the neutrino mass hierarchy in the JUNO experiment.

First, it is necessary to reproduce the oscillated antineutrino spectrum as seen by the JUNO experiment. This is done in the AntineutrinoSpectrum directory, see below for details.

Secondly, the sensitivity study will be performed through a chi-squared analysis. Codes will be uploaded once ready.

## How to run the codes

To run the codes, copy the entire content of the repository. 

Run either ```ipython main.py -i``` or ```python main.py``` on your shell.
In the former case, ```-i``` keeps the interactive shell open and a plot of the antineutrino spectrum with finite energy resolution should be returned (it may take a little since it is done via numerical convolution).
In the latter case, there will be no interactive shell and the plot will not be displayed. 

To select which plots to see and save, change from ```plot_this=False``` to ```plot_this=True```, or vice versa, where necessary in 'main.py' (see the file for more details).
All plots are saved in pdf format.

## Brief description of the content of the repository

* 'main.py': run this file to obtain plots and results

* 'latex.py': necessary to use LaTeX in the plots

### AntineutrinoSpectrum

* 'reactor.py' contains the class ```Reactor```.
This class has methods to create the reactor's flux, the IBD cross section and the reactor's spectrum.

* 'oscillation.py' contains the class ```Oscillation_Prob```.
This class has a method to obtain the survival probability both as a function of L/E and of E.

* 'spectrum.py' contains the class ```Oscillated_Spectrum``` which has ```Reactor``` and ```Oscillation_Prob``` as parent classes.
This class has a method to obtain the oscillated spectrum both with and without including the finite energy resolution.

* 'convolution.py' contains the class ```Convolution```.
This class has a method which executes the numerical convolution of a function f with a gaussian of fixed or variable width.
It is used in a method of the class ```Oscillated_Spectrum``` to introduce the finite experimental resolution in the oscillated spectrum.

* 'resolution.py' plots the energy resolution as a function of the deposited energy for given stochastic term, _a_, and constant term, _b_.

* '3d_graph.py' plots a 3-dimensional graphical representation of the numerical convolution of a function f with a gaussian of fixed or variable width.

Note: 'resolution.py' and '3d_graph.py' are independent of the main file, so they need to be run separately.

### IBD_events

The code to estimate the total number of IBD events in JUNO due to reactor antineutrinos given the number of year of data taking and the efficiencies can be found in this directory.

To run the code simply type either ```ipython IBD_main.py -i``` or ```python IBD_main.py``` on your shell.
The total number of events and the single contribution of the 12 reactors should be printed on shell.

In 'IBD_main.py', one can modify the number of year (```years```), the IBD detection efficiency (```epsilon_IBD```) and the time efficiency, i.d. the effective runtime of a year (```eff_T```).
