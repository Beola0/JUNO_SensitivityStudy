# Sensitivity study to the neutrino mass hierarchy in the JUNO experiment

This repository contains the codes to reproduce the antineutrino spectrum as seen by the JUNO experiment.

## How to run the codes

To run the codes, copy the entire content of the repository. 

Run either ```ipython main.py -i``` or ```python main.py``` on your shell.
In the former case, ```-i``` keeps the interactive shell open and a plot of the antineutrino spectrum with finite energy resolution should be returned (it may take a little since it is done via numerical convolution).
In the latter case, there will be no interactive shell and the plot will not be displayed. 

To select which plots to see, change from ```plot_this=False``` to ```plot_this=True```, or vice versa, where necessary in 'main.py' (see the file for more details).
All plots are saved in pdf format.

## Brief description of the content of the repository

* Reactor: 
'reactor.py' contains the class ```Reactor```.
This class has methods to create the reactor's flux, the IBD cross section and the reactor's spectrum.

* Oscillation: 
'oscillation.py' contains the class ```Oscillation_Prob```.
This class has a method to obtain the survival probability both as a function of L/E and of E.

* Spectrum: 
'spectrum.py' contains the class ```Oscillated_Spectrum``` which has ```Reactor``` and ```Oscillation_Prob``` as parent classes.
This class has a method to obtain the oscillated spectrum without including energy resolution.

* Exp_resolution: 
  * 'convolution.py' contains the class ```Convolution```.
This class has a method which executes the numerical convolution of a function f with a gaussian of fixed or variable width.
It is used in a method of the class ```Oscillated_Spectrum``` to introduce the finite experimental resolution in the oscillated spectrum
  * 'resolution.py' plots the energy resolution as a function of energy for given stochastic term, _a_, and constant term, _b_
  * '3d_graph.py' plots a 3-dimensional graphical representation of the convolution

* 'latex.py' is necessary to use LaTeX in the plots
