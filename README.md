# JUNO_codes
 
Sensitivity study to the neutrino mass hierarchy in the JUNO experiment.

To run the codes, copy the entire content of the repository except for the file 'test_main.py' and the directories 'prova' and 'FIT'.
Change the paths in the lines 2-5 of the file 'main.py'. For example: sys.path.insert(0,'your_path/JUNO_codes/Reactor').
Run 'ipython main.py -i' in your shell. 
A plot of the antineutrino spectrum with finite energy resolution should be returned (it may take a little since it is done via numerical convolution).
To obtain other plots, change from 'plot_this=False' to 'plot_this=True' where necessary in 'main.py' (see the file for more details).



Brief description of the content of this repository.

Reactor: 
'reactor.py' contains the class Reactor.
This class has methods to create the reactor's flux, the IBD cross section and the reactor's spectrum.

Oscillation: 
'oscillation.py' contains the class Oscillation_Prob.
This class has a method to obtain the survival probability both as a function of L/E and of E.

Spectrum: 
'spectrum.py' contains the class Oscillated_Spectrum which has Reactor and Oscillation_Prob as parent classes.
This class has a method to obtain the oscillated spectrum without including energy resolution.

Exp_resolution: 
- 'convolution.py' contains the class Convolution.
This class has a method which executes the numerical convolution of a function f with a gaussian of fixed or variable width.
It is used to introduce the finite experimental resolution in the oscillated spectrum.
(It will be included in a method in the class Oscillated_Spectrum)
- 'resolution.py' plots the energy resolution as a function of energy for given a and b
- '3d_graph.py' plots a 3-dimensional graphical representation of the convolution

'latex.py' is necessary to use LaTeX in the plots

Please ignore 'test_main.py' and the directory named 'prova'.
Please ignore also the directory named 'FIT' since the implemention is not ready yet.

