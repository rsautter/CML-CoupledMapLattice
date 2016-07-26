# 2D - CML

A simple python model of 2D Coupled Map Lattice Model (reference: https://en.wikipedia.org/wiki/Coupled_map_lattice).

In curently version: 
* Von Neumann neighborhood 
* 3 maps
* 2 initial conditions

###Required python libraries:
  * Matplotlib
  * Numpy
  * Scipy

###Running first example:
	python main.py

###Syntax
python main.py commands


Command | Arguments description
------------ | -------------
-d | Display CML last iteration
-mat type matlen | initial condition type (gaussian or bessel), matlen(integer>0)
-map map | type of map (logistic, doubling, som)
-nit nit | number of iterations
-c coupling | coupling factor(float)
-o | if set, saves each iteration in png files
-h or --help | display help


####Example
	python main.py -d -mat gaussian 256 -o -nit 20 -c 0.1 -map logistic


Outputs:


![mapExampleIt19](/cml/output/it19.png)

