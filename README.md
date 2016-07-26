# 2D - CML

A simple python model of 2D Coupled Map Lattice Model (reference: https://en.wikipedia.org/wiki/Coupled_map_lattice).

###Required python libraries:
  * Matplotlib
  * Numpy
  * Scipy

###Running first example:
	python main.py

###Syntax
python main.py commands


Command | What does
------------ | -------------
-d | display the last iteration of CML
-mat type matlen | initial condition type(gaussian or bessel), matlen(integer>0)
-map map | initial condition type(gaussian or bessel), matlen(integer>0)
-nit nit | number of iterations
-c coupling | coupling factor(float)
-o | if set  saves each iteration png files
-h or --help | display help


####Example
	python main.py -d -mat gaussian 256 -o -nit 20 -c 0.1 -map logistic


Outputs:


![mapExampleIt19](/cml/output/it19.png)

