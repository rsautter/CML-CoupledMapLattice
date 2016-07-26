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


First Header | Second Header
------------ | -------------
Content from cell 1 | Content from cell 2
Content in the first column | Content in the second column


####Example
	python main.py -d -mat gaussian 256 -o -nit 20 -c 0.1 -map logistic


Outputs:


![mapExampleIt19](/cml/output/it19.png)

