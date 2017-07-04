# PESViewer UNDER CONSTRUCTION
Visualize a potential energy surface

## Developer: Ruben Van de Vijver, Laboratory for Chemical Technology, Ghent University

The PESViewer is a code to depict and analyze a potential energy surface 
characterized by wells, bimolecular products, transition states and barrierless reactions.
Their energy is needed to plot the potential energy surface. 
Written values of the energies and 2D plots (structural formulas) of the wells and products can be added to the figure 

To run this code, you need python version 2.7, matplotlib, numpy and optionally OpenBabel or RDKit to create 2D plots

## INPUT

A text file containing the stationary points, energies, names identifiers and options

Optionally a folder xyz/ containing the xyz coordinates of the stationary points ($name.xyz)
(for bimolecular products, use several xyz coordinates files ($name$index.xyz) )


## OUTPUT

The output is a modifiable matplotlib figure.

2 modifications are possible
- modifing the x-position of a stationary point by a 'drag and drop' of the energy value
- modifing the position of 2D structure images by a 'drag and drop' of the image

## OPTIONS

The options in the input file include
-rescaling of the energies to a stationary point

## TESTING

UNDER CONSTRUCTION

