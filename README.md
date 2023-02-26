# PESViewer
Visualize a potential energy surface

This code was originally developed by [Ruben Van de Vijver](https://github.com/rubenvdvijver). 
Our fork has been developed and improved on the top of the original to maintain compatibility with KinBot at integrate new features.

The PESViewer is a code to depict and analyze a potential energy surface 
characterized by wells, bimolecular products, transition states and barrierless reactions.
Their energy is needed to plot the potential energy surface. 
Written values of the energies and 2D plots (structural formulas) of the wells and products can be added to the figure 

To run this code, you need python version 3.7, matplotlib, numpy and optionally OpenBabel or RDKit to create 2D plots

## INSTALL

Clone the project onto your machine and go to the PESViewer directory. Type:

    python setup.py build
    python setup.py install 

## INPUT

A text file containing the stationary points, energies, names identifiers and options, see the input.txt as example.

For wells, write each well on a separate line as follows (smiles are optional to depict the molecule, if no xyz's are available):

    name energy smiles

For bimolecular products (again smiles are optional):

    name energy smiles

For reactions (colors are optional, and if available the line of the reaction will be given that color):

    name energy reactant product color

For bimolecular products (again colors are optional)

    name reactant product color


The plotting options (to be written in the input file) are:


| option | default | Description |
| ------- | ------- | ------- |
| title | 0 | print a title (1) or not (0) |
| units | kcal/mol | energy units |
| use_xyz | 1 |use xyz, put 0  to switch off |
| rescale | 0 | no rescale , put the well or bimolecular name here to rescale to that value |
| fh | 9. | figure height |
| fw | 18. | figure width |
| margin | 0.2 | margin fraction on the x and y axis |
| dpi | 120 | dpi of the molecule figures |
| save | 0 | does the plot need to be saved (1) or displayed (0) |
| write_ts_values | 1 | booleans tell if the ts energy values should be written |
| write_well_values | 1 | booleans tell if the well and bimolecular energy values should be written |
| bimol_color | red | color of the energy values for the bimolecular products |
| well_color | blue | color of the energy values of the wells |
| ts_color | green | color or the energy values of the ts, put to 'none' to use same color as line |
| show_images | 1 | boolean tells whether the molecule images should be shown on the graph |
| rdkit4depict | 1 | boolean that specifies which code to use for the 2D depiction |
| axes_size | 10 | font size of the axes |
| text_size | 10 | font size of the energy values on the graph |
| linear_lines | 0 | plot polynomials (0) or linear lines (1) between de stationary points |
| graph_edge_color | None | color of graph edge, if set to 'energy', will be scaled accordingly

Optionally a folder xyz/ containing the xyz coordinates of the stationary points ($name.xyz)
(for bimolecular products, use several xyz coordinates files ($name$index.xyz) )


## RUN

With the input file input.inp, type:

    pesviewer input.inp

## OUTPUT

The output is a modifiable matplotlib figure.

2 modifications are possible
- modifing the x-position of a stationary point by a 'drag and drop' of the energy value
- modifing the position of 2D structure images by a 'drag and drop' of the image

