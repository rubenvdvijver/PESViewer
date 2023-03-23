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

Clone the project onto your machine and go to the PESViewer directory. There you need these commands:

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


The plotting options (to be written in the input file) are listed here. 


| option | default | Description | Scope* |
| ------- | ------- | ------- | ------- |
| `title` | `0` | print a title (1) or not (0) on the PES | T |
| `units` | `kcal/mol` | energy units | TG | 
| `use_xyz` | `1` | use xyz, put 0 to switch off | TG | 
| `rescale` | `0` | no rescale, put the well or bimolecular name here to set that/those species as the zero of energy | TG |
| `fh` | `9.` | figure height | T | 
| `fw` | `18.` | figure width | T |
| `fs` | `1` | scale factor for 2D depictions - they need to be deleted first to be regenerated with a new scaling | T |
| `lw` | `1.5` | line width of potential energy lines | T |
| `margin` | `0.2` | margin fraction on the x and y axis | T |
| `dpi` | `120` | dpi of the molecule figures | TG |
| `save` | `0` | does the plot need to be saved (1) or displayed (0) | T |
| `plot` | `1` | display (1) or not (0) traditional PES | T |
| `write_ts_values` | `1` | booleans tell if the ts energy values should be written | T |
| `write_well_values` | `1` | booleans tell if the well and bimolecular energy values should be written | T |
| `bimol_color` | `red` | color of the energy values for the bimolecular products | T |
| `well_color` | `blue` | color of the energy values of the wells | T |
| `ts_color` | `green` | color or the energy values of the ts, put to 'none' to use same color as line | T |
| `show_images` | `1` | boolean tells whether the molecule images should be shown on the graph | T |
| `rdkit4depict` | `1` | boolean that specifies which code to use for the 2D depiction | TG |
| `axes_size` | `10` | font size of the axes | T |
| `text_size` | `10` | font size of the energy values on the graph | T |
| `linear_lines` | `0` | plot polynomials (0) or linear lines (1) between de stationary points | T |
| `interpolation` | `hanning` | image interpolation method | T |
| `graph_edge_color` | `black` | color of graph edge, if set to `energy`, color is scaled accordingly | G |
| `reso_2D` | `1` | enable (1) or disable (2) generation of 2D depictions for resonant structures | TG |

* This column shows whether the parameter impacts the traditional (T) PES depiction, the graph (G), or both.

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

