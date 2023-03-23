# PESViewer
Visualize a potential energy surface

This code was originally developed by [Ruben Van de Vijver](https://github.com/rubenvdvijver). 
Our fork has been developed and improved on the top of the original to maintain compatibility with KinBot at integrate new features.

The PESViewer is a code to depict and analyze a potential energy surface 
characterized by wells, bimolecular products, transition states and barrierless reactions.
Their energy is needed to plot the potential energy surface. 
Written values of the energies and 2D plots (structural formulas) of the wells and products can be added to the figure 

To run this code, you need python version 3.7, matplotlib, numpy and optionally OpenBabel or RDKit to create 2D plots.

While the code can be used as a standalone application, it is designed to work with KinBot, which automatically generates the input files for PESViewer.

## INSTALL

Clone the project onto your machine and go to the PESViewer directory. There you need these commands:

    python setup.py build
    python setup.py install 

## INPUT

A text file containing the stationary points, energies, names identifiers and options, see the input.txt as example.

For wells, write each well on a separate line as follows (smiles are optional to depict the molecule, if no xyz's are available):

    name energy [smiles]

For bimolecular products:

    name energy [smiles]

For reactions (color is optional, and colors the pathway):

    name energy reactant product [color]

For barrierless reactions:

    name reactant product [color]


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

The other input is a folder called `xyz/` containing the xyz coordinates of the stationary points (`name.xyz`)
(for bimolecular products, use several xyz coordinates files, `name_index.xyz`). These 


## RUN

With the input file input.inp, type:

    pesviewer input.inp

## OUTPUT

When the code runs and the depiction of the molecules is requested, first, all depictions are generated based on the content of the `/xyz` folder. On repeated runs the depictions are reused, unless missing. 

**Traditional PES depiction** 

The output is a modifiable matplotlib figure, which can be displayed and arranged manually, or saved. 
The possible modifications are:
- modifing the x-position of a stationary point by draggning the energy value
- modifing the position of 2D structure images by dragging the image

Helper files with `.txt` extension are also generated, which, in certain cases need to be deleted to recreate the plot on subsequent runs (e.g., when changing the `fs` parameter). However, the information about the adjustmens made are stored here.

This is an example of a nicely arranged traditional PES plot from our recent paper:

![image](https://user-images.githubusercontent.com/40675474/227331800-373cf4b7-5d17-4f7a-8347-06544badc5b8.png)

* Martí, C., Michelsen, H. A., Najm, H. N., Zádor, J.: _Comprehensive kinetics on the C7H7 potential energy surface under combustion conditions._ J. Phys. Chem. A, **2023**, 127, 1941–1959. https://pubs.acs.org/doi/full/10.1021/acs.jpca.2c08035  

**Interactive graph representation**

Here, wells (black circles) and bimolecular products (blue circles) are shown as nodes of a graph. Their energies are printed as well, and the name of the species can be read by hovering over their depiction. Transition states are shown as edges. Edges representing saddle points are black, barrierless reaction pathways are gray. The thickness of the edge is inversely proportional to the absolute height of the barrier, and the energy can be read when hovering over the edge. Optionally, the color of the edges can be controlled using the 

Our code, built on the pyvis package, generates a `html` file, which can be opened in a browser. Note the settings below the graph to help arrange the graph using intuitive physics analogies, which can also be turned off.

![image](https://user-images.githubusercontent.com/40675474/227336672-c7448207-fc3b-42c3-ad89-7da45f84e985.png)

