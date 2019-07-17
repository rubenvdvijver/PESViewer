"""
This code reads in an input files containing the wells,
bimolecular products, transition states and
barrierless reactions and creates a PES plot
"""
from __future__ import print_function, division
import os
import sys
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pylab as plt
import matplotlib.image as mpimg
from PIL import Image
import numpy as np
import numpy.linalg as la
import math
# try import RDKit
try:
    import rdkit.Chem as Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import AllChem
except ImportError:
    pass
# end try

# try import pybel
try:
    import pybel
except ImportError:
    pass
# end try

# contains all the options for this PES
options = {}

# global parameters for the plot
xlow = 0.0  # lowest x value on x axis
xhigh = 0.0  # highest x value on x axis
xmargin = 0.0  # margin on x axis
ylow = 0.0  # lowest y value on x axis
yhigh = 0.0  # highest y value on x axis
ymargin = 0.0  # margin on y axis

wells = []  # list of wells
bimolecs = []  # list of bimolecular products
tss = []  # list of transition states
barrierlesss = []  # list of barrierless reactions

# text dictionary: key is the chemical structure, value is the text
textd = {}
# lines dictionary: key is the chemical structure,
# value is a list of lines for that structure
linesd = {}
# figures dictionary: key is the chemical structure,
# value is the figure of that structure
imgsd = {}
# extents of the images
extsd = {}


class dragimage(object):
    """
    Class to drag an image
    """
    def __init__(self, figure=None):
        if figure is None:
            figure = plt.gcf()
        # simple attibute to store the dragged text object
        self.struct = None
        self.img = None
        # Connect events and callbacks
        figure.canvas.mpl_connect("button_press_event", self.pick_image)
        figure.canvas.mpl_connect("button_release_event", self.release_image)
        figure.canvas.mpl_connect("motion_notify_event", self.move_image)
    # end def

    def pick_image(self, event):
        for key in imgsd.keys():
            self.struct = None
            self.img = None
            if (imgsd[key].get_extent()[0] < event.xdata and
                    imgsd[key].get_extent()[1] > event.xdata and
                    imgsd[key].get_extent()[2] < event.ydata and
                    imgsd[key].get_extent()[3] > event.ydata):
                self.struct = key
                self.img = imgsd[key]
                self.current_pos = (event.xdata, event.ydata)
                break
            # end if
        # end for
    # end def

    def move_image(self, event):
        if self.img is not None:
            old_extent = self.img.get_extent()
            xchange = event.xdata-self.current_pos[0]
            ychange = event.ydata-self.current_pos[1]
            extent_change = (xchange, xchange, ychange, ychange)
            extent = [old_extent[i] + extent_change[i] for i in range(0, 4)]
            self.img.set_extent(extent=extent)
            self.current_pos = (event.xdata, event.ydata)
            plt.draw()
    # end def

    def release_image(self, event):
        if self.img is not None:
            self.struct = None
            self.img = None
            save_im_extent()
            # end if
        # end if
    # end def
# end class


class selecthandler(object):
    """
    Class to select and move stationary points, which highlight its reactions
    and in the future (TODO) renders the 3D images
    """
    def __init__(self, figure=None):
        if figure is None:
            figure = plt.gcf()
        self.struct = None  # stationary point that is selected
        figure.canvas.mpl_connect("button_press_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event",
                                  self.on_release_event)
        figure.canvas.mpl_connect("motion_notify_event",
                                  self.motion_notify_event)
    # end def

    def on_pick_event(self, event):
        self.struct = None
        # TODO: find more efficient way to iterate
        # all stationary points in one loop?
        # create a new list?
        for w in wells:
            if self.is_close(w, event):
                self.struct = w
                highlight_structure(self.struct)
            # end if
        # end for
        for b in bimolecs:
            if self.is_close(b, event):
                self.struct = b
                highlight_structure(self.struct)
            # end if
        # end for
        for t in tss:
            if self.is_close(t, event):
                self.struct = t
                highlight_structure(self.struct)
            # end if
        # end for
        if self.struct is None:
            highlight_structure()
        # end if
    # end def

    def motion_notify_event(self, event):
        if self.struct is not None:
            # a stationary point got selected
            # current position of the stationary point
            old_pos = (self.struct.x, self.struct.y)
            self.struct.x = event.xdata  # set the new position
            # move all the elements(image, text and lines)
            updateplot(self.struct, (event.xdata-old_pos[0]))
            self.current_pos = old_pos
    # end def

    def on_release_event(self, event):
        if self.struct is not None:
            self.struct = None
            # save the x-values of the startionary points to a file
            save_x_values()
            # save the image extents (x and y coordinates) to a file
            save_im_extent()
        # end if
        return True
    # end def

    def is_close(self, struct, event):
        """
        An event is close if it comes within 2% of the stationary point
        both in the x and in the y direction
        """
        xc = math.fabs(event.xdata - struct.x) < (xhigh-xlow)*0.02
        yc = math.fabs(event.ydata - struct.y) < (yhigh-ylow)*0.02
        return xc and yc
# end class


class line:
    """
    A line contains information about the line on the graph
    it is either a line between a reactant and ts, between a ts and product
    or between a reactant and product (for barrierless reactions)
    the chemstructs are the reactant and ts (f orward),
    ts and product (reverse), or reactant and product (barrierless)
    """
    def __init__(self, x1, y1, x2, y2, chemstruct=None, col='black'):
        if x1 <= x2:
            self.xmin = x1
            self.y1 = y1  # y1 corresponds to xmin
            self.xmax = x2
            self.y2 = y2  # y2 corresponds to xmax
        else:
            self.xmin = x2
            self.y1 = y2  # y1 corresponds to xmin
            self.xmax = x1
            self.y2 = y1  # y2 corresponds to xmax
        # end if
        if x1 == x2 or y1 == y2:
            self.straight_line = True
            self.coeffs = []
        else:
            self.straight_line = False
            self.coeffs = get_polynomial(self.xmin, self.y1,
                                         self.xmax, self.y2)
        # end if
        if chemstruct is None:
            self.chemstruct = []
        else:
            self.chemstruct = chemstruct
        self.color = col
    # end def
# end class


class well:
    # Well class, contains the name, smiles and energy of a well
    def __init__(self, name, energy, smi=None):
        self.name = name
        self.energy = energy
        self.smi = smi
        self.x = 0.
        self.y = 0.
        self.xyz_files = []
        fn = 'xyz/{name}.xyz'.format(name=name)
        if os.path.exists(fn):
            self.xyz_files.append(fn)
    # end def
# end class


class bimolec:
    # Bimolec class, contains the name,
    # both smiles and energy of a bimolecular product
    def __init__(self, name, energy, smi=None):
        self.name = name
        self.energy = energy
        if smi is None:
            self.smi = []
        else:
            self.smi = smi
        self.x = 0.
        self.y = 0.
        self.xyz_files = []
        # this bimolecular product is placed on the right side of the graph
        self.right = False
        i = 1
        fn = 'xyz/{name}{i}.xyz'.format(name=name, i=i)
        while os.path.exists(fn):
            self.xyz_files.append(fn)
            i += 1
            fn = 'xyz/{name}{i}.xyz'.format(name=name, i=i)
        # end for
    # end def
# end class


class ts:
    """
    TS class, contains the name, the names of the
    reactant and product and the energy of the ts
    """
    def __init__(self, name, names, energy, col='black'):
        self.name = name
        self.energy = energy
        self.color = col
        self.xyz_files = []
        fn = 'xyz/{name}.xyz'.format(name=name)
        if os.path.exists(fn):
            self.xyz_files.append(fn)
        self.reactant = next((w for w in wells if w.name == names[0]), None)
        if self.reactant is None:
            list = (b for b in bimolecs if b.name == names[0])
            self.reactant = next(list, None)
        if self.reactant is None:
            e = exceptions.not_recognized('reactant', names[0], name)
            raise Exception(e)
        self.product = next((w for w in wells if w.name == names[1]), None)
        if self.product is None:
            list = (b for b in bimolecs if b.name == names[1])
            self.product = next(list, None)
        if self.product is None:
            e = exceptions.not_recognized('product', names[1], name)
            raise Exception(e)
        self.lines = []
        self.x = 0.
        self.y = 0.
    # end def
# end class


class barrierless:
    """
    Barrierless class, contains the name and the
    names of the reactant and product
    """
    def __init__(self, name, names, col='black'):
        self.name = name
        self.xyz_files = []
        self.color = col
        fn = 'xyz/{name}.xyz'.format(name=name)
        if os.path.exists(fn):
            self.xyz_files.append(fn)
        self.reactant = next((w for w in wells if w.name == names[0]), None)
        if self.reactant is None:
            list = (b for b in bimolecs if b.name == names[0])
            self.reactant = next(list, None)
        if self.reactant is None:
            e = exceptions.not_recognized('reactant', names[0], name)
            raise Exception(e)
        self.product = next((w for w in wells if w.name == names[1]), None)
        if self.product is None:
            list = (b for b in bimolecs if b.name == names[1])
            self.product = next(list, None)
        if self.product is None:
            e = exceptions.not_recognized('product', names[1], name)
            raise Exception(e)
        self.line = None
    # end def
# end class


class exceptions(object):
    """
    Class that stores the exception messages
    """
    def not_recognized(role, ts_name, species_name):
        s = 'Did not recognize {role}'.format(role=role)
        s += ' {prod} '.format(prod=species_name)
        s += 'for the transition state {ts}'.format(ts=ts_name)
        return s


def get_polynomial(x1, y1, x2, y2):
    """
    Method fits a third order polynomial through two points as such
    that the derivative in both points is zero
    This method should only be used if x1 is not equal to x2
    """
    if x1 == x2:
        print('Error, cannot fit a polynomial if x1 equals x2')
        sys.exit()
    else:
        y = np.matrix([[y1], [y2], [0], [0]])
        x = np.matrix([[x1**3, x1**2, x1, 1],
                      [x2**3, x2**2, x2, 1],
                      [3*x1**2, 2*x1, 1, 0],
                      [3*x2**2, 2*x2, 1, 0]])
        xinv = la.inv(x)
        a = np.dot(xinv, y)
        return np.transpose(a).tolist()[0]
# end def


def read_input(fname):
    """
    Method to read the input file
    """
    if not os.path.exists(fname):  # check if file exists
        raise Exception(fname + ' does not exist')
    # end if
    f = open(fname, 'r')
    a = f.read()
    f.close()
    a = a.replace('\r\n', '\n').replace('\r', '\n').replace('\n\n', '\n')
    inputs = get_sd_prop(a)
    options['id'] = inputs['id'][0]
    # by default, print the graph title
    options['title'] = 1
    # default units
    options['units'] = 'kJ/mol'
    # use xyz by default, put 0  to switch off
    options['use_xyz'] = 1
    # no rescale as default, put the well or
    # bimolecular name here to rescale to that value
    options['rescale'] = 0
    # default figure height
    options['fh'] = 9.
    # default figure width
    options['fw'] = 18.
    # default margin on the x and y axis
    options['margin'] = 0.2
    # default dpi of the molecule figures
    options['dpi'] = 120
    # does the plot need to be saved (1) or displayed (0)
    options['save'] = 0
    # booleans tell if the ts energy values should be written
    options['write_ts_values'] = 1
    # booleans tell if the well and bimolecular energy values should be written
    options['write_well_values'] = 1
    # color of the energy values for the bimolecular products
    options['bimol_color'] = 'red'
    # color of the energy values of the wells
    options['well_color'] = 'blue'
    # color or the energy of the ts, put to 'none' to use same color as line
    options['ts_color'] = 'none'
    # boolean tells whether the molecule images should be shown on the graph
    options['show_images'] = 1
    # boolean that specifies which code was used for the 2D depiction
    options['rdkit4depict'] = 1
    # font size of the axes
    options['axes_size'] = 10
    # font size of the energy values
    options['text_size'] = 10

    if 'options' in inputs:
        for line in inputs['options']:
            if line.startswith('title'):
                options['title'] = int(line.split()[1])
            elif line.startswith('units'):
                options['units'] = line.split()[1]
            elif line.startswith('use_xyz'):
                options['use_xyz'] = int(line.split()[1])
            elif line.startswith('rescale'):
                options['rescale'] = line.split()[1]
            elif line.startswith('fh'):
                options['fh'] = float(line.split()[1])
            elif line.startswith('fw'):
                options['fw'] = float(line.split()[1])
            elif line.startswith('margin'):
                options['margin'] = float(line.split()[1])
            elif line.startswith('dpi'):
                options['dpi'] = int(line.split()[1])
            elif line.startswith('save'):
                if not options['save_from_command_line']:
                    options['save'] = int(line.split()[1])
            elif line.startswith('write_ts_values'):
                options['write_ts_values'] = int(line.split()[1])
            elif line.startswith('write_well_values'):
                options['write_well_values'] = int(line.split()[1])
            elif line.startswith('bimol_color'):
                options['bimol_color'] = line.split()[1]
            elif line.startswith('well_color'):
                options['well_color'] = line.split()[1]
            elif line.startswith('ts_color'):
                options['ts_color'] = line.split()[1]
            elif line.startswith('show_images'):
                options['show_images'] = int(line.split()[1])
            elif line.startswith('rdkit4depict'):
                options['rdkit4depict'] = int(line.split()[1])
            elif line.startswith('axes_size'):
                options['axes_size'] = float(line.split()[1])
            elif line.startswith('text_size'):
                options['text_size'] = float(line.split()[1])
            elif line.startswith('#'):
                # comment line, don't do anything
                continue
            else:
                if len(line) > 0:
                    print('Cannot recognize input line:')
                    print(line)
                # end if
            # end if
        # end for
    else:
        print('Warning, the input file arcitecture has changed,' +
              'use an "options" input tag to put all the options')
    # end if

    for w in inputs['wells']:
        w = w.split()
        name = w[0]
        energy = eval(w[1])
        smi = None
        if len(w) > 2:
            smi = w[2]
        # end if
        w = well(name, energy, smi)
        wells.append(w)
    # end for
    for b in inputs['bimolec']:
        b = b.split()
        name = b[0]
        energy = eval(b[1])
        smi = []
        if len(b) > 2:
            smi = b[2:]
        b = bimolec(name, energy, smi)
        bimolecs.append(b)
    # end for

    # it is important that the wells and bimolecular products
    # are read prior to the transition states and barrierless
    # reactions because they need to be added as reactants
    # and products
    for t in inputs['ts']:
        t = t.split()
        name = t[0]
        energy = eval(t[1])
        names = [t[2], t[3]]
        col = 'black'
        if len(t) > 4:
            col = t[4]
        t = ts(name, names, energy, col=col)
        tss.append(t)
    # end for
    for b in inputs['barrierless']:
        b = b.split()
        name = b[0]
        names = [b[1], b[2]]
        col = 'black'
        if len(b) > 3:
            col = b[3]
        b = barrierless(name, names, col=col)
        barrierlesss.append(b)
    # end for
# end def


def get_sd_prop(all_lines):
    """
    The get_sd_prop method interprets the input file
    which has a structure comparable to sdf molecular files
    """
    # split all the lines according to the keywords
    inputs = all_lines.split('> <')
    # do not consider the first keyword, this contains the comments
    inputs = inputs[1:]
    ret = {}
    for inp in inputs:
        inp = inp .split('>', 1)
        kw = inp[0].lower()  # keyword in lower case
        val = inp[1].strip().split('\n')  # values of the keywords
        val = [vi.strip() for vi in val]
        val = [vi for vi in val if vi]
        ret[kw] = val
    # end for
    return ret
# end def


def position():
    """
    This method find initial position for all the wells, products and
    transition states. Initially, the wells are put in the middle and
    the products are divided on both sides. The transition states are
    positioned inbetween the reactants and products of the reaction
    """
    # do the rescaling, i.e. find the y values
    # y0 is the energy of the species to which rescaling is done
    y0 = next((w.energy for w in wells if w.name == options['rescale']), 0.)
    if y0 == 0.:
        list = (b.energy for b in bimolecs if b.name == options['rescale'])
        y0 = next(list, 0.)
    for w in wells:
        w.y = w.energy - y0
    # end for
    for b in bimolecs:
        b.y = b.energy - y0
    # end for
    for t in tss:
        t.y = t.energy - y0
    # end for

    # find the x values
    file_name = '{id}_xval.txt'.format(id=options['id'])
    if os.path.exists(file_name):
        # if file exists, read the x values
        fi = open(file_name, 'r')
        a = fi.read()
        fi.close()
        a = a.split('\n')
        for entry in a:
            if len(entry) > 0:
                name = entry.split(' ')[0]
                xval = eval(entry.split(' ')[1])
                for w in wells:
                    if w.name == name:
                        w.x = xval
                # end for
                for b in bimolecs:
                    if b.name == name:
                        b.x = xval
                # end for
                for t in tss:
                    if t.name == name:
                        t.x = xval
                # end for
            # end if
        # end for
    else:
        n = len(bimolecs)  # n is the number of bimolecular products
        n2 = n // 2

        for i, w in enumerate(wells):
            # wells are put in the middle
            w.x = n2 + 1 + i + 0.

        # end for
        for i, b in enumerate(bimolecs):
            # bimolecular products are put on both sides
            if i < n2:
                b.x = 1 + i + 0.
            else:
                b.x = len(wells) + 1 + i + 0.
                b.right = True
            # end if

        # end for
        for i, t in enumerate(tss):
            # transition states are put inbetween the reactant and proudct
            x1 = t.reactant.x
            x2 = t.product.x
            x3 = (x1+x2)/2.
            t.x = x3
        # end for
        save_x_values()  # write the x values to a file
    # end if
# end def


def generate_lines():
    """
    The method loops over the transition states and barrierless reactions
    and creates lines accordingly depending on the x and y coordinates
    """
    for t in tss:
        line1 = line(t.x,
                     t.y,
                     t.reactant.x,
                     t.reactant.y,
                     [t, t.reactant],
                     col=t.color)
        line2 = line(t.x,
                     t.y,
                     t.product.x,
                     t.product.y,
                     [t, t.product],
                     col=t.color)
        t.lines.append(line1)
        t.lines.append(line2)
    # end for
    for b in barrierlesss:
        b.line = line(b.reactant.x,
                      b.reactant.y,
                      b.product.x,
                      b.product.y,
                      [b, b.reactant, b.product],
                      col=b.color)
    # end for
# end def


def get_sizes():
    """
    Get the axis lengths and the sizes of the images
    """
    global xlow, xhigh, xmargin, ylow, yhigh, ymargin
    # TODO: what if wells or bimoleculs is empty,
    # min and max functions will give error
    if len(bimolecs) > 0:
        # x coords of tss are always inbetween stable species
        xlow = min(min([w.x for w in wells]), min([b.x for b in bimolecs]))
        xhigh = max(max([w.x for w in wells]), max([b.x for b in bimolecs]))
        # all tss lie above the lowest well
        ylow = min(min([w.y for w in wells]), min([b.y for b in bimolecs]))
    else:
        # x coords of tss are always inbetween stable species
        xlow = min([w.x for w in wells])
        xhigh = max([w.x for w in wells])
        # all tss lie above the lowest well
        ylow = min([w.y for w in wells])
    xmargin = options['margin']*(xhigh-xlow)
    try:
        yhigh = max([t.y for t in tss])
    except ValueError:
        yhigh = max([b.y for b in bimolecs])
    yhigh = max(yhigh, max([w.x for w in wells]))
    if len(bimolecs) > 0:
        yhigh = max(yhigh, max([b.y for b in bimolecs]))
    ymargin = options['margin']*(yhigh-ylow)
# end def


def plot():
    """
    Plotter method takes all the lines and plots them in one graph
    """
    global xlow, xhigh, xmargin, ylow, yhigh, ymargin

    def showimage(s):
        """
        Get the extent and show the image on the plot
        s is the structure of which the image is to be
        plotted
        """
        global ymargin
        fn = '{id}_2d/{name}_2d.png'.format(id=options['id'], name=s.name)
        if os.path.exists(fn):
            img = mpimg.imread(fn)
            extent = None
            if s.name in extsd:
                extent = extsd[s.name]
            else:
                if not options['rdkit4depict']:
                    options['dpi'] = 120
                # end if

                imy = len(img) + 0.
                imx = len(img[0]) + 0.
                imw = (xhigh-xlow+0.)/(options['fw']+0.)*imx/options['dpi']
                imh = (yhigh-ylow+0.)/(options['fh']+0.)*imy/options['dpi']

                if isinstance(s, bimolec):
                    if s.x > (xhigh-xlow)/2.:
                        extent = (s.x + xmargin/5.,
                                  s.x + xmargin/5. + imw,
                                  s.y-imh/2.,
                                  s.y+imh/2.)
                    else:
                        extent = (s.x - xmargin/5. - imw,
                                  s.x - xmargin/5.,
                                  s.y-imh/2.,
                                  s.y+imh/2.)
                    # end if
                else:
                    extent = (s.x - imw/2.,
                              s.x + imw/2.,
                              s.y-ymargin/5. - imh,
                              s.y-ymargin/5.)
                # end if
            # end if
            im = ax.imshow(img, aspect='auto', extent=extent, zorder=-1)
            # add to dictionary with the well it belongs to as key
            imgsd[s] = im
        # end if
    # end def
    lines = []
    for t in tss:
        lines.append(t.lines[0])
        lines.append(t.lines[1])
    # end for
    for b in barrierlesss:
        lines.append(b.line)
    # end for
    get_sizes()
    plt.rcParams["figure.figsize"] = [options['fw'], options['fh']]
    plt.rcParams["font.size"] = options['axes_size']

    matplotlib.rc("figure", facecolor="white")
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.set_xticklabels([])

    if options['show_images']:
        for w in wells:
            showimage(w)
        # end for
        for b in bimolecs:
            showimage(b)
        # end for
    save_im_extent()  # save the positions of the images to a file
    for i, line in enumerate(lines):
        lw = 1.5
        alpha = 1.0
        ls = 'solid'
        if line.color == 'gray':
            ls = 'dotted'
        elif line.color == 'blue' or line.color == 'b':
            ls = 'dashed'
        if line.straight_line:
            if line.xmin == line.xmax:  # plot a vertical line
                ymin = min(line.y1, line.y2)
                ymax = max(line.y1, line.y2)
                a = ax.vlines(x=line.xmin,
                              ymin=ymin,
                              ymax=ymax,
                              color=line.color,
                              ls=ls,
                              linewidth=lw,
                              alpha=alpha)
                for struct in line.chemstruct:
                    # add to the lines dictionary
                    linesd[struct].append(a)
                # end for
            else:  # plot a horizontal line
                a = ax.hlines(y=line.y1,
                              xmin=line.xmin,
                              xmax=line.xmax,
                              color=line.color,
                              linewidth=lw,
                              alpha=alpha)
                for struct in line.chemstruct:
                    # add to the lines dictionary
                    linesd[struct].append(a)
                # end for
            # end if
        else:
            xlist = np.arange(line.xmin,
                              line.xmax,
                              (line.xmax-line.xmin) / 1000)
            a = line.coeffs
            y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
            pl = ax.plot(xlist,
                         y,
                         color=line.color,
                         ls=ls,
                         linewidth=lw,
                         alpha=alpha)
            for struct in line.chemstruct:
                # add to the lines dictionary
                linesd[struct].append(pl)
            # end for
        # end if
    # end for
    ax.set_xlim([xlow-xmargin, xhigh+xmargin])
    ax.set_ylim([ylow-ymargin, yhigh+ymargin])

    # write the name and energies to the plot
    for w in wells:
        if options['write_well_values']:
            t = ax.text(w.x,
                        w.y-ymargin/10,
                        '{:.1f}'.format(w.y),
                        fontdict={'size': options['text_size']},
                        ha='center', va='top',
                        color=options['well_color'],
                        picker=True)
            # add to dictionary with the well it belongs to as key
            textd[w] = t
        # end if
    # end for
    for b in bimolecs:
        if options['write_well_values']:
            # write the text values below the line:
            t = ax.text(b.x,
                        b.y-ymargin/10,
                        '{:.1f}'.format(b.y),
                        fontdict={'size': options['text_size']},
                        ha='center',
                        va='top',
                        color=options['bimol_color'],
                        picker=True)
            # add to dictionary with the well it belongs to as key
            textd[b] = t
        # end if
    # end for
    for t in tss:
        if options['write_ts_values']:
            color = t.color
            if options['ts_color'] != 'none':
                color = options['ts_color']
            te = ax.text(t.x,
                         t.y+ymargin/30,
                         '{:.1f}'.format(t.y),
                         fontdict={'size': options['text_size']},
                         ha='center',
                         va='bottom',
                         color=color,
                         picker=True)
            # add to dictionary with the well it belongs to as key
            textd[t] = te
        # end if
    # end for

    sel = selecthandler()
    dr = dragimage()

    if options['title']:
        plt.title('Potential energy surface of {id}'.format(id=options['id']))
    plt.ylabel('Energy ({units})'.format(units=options['units']))
    if options['save']:
        plt.savefig('{id}_pes_plot.png'.format(id=options['id']),
                    bbox_inches='tight')
    else:
        plt.show()
# end def


def generate_2d_depiction():
    """
    2D depiction is generated (if not yet available) and
    stored in the directory join(input_id, '_2d')
    This is only done for the wells and bimolecular products,
    2D of tss need to be supplied by the user
    """
    def get_smis(m, smis, files):
        # name and path of png file
        if len(smis) > 0:
            return smis
        elif all([os.path.exists(f) for f in files]):
            smis = []
            for f in files:
                try:
                    # weird syntax to allow python2 and python3
                    # python2: obmol = pybel.readfile('xyz', f).next()
                    # python3: obmol = pybel.readfile('xyz', f).__next__()
                    obmol = [mol for mol in pybel.readfile('xyz', f)][-1]
                    smis.append(obmol.write("smi").split()[0])
                except NameError:
                    print('Could not generate smiles for {n}'.format(n=m.name))
                # end try
            # end for
            return smis
    # end def

    def generate_2d(m, smis):
        # name and path of png file
        png = '{id}_2d/{name}_2d.png'.format(id=options['id'], name=m.name)
        if not os.path.exists(png):
            if len(smis) > 0:
                smi = '.'.join(smis)
                try:
                    options['rdkit4depict'] = 0
                    obmol = pybel.readstring("smi", smi)
                    obmol.draw(show=False, filename=png)
                except NameError:
                    try:
                        mol = Chem.MolFromSmiles(smi)
                        AllChem.Compute2DCoords(mol)
                        cc = mol.GetConformer()
                        xx = []
                        yy = []
                        for i in range(cc.GetNumAtoms()):
                            pos = cc.GetAtomPosition(i)
                            xx.append(pos.x)
                            yy.append(pos.y)
                        # end for
                        sc = 50
                        dx = (max(xx)-min(xx))*sc+30
                        dy = (max(yy)-min(yy))*sc+30
                        if not isinstance(dx, int):
                            dx = int(dx)
                        if not isinstance(dy, int):
                            dy = int(dy)
                        Draw.MolToFile(mol, png, size=(dx, dy), kekulize=False)
                    except NameError:
                        print('Could not generate 2d for {n}'.format(n=m.name))
                        return
                    # end try
                # end try

                im = Image.open(png)
                im.load()
                ix_low = im.getbbox()[0]
                ix_high = im.getbbox()[2]
                iy_low = im.getbbox()[1]
                iy_high = im.getbbox()[3]

                ar = np.asarray(im)

                for i, row in enumerate(ar):
                    if not all([all([ci == 255 for ci in c]) for c in row]):
                        if i > 10:
                            iy_low = i-10
                        else:
                            iy_low = i
                        # end if
                        break
                    # end if
                # end for
                for j, col in enumerate(ar.swapaxes(0, 1)):
                    if not all([all([ci == 255 for ci in c]) for c in col]):
                        if j > 10:
                            ix_low = j-10
                        else:
                            ix_low = j
                        # end if
                        break
                    # end if
                # end for
                for k, revrow in reversed(list(enumerate(ar))):
                    if not all([all([ci == 255 for ci in c]) for c in revrow]):
                        if k < iy_high - 10:
                            iy_high = k+10
                        else:
                            iy_high = k
                        # end if
                        break
                    # end if
                # end for
                for l, revcol in reversed(list(enumerate(ar.swapaxes(0, 1)))):
                    if not all([all([ci == 255 for ci in c]) for c in revcol]):
                        if k < ix_high - 10:
                            ix_high = l+10
                        else:
                            ix_high = l
                        # end if
                        break
                    # end if
                # end for
                im.crop((ix_low, iy_low, ix_high, iy_high)).save(png)
            else:
                # (TODO) add warning messages
                print('Could not generate 2d for {name}'.format(name=m.name))
                return
            # end if
        # end if
    # end def
    # make the directory with the 2d depictions, if not yet available
    dir = options['id'] + '_2d'
    try:
        os.stat(dir)
    except:
        os.mkdir(dir)
    for w in wells:
        s = []
        if w.smi is not None:
            s.append(w.smi)
        generate_2d(w, get_smis(w, s, w.xyz_files))
    # end for
    for b in bimolecs:
        generate_2d(b, get_smis(b, b.smi, b.xyz_files))
    # end for
# end def


def updateplot(struct, x_change):
    """
    Update the plot after the drag event of a stationary point (struct),
    move all related objects by x_change in the x direction,
    and regenerate the corresponding lines
    """
    global xlow, xhigh, xmargin, ylow, yhigh, ymargin
    # set the new sizes of the figure
    get_sizes()
    plt.gca().set_xlim([xlow-xmargin, xhigh+xmargin])
    plt.gca().set_ylim([ylow-ymargin, yhigh+ymargin])
    # generate new coordinates for the images
    if struct in imgsd:
        old_extent = imgsd[struct].get_extent()
        extent_change = (x_change, x_change, 0, 0)
        extent = [old_extent[i] + extent_change[i] for i in range(0, 4)]
        imgsd[struct].set_extent(extent=extent)
    # end if
    # generate new coordinates for the text
    if struct in textd:
        old_pos = textd[struct].get_position()
        new_pos = (old_pos[0]+x_change, old_pos[1])
        textd[struct].set_position(new_pos)
    # generate new coordinates for the lines
    for t in tss:
        if (struct == t or struct == t.reactant or struct == t.product):
            t.lines[0] = line(t.x,
                              t.y,
                              t.reactant.x,
                              t.reactant.y,
                              [t, t.reactant],
                              col=t.color)
            t.lines[1] = line(t.x,
                              t.y,
                              t.product.x,
                              t.product.y,
                              [t, t.product],
                              col=t.color)
            for i in range(0, 2):
                li = t.lines[i]
                if li.straight_line:
                    print('straight line')
                else:
                    xlist = np.arange(li.xmin,
                                      li.xmax,
                                      (li.xmax-li.xmin) / 1000)
                    a = li.coeffs
                    y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
                    linesd[t][i][0].set_xdata(xlist)
                    linesd[t][i][0].set_ydata(y)
                # end if
            # end for
        # end if
    # end for
    for b in barrierlesss:
        if (struct == b.reactant or struct == b.product):
            b.line = line(b.reactant.x,
                          b.reactant.y,
                          b.product.x,
                          b.product.y,
                          [b.reactant, b.product],
                          col=b.color)
            li = b.line
            if li.straight_line:
                print('straight line')
            else:
                xlist = np.arange(li.xmin, li.xmax, (li.xmax-li.xmin) / 1000)
                a = li.coeffs
                y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
                linesd[b][0][0].set_xdata(xlist)
                linesd[b][0][0].set_ydata(y)
            # end if
        # end if
    # end for
    plt.draw()
# end def


def highlight_structure(struct=None):
    """
    for all the lines, text and structures that are not
    directly connected to struct, set alpha to 0.15
    """
    if struct is None:
        alpha = 1.
    else:
        alpha = 0.15
    # end if
    # get all the tss and barrierlesss with this struct
    highlight = []
    lines = []
    for t in tss:
        if t == struct or t.reactant == struct or t.product == struct:
            highlight.append(t)
            if t.reactant not in highlight:
                highlight.append(t.reactant)
            if t.product not in highlight:
                highlight.append(t.product)
            lines = lines + linesd[t]
        # end if
    # end for
    for b in barrierlesss:
        if b.reactant == struct or b.product == struct:
            highlight.append(b)
            if b.reactant not in highlight:
                highlight.append(b.reactant)
            if b.product not in highlight:
                highlight.append(b.product)
            lines = lines + linesd[b]
        # end if
    # end for
    for struct in linesd:
        for li in linesd[struct]:
            if li in lines:
                li[0].set_alpha(1.)
            else:
                li[0].set_alpha(alpha)
            # end if
        # end for
    # end for
    for struct in textd:
        if struct in highlight:
            textd[struct].set_alpha(1.)
        else:
            textd[struct].set_alpha(alpha)
        # end if
    # end for
    for struct in imgsd:
        if struct in highlight:
            imgsd[struct].set_alpha(1.)
        else:
            imgsd[struct].set_alpha(alpha)
        # end if
    # end for
    plt.draw()
# end def


def save_x_values():
    """
    save the x values of the stationary points to an external file
    """
    fi = open('{id}_xval.txt'.format(id=options['id']), 'w')
    if len(wells) > 0:
        lines = ['{n} {v:.2f}'.format(n=w.name, v=w.x) for w in wells]
        fi.write('\n'.join(lines)+'\n')
    if len(bimolecs) > 0:
        lines = ['{n} {v:.2f}'.format(n=b.name, v=b.x) for b in bimolecs]
        fi.write('\n'.join(lines)+'\n')
    if len(tss) > 0:
        lines = ['{n} {v:.2f}'.format(n=t.name, v=t.x) for t in tss]
        fi.write('\n'.join(lines))
    fi.close()
# end def


def save_im_extent():
    """
    save the x values of the stationary points to an external file
    """
    fi = open('{id}_im_extent.txt'.format(id=options['id']), 'w')
    for key in imgsd:
        e = imgsd[key].get_extent()
        vals = '{:.2f} {:.2f} {:.2f} {:.2f}'.format(e[0], e[1], e[2], e[3])
        fi.write('{name} {vals}\n'.format(name=key.name, vals=vals))
    fi.close()
# end def


def read_im_extent():
    """
    Read the extents of the images if they are present in a file_name
    """
    fname = '{id}_im_extent.txt'.format(id=options['id'])
    if os.path.exists(fname):
        fi = open(fname, 'r')
        a = fi.read()
        fi.close()
        a = a.split('\n')
        for entry in a:
            pieces = entry.split(' ')
            if len(pieces) == 5:
                extsd[pieces[0]] = [eval(pieces[i]) for i in range(1, 5)]
            # end if
        # end for
    # end if
# end def


def main(argv):
    """
    Main method to run the PESViewer
    """
    if len(argv) > 1:  # read the arguments
        fname = argv[1]
        options['save'] = 0
        options['save_from_command_line'] = 0
        if len(argv) > 2:
            # argument to specify whether plot needs to be saved or displayed
            if argv[2] == 'save':
                options['save'] = 1
                options['save_from_command_line'] = 1
    # end if
    read_input(fname)  # read the input file
    # initialize the dictionaries
    for w in wells:
        linesd[w] = []
    # end for
    for b in bimolecs:
        linesd[b] = []
    # end for
    for t in tss:
        linesd[t] = []
    # end for
    for b in barrierlesss:
        linesd[b] = []
    # end for
    read_im_extent()  # read the position of the images, if known
    position()  # find initial positions for all the species on the graph
    generate_lines()  # generate all the line
    # generate 2d depiction from the smiles or 3D structure,
    # store them in join(input_id, '_2d')
    generate_2d_depiction()
    plot()  # plot the graph
# end def


main(sys.argv)
