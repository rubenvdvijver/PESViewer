"""
This code reads in an input files containing the wells, bimolecular products, transition states
and barrierless reactions and creates a PES plot
"""

import os,sys,string
import subprocess
import matplotlib.pyplot
from matplotlib import pylab as plt
from matplotlib.text import Text
import matplotlib.image as mpimg
import numpy as np
import numpy.linalg as la
try:
    #sys.path.insert(0,os.environ['CHEMINFO'])
    import rdkit.Chem as Chem
    from rdkit.Chem import Draw
    from rdkit.Chem import AllChem
except ImportError:
    pass
#end try
import pybel
from PIL import Image

fname = 'input.txt'
id = 'pes_viewer_id'
units = ''
use_xyz = False
rescale = ''

fh = 9. # figure height d efault
fw = 18. # figure width d efault
xlow = 0.0 # lowest x value on x axis
xhigh = 0.0 # highest x value on x axis
xmargin = 0.0 # margin on x axis
ylow = 0.0 # lowest y value on x axis
yhigh = 0.0 # highest y value on x axis
ymargin = 0.0 # margin on y axis
imh = 0.0 # height of the images
imw = 0.0 # width of the images
    
wells = [] # list of wells
bimolecs = [] # list of bimolecular products
tss = [] #list of transition states
barrierlesss = [] # list of barrierless reactions

#text dictionary: key is the text of a chemical structure, value is the structure
textd = {}
#lines dictionary: key is the chemical structure, value is a list of lines for that structure
linesd = {}
#figures dictionary: key is the chemical structure, value is the figure of that structure
imgsd = {}
#extents of the images
extsd = {}

class dragimage(object):
    """
    Class to drag an image 
    """
    def __init__(self, figure=None) :
        if figure is None : figure = plt.gcf()
        # simple attibute to store the dragged text object
        self.struct = None
        self.img = None
        # Connect events and callbacks
        figure.canvas.mpl_connect("button_press_event", self.pick_image)
        figure.canvas.mpl_connect("button_release_event", self.release_image)
        figure.canvas.mpl_connect("motion_notify_event", self.move_image)
    #end def
    
    def pick_image(self,event):
        for key in imgsd.keys():
            self.struct=None
            self.img=None
            if (imgsd[key].get_extent()[0] < event.xdata and
            imgsd[key].get_extent()[1] > event.xdata and
            imgsd[key].get_extent()[2] < event.ydata and
            imgsd[key].get_extent()[3] > event.ydata):
                self.struct=key
                self.img=imgsd[key]
                self.current_pos = (event.xdata,event.ydata)
                break
            #end if
        #end for
        if self.struct != None:
            #xc =self.struct.x 
            #yc =self.struct.y
            #self.img.set_extent(extent=(xc - imw, xc + imw, yc-ymargin/20 - 2*imh, yc-ymargin/20))
            highlight_structure(self.struct)
        #end if
    #end def
    
    def move_image(self, event):
        if self.img is not None :
            old_extent = self.img.get_extent()
            xchange =event.xdata-self.current_pos[0]
            ychange =event.ydata-self.current_pos[1]
            extent_change = (xchange,xchange,ychange,ychange)
            extent =  [old_extent[i] + extent_change[i] for i in range(0,4)]
            self.img.set_extent(extent=extent)
            self.current_pos = (event.xdata,event.ydata)
            plt.draw()
    #end def
    
    def release_image(self,event):
        if self.img!= None:
            self.struct=None
            self.img=None
            restore()
            save_im_extent()
            #end if
        #end if
    #end def
    
#end class

class draghandler(object):
    """ 
    Class to handle Drag n Drop.
    This is a simple case, which works for Text objects only
    """
    def __init__(self, figure=None) :
        if figure is None : figure = plt.gcf()
        # simple attibute to store the dragged text object
        self.dragged = None
        # Connect events and callbacks
        figure.canvas.mpl_connect("pick_event", self.on_pick_event)
        figure.canvas.mpl_connect("button_release_event", self.on_release_event)
        figure.canvas.mpl_connect("motion_notify_event", self.motion_notify_event)
    #end def
    
    def on_pick_event(self, event):
        " Store which text object was picked and were the pick event occurs."
        if isinstance(event.artist, Text):
            self.dragged = event.artist
            self.current_pos = (event.mouseevent.xdata,event.mouseevent.ydata)
            highlight_structure(textd[self.dragged])
        #end if
        return True
    #end def
    
    def motion_notify_event(self, event):
        if self.dragged is not None :
            old_pos = self.dragged.get_position()
            new_pos = (event.xdata,old_pos[1])
            self.dragged.set_position(new_pos)
            updateplot(self.dragged,(event.xdata-old_pos[0]))
            self.current_pos = old_pos
    #end def
    
    def on_release_event(self, event):
        " Update text position and redraw"
        if self.dragged is not None :
            self.dragged = None
            restore()
            save_x_values()
            save_im_extent()
        #end if
        return True
    #end def

class line:
    """
    A line contains information about the line on the graph
    it is either a line between a reactant and ts, between a ts and product
    or between a reactant and product (for barrierless reactions)
    the chemstructs are the reactant and ts (f orward), 
    ts and product (reverse), or reactant and product (barrierless)
    """
    def __init__(self,x1,y1,x2,y2,chemstruct=[]):
        if x1 <= x2:
            self.xmin = x1
            self.y1 = y1 # y1 corresponds to xmin
            self.xmax = x2
            self.y2 = y2 # y2 corresponds to xmax
        else:
            self.xmin = x2
            self.y1 = y2  # y1 corresponds to xmin
            self.xmax = x1
            self.y2 = y1 # y2 corresponds to xmax
        #end if
        if x1 == x2 or y1 == y2:
            self.straight_line = True
            self.coeffs=[]
        else:
            self.straight_line = False
            self.coeffs=get_polynomial(self.xmin,self.y1,self.xmax,self.y2)
        #end if
        self.chemstruct = chemstruct
    #end def
#end class

class well:
    # Well class, contains the name, smiles and energy of a well
    def __init__(self,name,energy,smi=None):
        self.name = name
        self.energy = energy
        self.smi = smi
        self.x=0.
        self.y=0.
        self.xyz_file='xyz/%s.xyz'%name
    #end def
#end class

class bimolec:
    # Bimolec class, contains the name, both smiles and energy of a bimolecular product
    def __init__(self,name,energy,smi=[],nr=2):
        self.name = name
        self.energy = energy
        self.smi = smi
        self.x=0.
        self.y=0.
        self.xyz_files=[]
        self.right=False # this bimolecular product is placed on the right side of the graph
        for i in range(1,(nr+1)):
            self.xyz_files.append('xyz/%s%i.xyz'%(name,i))
        #end for
    #end def
#end class

class ts:
    """
    TS class, contains the name, the names of the reactant and product and the energy of the ts
    """
    def __init__(self,name,names,energy):
        self.name = name
        self.energy = energy
        self.reactant = next((w for w in wells if w.name == names[0]), None)
        if self.reactant == None: self.reactant = next((b for b in bimolecs if b.name == names[0]), None)
        if self.reactant == None: raise Exception('Did not recognize reactant %s for the transition state %s'%(names[0],name))
        self.product = next((w for w in wells if w.name == names[1]), None)
        if self.product == None: self.product = next((b for b in bimolecs if b.name == names[1]), None)
        if self.product == None: raise Exception('Did not recognize product %s for the transition state %s'%(names[1],name))
        self.lines = []
        self.x=0.
        self.y=0.
    #end def
#end class

class barrierless:
    """
    Barrierless class, contains the name and the names of the reactant and product 
    """
    def __init__(self,name,names):
        self.name = name
        self.reactant = next((w for w in wells if w.name == names[0]), None)
        if self.reactant == None: self.reactant = next((b for b in bimolecs if b.name == names[0]), None)
        if self.reactant == None: raise Exception('Did not recognize reactant %s for the transition state %s'%(names[0],name))
        self.product = next((w for w in wells if w.name == names[1]), None)
        if self.product == None: self.product = next((b for b in bimolecs if b.name == names[1]), None)
        if self.product == None: raise Exception('Did not recognize product %s for the transition state %s'%(names[1],name))
        self.line = None
    #end def
#end class

"""
Method fits a third order polynomial through two points as such
that the derivative in both points is zero
This method should only be used if x1 is not equal to x2
"""
def get_polynomial(x1,y1,x2,y2):
    y = np.matrix([[y1],[y2],[0],[0]])
    x = np.matrix([[x1**3 ,x1**2, x1, 1],[x2**3,x2**2,x2,1],[3*x1**2,2*x1,1,0],[3*x2**2,2*x2,1,0]])
    xinv = la.inv(x)
    a = np.dot(xinv,y)
    return np.transpose(a).tolist()[0]
#end def

def read_input():
    """
    Method to read the input file
    """
    global id, rescale, units, use_xyz
    if not os.path.exists(fname): # check if file exists
        raise Exception(fname + ' does not exist')
    #end if
    f = open(fname, 'r')
    a = f.read()
    f.close()
    a = a.replace('\r\n', '\n').replace('\r', '\n').replace('\n\n', '\n')
    inputs= get_sd_prop(a)
    id = inputs['id'][0]
    units = inputs['units'][0]
    if len(inputs['rescale']) > 0: rescale = inputs['rescale'][0]
    if len(inputs['use_xyz']) > 0: use_xyz = inputs['use_xyz'][0]=='yes'
    for w in inputs['wells']:
        w = w.split()
        name = w[0];
        energy = eval(w[1])
        smi = None
        if len(w) > 2:
            smi = w[2]
        #end if
        w = well(name,energy,smi)
        wells.append(w)
    #end for
    for b in inputs['bimolec']:
        b = b.split()
        name = b[0];
        energy = eval(b[1])
        smi = []
        if len(b) > 2:
            smi = b[2:]
        b = bimolec(name,energy,smi)
        bimolecs.append(b)
    #end for
    # it is important that the wells and bimolecular products are read prior to the
    # transition states and barrierless reactions because they need to be added
    # as reactants and products
    for t in inputs['ts']:
        t = t.split()
        name = t[0];
        energy = eval(t[1])
        names = [t[2],t[3]]
        t = ts(name,names,energy)
        tss.append(t)
    #end for
    for b in inputs['barrierless']:
        b = b.split()
        name = b[0];
        names = [b[1],b[2]]
        b = barrierless(name,names)
        barrierlesss.append(b)
    #end for
#end def

def get_sd_prop(all_lines):
    """
    The get_sd_prop method interprets the input file
    which has a structure comparable to sdf molecular files
    """
    inputs = all_lines.split('> <') # split all the lines according to the keywords
    inputs = inputs [1:] # do not consider the first keyword, this contains the comments
    ret = {}
    for inp in inputs:
        inp = inp .split('>', 1)
        kw = inp[0].lower() # keyword in lower case
        val = inp[1].strip().split('\n') # values of the keywords 
        val = map(string.strip, val)
        val = filter(None, val)
        ret[kw] = val
    #end for
    return ret
#end def

def position():
    """
    This method find initial position for all the wells, products and transition states
    Initially, the well are put in the middle and the products are divided on both sides
    The transition states are positioned inbetween the reactants and products of the reaction
    """
    #do the rescaling, i.e. find the y values
    # y0 is the energy of the species to which rescaling is done
    y0 = next((w.energy for w in wells if w.name == rescale), 0.)
    if y0 == 0.: y0 = next((b.energy for b in bimolecs if b.name == rescale), 0.)
    for w in wells:
        w.y = w.energy - y0
    #end for
    for b in bimolecs:
        b.y = b.energy - y0
    #end for
    for t in tss:
        t.y = t.energy - y0
    #end for
    # find the x values
    file_name = '%s_xval.txt'%id
    if os.path.exists(file_name): 
        # if file exists, read the x values
        fi = open(file_name,'r')
        a = fi.read()
        fi.close()
        a = a.split('\n')
        for entry in a:
            name = entry.split(' ')[0]
            xval = eval(entry.split(' ')[1])
            for w in wells:
                if w.name == name: w.x = xval
            #end for
            for b in bimolecs:
                if b.name == name: b.x = xval
            #end for
            for t in tss:
                if t.name == name: t.x = xval
            #end for
        #end for
    else:
        n = len(bimolecs) # n is the number of bimolecular products
        n2 = n/2
        
        for i,w in enumerate(wells):
            w.x = n2 + 1 + i + 0.
            
        #end for
        for i,b in enumerate(bimolecs):
            if i < n2:
                b.x = 1 + i + 0.
            else:
                b.x = len(wells) + 1 + i + 0.
                b.right = True
            #end if
            
        #end for
        for i,t in enumerate(tss):
            x1 = t.reactant.x
            x2 = t.product.x
            y1 = t.reactant.y
            y2 = t.product.y
            y3 = t.energy
            x3 = 0.
            if x1 == x2 or y1==y3 or y2==y3:
                x3 = (x1+x2)/2.
            else:
                # as such that the slopes of two straigt lines through the ts and both wells are inverse
                x3 = (((y3-y2)*x1+(y3-y1)*x2)/(2.*y3-y1-y2)) 
            #end if
            x3 = (x1+x2)/2.
            t.x = x3
            
        #end for
        save_x_values() # write the x values to a file
    #end if
#end def

def generate_lines():
    """
    The method loops over the transition states and barrierless reactions
    and creates lines accordingly depending on the x and y coordinates
    """
    for t in tss:
        line1 = line(t.x,t.y,t.reactant.x,t.reactant.y,[t,t.reactant])
        line2 = line(t.x,t.y,t.product.x,t.product.y,[t,t.product])
        t.lines.append(line1)
        t.lines.append(line2)
    #end for
    for b in barrierlesss:
        b.line = line(b.reactant.x,b.reactant.y,b.product.x,b.product.y,[b,b.reactant,b.product])
    #end for
#end def

def get_sizes():
    """
    Get the ax lengths and the sizes of the images
    """
    global fh,fw, xlow, xhigh,xmargin, ylow, yhigh, ymargin, imh, imw
    xlow = min(min([w.x for w in wells]),min([b.x for b in bimolecs])) # x coords of tss are always inbetween stable species
    xhigh = max(max([w.x for w in wells]),max([b.x for b in bimolecs])) 
    xmargin = 0.20*(xhigh-xlow) # take a 20% margin
    ylow = min(min([w.y for w in wells]),min([b.y for b in bimolecs])) #  all tss lie above the lowest well
    yhigh = max([t.y for t in tss]) #  all wells and bimolec products lie below the highest ts
    ymargin = 0.20*(yhigh-ylow) # take a 20% margin
    imh = (yhigh-ylow)/6. # height of the images
    imw = fh/fw*imh*(xhigh-xlow)/(yhigh-ylow) # width of the images
#end def

def plot():
    """
    Plotter method takes all the lines and plots them in one graph
    """
    global fh,fw, xlow, xhigh,xmargin, ylow, yhigh, ymargin, imh, imw
    def showimage(struct):
        global ymargin, imh, imw
        img=mpimg.imread('%s_2d/%s_2d.png'%(id,struct.name))
        extent = None
        if struct.name in extsd:
            extent = extsd[struct.name]
        else:
            if isinstance(struct,bimolec):
                if struct.x > (xhigh-xlow)/2.:
                    extent=(struct.x + xmargin/5 , struct.x +  xmargin/5 + imw, struct.y-imh/2, struct.y+imh/2)
                else:
                    extent=(struct.x - xmargin/5 - imw, struct.x -  xmargin/5, struct.y-imh/2, struct.y+imh/2)
                #end if
            else:
                extent=(struct.x - imw/2., struct.x + imw/2., struct.y-ymargin/20 - imh, struct.y-ymargin/20)
            #end if
        #end if  
        im=ax.imshow(img, aspect='auto', extent=extent, zorder=-1)
        imgsd[struct] = im # add to dictionary with the well it belongs to as key
    #end def
    lines = []
    for t in tss:
        lines.append(t.lines[0])
        lines.append(t.lines[1])
    #end for
    for b in barrierlesss:
        lines.append(b.line)
    #end for
    get_sizes()
    plt.rcParams["figure.figsize"]=[fw,fh]
    plt.rcParams["font.size"]=10
    matplotlib.rc("figure",facecolor="white")
    fig, ax = plt.subplots()
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('left')
    ax.set_xticklabels([])
    for w in wells:
        showimage(w)
    #end for
    for b in bimolecs:
        showimage(b)
    #end for
    save_im_extent() # save the positions of the images to a file
    for i,line in enumerate(lines):
        if line.straight_line:
            if line.xmin == line.xmax: # plot a vertical line
                ymin = min(line.y1,line.y2)
                ymax = max(line.y1,line.y2)
                a=ax.vlines(x = line.xmin,ymin=ymin,ymax=ymax,color='black')
                for struct in line.chemstruct:
                    # add to the lines dictionary
                    linesd[struct].append(a)
                #end for
            else: # plot a horizontal line
                a=ax.hlines(y = line.y1,xmin=line.xmin,xmax=line.xmax,color='black')
                for struct in line.chemstruct:
                    # add to the lines dictionary
                    linesd[struct].append(a)
                #end for
            #end if
        else:
            xlist = np.arange(line.xmin, line.xmax, (line.xmax-line.xmin)/1000)
            a = line.coeffs
            y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
            pl=ax.plot(xlist,y,color='black')
            for struct in line.chemstruct:
                # add to the lines dictionary
                linesd[struct].append(pl)
            #end for
        #end if
    #end for
    ax.set_xlim([xlow-xmargin,xhigh+xmargin])
    ax.set_ylim([ylow-ymargin,yhigh+ymargin])
    
    # write the name and energies to the plot
    for w in wells:
        t = ax.text(w.x,w.y-ymargin/30, '%.1f'%(w.y),ha='center',va='top',color='blue', picker=True)
        textd[t] = w # add to dictionary with the well it belongs to as key
    #end for
    for b in bimolecs:
        t=None
         # at this point, put the energy value below the actual energy, 
         # with this loop it is possible to put it next to it (on the correct side)
        if b.x > (xhigh-xlow)/2.:
            t=ax.text(b.x+xmargin/30,b.y, '%.1f'%(b.y),ha='left',va='center',color='red', picker=True)
        else: 
            t=ax.text(b.x-xmargin/30,b.y, '%.1f'%(b.y),ha='right',va='center',color='red', picker=True)
        #end if
        textd[t] = b # add to dictionary with the well it belongs to as key
    #end for
    for t in tss:
        te=ax.text(t.x,t.y+ymargin/30, '%.1f'%(t.y),ha='center',va='bottom',color='green',fontweight='bold', picker=True)
        textd[te]=t # add to dictionary with the well it belongs to as key
    #end for
    
    #plt.imshow(img)
    dragh = draghandler()
    dragi = dragimage()
    plt.title('Potential energy surface of %s'%id)
    plt.ylabel('Energy (%s)'%units)
    plt.show()
    #plt.savefig('plot.png')
#end def

def generate_2d_depiction():
    """
    2D depiction is generated (if not yet available) and stored in the directory join(input_id,'_2d')
    This is only done for the wells and bimolecular products, 2D of tss need to be supplied by the user
    """
    def generate_2d(m,smis,files):
        # name and path of png file
        png = '%s_2d/%s_2d.png'%(id,m.name)
        if not os.path.exists(png):
            if len(smis) > 0:
                smi = string.join(smis,'.')
                if module_exists(Chem):
                    mol = Chem.MolFromSmiles(smi)
                    AllChem.Compute2DCoords(mol)
                    Draw.MolToFile(mol,png,size=(120,120),fitImage=True)
                else:
                    obmol = pybel.readstring("smi",smi)
                    obmol.draw(show=False,filename=png)
                #end if
            elif all([os.path.exists(f) for f in files]):
                smis = []
                for f in files:
                    obmol = pybel.readfile(f)
                    smis.append[obmol.write]
                #end for
                if module_exists(Chem):
                    mol = Chem.MolFromSmiles(smi)
                    AllChem.Compute2DCoords(mol)
                    Draw.MolToFile(mol,png,size=(120,120),fitImage=True)
                else:
                    obmol = pybel.readstring("smi",smi)
                    obmol.draw(show=False,filename=png)
                #end if
            else:
                # (TODO) add warning messages
                print 'Could not generate 2d for %s'%m.name
            #end if
        #end if
    #end def
    # make the directory with the 2d depictions, if not yet available
    dir = id + '_2d'
    try:
        os.stat(dir)
    except:
        os.mkdir(dir)
    for w in wells:
        s = []
        if w.smi != None: s.append(w.smi)
        generate_2d(w,s,[w.xyz_file])
    #end for
    for b in bimolecs:
        generate_2d(b,b.smi,b.xyz_files)
    #end for
#end def

def updateplot(dragged,x_change):
    """
    Update the plot after the drag event of a text item (dragged), move all related objects
    by x_change in the x direction, and regenerate the corresponding lines
    """
    global fh,fw, xlow, xhigh,xmargin, ylow, yhigh, ymargin, imh, imw
    
    textd[dragged].x += x_change
    """
    if isinstance(textd[dragged],bimolec):
        # (TODO) add the text width to the positioning of the line and the image
        if dragged.get_position()[0] > (xhigh-xlow)/2:
            bbox = dragged.set_ha('left')
            textd[dragged].x = dragged.get_position()[0]-xmargin/30
        else:
            bbox = dragged.set_ha('right')
            textd[dragged].x = dragged.get_position()[0]+xmargin/30
        #end if
    else:
        textd[dragged].x = dragged.get_position()[0]
    #end if
    """
    # set the new sizes of the figure
    get_sizes()
    plt.gca().set_xlim([xlow-xmargin,xhigh+xmargin])
    plt.gca().set_ylim([ylow-ymargin,yhigh+ymargin])
    # generate new coordinates for the images
    if textd[dragged] in imgsd:
        struct = textd[dragged]
        old_extent = imgsd[struct].get_extent()
        extent_change = (x_change,x_change,0,0)
        extent =  [old_extent[i] + extent_change[i] for i in range(0,4)]
        imgsd[struct].set_extent(extent=extent)
        """
        if isinstance(struct,bimolec):
            if struct.x > (xhigh-xlow)/2.:
                imgsd[struct].set_extent(extent=(struct.x + xmargin/10 , struct.x +  xmargin/10 + imw, struct.y-imh/2, struct.y+imh/2))
            else:
                imgsd[struct].set_extent(extent=(struct.x - xmargin/10 - imw, struct.x -  xmargin/10, struct.y-imh/2, struct.y+imh/2))
            #end if
        else:
            imgsd[struct].set_extent(extent=(struct.x - imw/2., struct.x + imw/2., struct.y-ymargin/20 - imh, struct.y-ymargin/20))
        #end if
        """
    #end if
    # generate new coordinates for the lines   
    for t in tss:
        if (textd[dragged] == t or 
        textd[dragged] == t.reactant or 
        textd[dragged] == t.product):
            t.lines[0] = line(t.x,t.y,t.reactant.x,t.reactant.y,[t,t.reactant])
            t.lines[1] = line(t.x,t.y,t.product.x,t.product.y,[t,t.product])
            for i in range(0,2):
                li=t.lines[i]
                if li.straight_line:
                    print 'straight line'
                else:
                    xlist = np.arange(li.xmin, li.xmax, (li.xmax-li.xmin)/1000)
                    a = li.coeffs
                    y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
                    linesd[t][i][0].set_xdata(xlist)
                    linesd[t][i][0].set_ydata(y)
                #end if
            #end for
        #end if 
    #end for
    for b in barrierlesss:
        if (textd[dragged] == b.reactant or 
        textd[dragged] == b.product):
            b.line = line(b.reactant.x,b.reactant.y,b.product.x,b.product.y,[b.reactant,b.product])
            li=b.line
            if li.straight_line:
                print 'straight line'
            else:
                xlist = np.arange(li.xmin, li.xmax, (li.xmax-li.xmin)/1000)
                a = li.coeffs
                y = a[0]*xlist**3 + a[1]*xlist**2 + a[2]*xlist + a[3]
                linesd[b][0][0].set_xdata(xlist)
                linesd[b][0][0].set_ydata(y)
            #end if
        #end if
    #end for
    plt.draw()
#end def

def restore():
    """
    restore the highlighting, set the alpha of every element back to 1.
    """
    for struct in linesd:
        for li in linesd[struct]:
            li[0].set_color('black')
        #end for
    #end for
    for te in textd:
        te.set_alpha(1.)
    #end for
    for struct in imgsd:
        imgsd[struct].set_alpha(1.)
    #end for
    plt.draw()
#end def

def highlight_structure(struct):
    """
    for all the lines, text and structures that are not directly connected to struct,
    set alpha to 0.2 and lines to light gray
    """
    # get all the tss and barrierlesss with this struct
    highlight = []
    lines = []
    for t in tss:
        if t == struct or t.reactant == struct or t.product == struct:
            highlight.append(t)
            if not t.reactant in highlight: highlight.append(t.reactant)
            if not t.product in highlight: highlight.append(t.product)
            lines = lines + linesd[t]
        #end if
    #end for
    for b in barrierlesss:
        if b.reactant == struct or b.product == struct:
            highlight.append(b)
            if not b.reactant in highlight: highlight.append(b.reactant)
            if not b.product in highlight: highlight.append(b.product)
            lines = lines + linesd[b]
        #end if
    #end for
    for struct in linesd:
        for li in linesd[struct]:
            if not li in lines:
                li[0].set_color('lightgray')
            #end if
        #end for
    #end for
    for te in textd:
        if not textd[te] in highlight:
            te.set_alpha(0.2)
        #end if
    #end for
    for struct in imgsd:
        if not struct in highlight:
            imgsd[struct].set_alpha(0.2)
        #end if
    #end for
    plt.draw()
#end def

def save_x_values():
    """
    save the x values of the stationary points to an external file
    """
    fi = open('%s_xval.txt'%id,'w')
    fi.write(string.join(['%s %.2f'%(w.name,w.x) for w in wells],'\n')+'\n')
    fi.write(string.join(['%s %.2f'%(b.name,b.x) for b in bimolecs],'\n')+'\n')
    fi.write(string.join(['%s %.2f'%(t.name,t.x) for t in tss],'\n'))
    fi.close()
#end def

def save_im_extent():
    """
    save the x values of the stationary points to an external file
    """
    fi = open('%s_im_extent.txt'%id,'w')
    for key in imgsd:
        e = imgsd[key].get_extent()
        fi.write('%s %.2f %.2f %.2f %.2f\n'%(key.name,e[0],e[1],e[2],e[3]))
    fi.close()
#end def

def read_im_extent():
    """
    Read the extents of the images if they are present in a file_name
    """
    fname = '%s_im_extent.txt'%id
    if os.path.exists(fname):
        fi = open(fname,'r')
        a = fi.read()
        fi.close()
        a = a.split('\n')
        for entry in a:
            pieces = entry.split(' ')
            if len(pieces) == 5:
                extsd[pieces[0]] = [eval(pieces[i]) for i in range(1,5)]
            #end if
        #end for
    #end if
#end def

def main(argv):
    """
    Main method to run the PESViewer
    """
    global fname
    if len(argv) > 1: # read the arguments 
        fname = argv[1]
    #end if
    read_input() # read the input file 
    # initialize the dictionaries
    for w in wells:
        linesd[w]=[]
    #end for
    for b in bimolecs:
        linesd[b]=[]
    #end for
    for t in tss:
        linesd[t]=[]
    #end for
    for b in barrierlesss:
        linesd[b]=[]
    #end for
    read_im_extent() # read the position of the images, if known
    position() # find initial positions for all the species on the graph
    generate_lines() # generate all the line
    generate_2d_depiction() # generate 2d depiction from the smiles or 3D structure, store them in join(input_id,'_2d')
    plot() # plot the graph
#end def
    

main(sys.argv)