from scipy.io import loadmat
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io as sio



def _check_keys(dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict

def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    
    for strg in matobj._fieldnames:
        #print(strg, 'ffff', matobj.__dict__[strg])

        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict

def loadmat_homebrew(filename):
    """
    this function should be called instead of direct scipy.io .loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    dict = _check_keys(data)
    dictCDF = []
    for key in dict:
        if key == 'CFD':
            dictCDF = dict[key]
  
    return dict, dictCDF

def read_mat_victor(name):

    matobj, matCFD = loadmat_homebrew(name)
        

    return matCFD
    #Sch, lm y PP
    #PSI, x, y
    
    '''
    con_list = [[element for element in upperElement] for upperElement in data]
    
    vals = data #<-- set the array you want to access. 
    keys = data.dtype.descr

    print(keys)
    # zip provides us with both the x and y in a tuple.
    newData = list(zip(con_list[0], con_list[1]))
    columns = ['Sch', 'lm']# 'PP']
    df = pd.DataFrame(newData, columns=columns)
    
    Sch = data['Sch']
    lm = data['lm']
    PP = data['PP']
    
    PSI = data['PSI']
    x = data['x']
    y = data['y']
    return PSI, x, y
    '''
  


def plotFig(filename, data, fignr=1):
    d = loadmat(filename, squeeze_me=True, struct_as_record=False)
    matfig = d['hgS_070000']
    
    childs = matfig.children

    ax1 = [c for c in childs if c.type == 'axes']
    if(len(ax1) > 0):
        ax1 = ax1[0]
    legs = [c for c in childs if c.type == 'scribe.legend']
    
    if(len(legs) > 0):
        legs = legs[0]
    else:
        legs=0


    dataX = [c for c in childs if c.type == 'XData']


    pos = matfig.properties.Position
    size = np.array([pos[2]-pos[0],pos[3]-pos[1]])/96
    plt.figure(fignr,figsize=size)
    plt.clf()

    counter = 0    

    for line in ax1.children:
        if line.type == 'graph2d.lineseries':
            if hasattr(line.properties,'Marker'):
                mark = "%s" % line.properties.Marker
                if(mark != "none"):
                    mark = mark[0]
            else:
                mark = '.'
            if hasattr(line.properties,'LineStyle'):
                linestyle = "%s" % line.properties.LineStyle
            else:
                linestyle = '-'
            if hasattr(line.properties,'Color'):
                r,g,b =  line.properties.Color
            else:
                r = 0
                g = 0
                b = 1
            if hasattr(line.properties,'MarkerSize'):
                marker_size = line.properties.MarkerSize
            else:
                marker_size = -1                
            x = line.properties.XData
            y = line.properties.YData
            if(mark == "none"):
                plt.plot(x,y,linestyle=linestyle,color=[r,g,b])
            elif(marker_size==-1):
                plt.plot(x,y,marker=mark,linestyle=linestyle,color=[r,g,b])
            else:
                plt.plot(x,y,marker=mark,linestyle=linestyle,color=[r,g,b],ms=marker_size)
        elif line.type == 'text':
            if counter == 0:
                plt.xlabel("$%s$" % line.properties.String,fontsize =16)
            elif counter == 1:
                plt.ylabel("$%s$" % line.properties.String,fontsize = 16)
            elif counter == 3:
                plt.title("$%s$" % line.properties.String,fontsize = 16)
            counter += 1  
    if hasattr(ax1.properties, 'XGrid'): plt.grid(ax1.properties.XGrid)      
    

    if(hasattr(ax1.properties,'XTick')):
        if(hasattr(ax1.properties,'XTickLabelRotation')):
            plt.xticks(ax1.properties.XTick,ax1.properties.XTickLabel,rotation=ax1.properties.XTickLabelRotation)
        else:
            plt.xticks(ax1.properties.XTick,ax1.properties.XTickLabel)
    if(hasattr(ax1.properties,'YTick')):
        if(hasattr(ax1.properties,'YTickLabelRotation')):
            plt.yticks(ax1.properties.YTick,ax1.properties.YTickLabel,rotation=ax1.properties.YTickLabelRotation)
        else:
            plt.yticks(ax1.properties.YTick,ax1.properties.YTickLabel)
    plt.xlim(ax1.properties.XLim)
    plt.ylim(ax1.properties.YLim)




    if legs:        
        leg_entries = tuple(['$' + l + '$' for l in legs.properties.String])
        py_locs = ['upper center','lower center','right','left','upper right','upper left','lower right','lower left','best','best']
        MAT_locs=['North','South','East','West','NorthEast', 'NorthWest', 'SouthEast', 'SouthWest','Best','none']
        Mat2py = dict(zip(MAT_locs,py_locs))
        location = legs.properties.Location
        plt.legend(leg_entries,loc=Mat2py[location])


    plt.show()

filename_fig = "../../Desktop/Dropbox/transfer_inormation_prague/data/imput/CFD_rats/CFD_S10_Ch3_PP.fig"
filename_mat = "../../Desktop/Dropbox/transfer_inormation_prague/data/imput/CFD_rats/CFD_S10.mat"
#data = read_mat_victor(filename_mat)
#plotFig(filename_fig, data, fignr=1)

#CFD
# _S10_Ch1_Sch.fig, 
#está Sch, lm y PP, que son las tres zonas del hipocampo.
#Para cada zona del hipocampo hay 3 variables PSI, x, y.
