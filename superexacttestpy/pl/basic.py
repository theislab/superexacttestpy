import upsetplot as upset
import pandas as pd 
from matplotlib.colors import Normalize,rgb2hex
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import superexacttestpy as stest
from math import log
import matplotlib

def color_map_color(value, cmap_name='YlOrRd', vmin=0, vmax=1):
    """
    Create a color map from a given value

    Parameters
    ----------
    value : float
        Value to map
    cmap_name : str
        Name of the color palette 
    vmin : float
        Minimum value of the color map
    vmax : float
        Maximum value of the color map
    
    Returns
    -------
    str
        Hexadecimal color code
    
    Example
    -------
    >>> color_map_color(0.5, cmap_name='YlOrRd', vmin=0, vmax=1)
    ... '#fd8c3c'

def basic_plot(adata: AnnData) -> int:
    """Generate a basic plot for an AnnData object."""
    print("Import matplotlib and implement a plotting function here.")
    return 0


def color_map_color(value, cmap_name='YlOrRd', vmin=0, vmax=1):
    if value == 0 : 
        return '#CACACA'
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = get_cmap(cmap_name)  # PiYG
    rgb = cmap(norm(abs(value)))[:3]  # will return rgba, we take only first 3 so we get rgb
    color = rgb2hex(rgb)
    return color

import numpy as np
def get_specific_color_gradient(colormap,inputList,**kwargs):
    vmin = kwargs.get('vmin','blaq')
    vmax = kwargs.get('vmax','blaq')
    cm = plt.get_cmap(colormap)
    if vmin=='blaq' or vmax=='blaq':
        if type(inputList)==list:
            cNorm = matplotlib.colors.Normalize(vmin=min(inputList), vmax=max(inputList))
        else:
            cNorm = matplotlib.colors.Normalize(vmin=inputList.min(), vmax=inputList.max())
    else:
        cNorm = matplotlib.colors.Normalize(vmin=vmin, vmax = vmax)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    scalarMap.set_array(inputList)
    colorList=scalarMap.to_rgba(inputList)
    return scalarMap,colorList

def plot(
    data:list,
    n:int,
    name:list,
    degree=-1,
    sort_by="degree",
    show_count=True,
    orientation:str="horizontal",
    color_p_val:bool=True,
    size:tuple=(10,5),
    background_color:str="seaborn-notebook"
    ) -> pd.DataFrame :
    """
    Plot the results of the superexact test

    Parameters
    ----------
    data : list
        List of data to compute with superexact test
    n : int
        Background size
    name : list
        List of names of each set 
    degree : int
        Intersection degree
    sort_by : str
        Sort the results by "degree" or "cardinality"
    show_count : bool
        Show the number of genes in each set
    orientation : str
        Orientation of the plot. "horizontal" or "vertical" 
    color_p_val : bool
        Color the bars by their p-value
    size : tuple
        Size of the plot
    background_color : str
        Background color of the plot. 
        Default = "seaborn-notebook" other possibility : 
            see `style.available` for list of available styles
    
    Returns
    -------
    df : pd.DataFrame
        Dataframe with the results of the superexact test 

    Example
    -------
    >>> plot(data, n, name, degree=-1, sort_by="degree", show_count=True, orientation="horizontal", color_p_val=True, size=(10,5), background_color="dark_background")
    ... None
    
    """ 
    if type(degree)==list : 
        df = stest.tl.supertest(data,n,name,degree=-1,lower_tail=True)
        df = df[df['degree'].isin(degree)]
    elif type(degree)==int : 
        df = stest.tl.supertest(data,n,name,degree,lower_tail=True)
    else : 
        print("degree should be a list or an int")
        return False

    res_intersect = []
    res_overlap = []

    for elem in list(df["Intersection"]):
        if " & " in elem:
            elem = elem.split()
            elem = list(filter(lambda a: a != "&", elem))  # remouve the &
            res_intersect += [elem]

        else:
            res_intersect += [[elem]]

    for nb in list(df["Observed_overlap"]):
        res_overlap.append(nb)
    plot_data = upset.from_memberships(res_intersect, res_overlap)
    for elem in list(df["intersection"]) :
        if " & " in elem: 
            elem = elem.split()
            elem = list(filter(lambda a: a != "&", elem)) # remouve the & 
            res_intersect += [elem]
            
        else : 
            res_intersect += [[elem]]
    
    for nb in list(df["observed_overlap"]) : 
        res_overlap.append(nb)
    plot_data = upset.from_memberships(res_intersect,res_overlap)
    fig = plt.figure(figsize=size)
    ax = plt.subplot()
    if show_count:
        show_count = '%d'
    else:
        show_count = None

    # Construct the plot
    if background_color != "dark_background":
        res = upset.UpSet(plot_data, orientation=orientation, sort_by=sort_by, show_counts=show_count,
                          intersection_plot_elements=20, subset_size="auto")
    else:
        res = upset.UpSet(plot_data, orientation=orientation, sort_by=sort_by, show_counts=show_count,
                          intersection_plot_elements=20, subset_size="auto", facecolor='white')

    if color_p_val:
        p_val = []
        for code in list(df.index):
            tmp = list(df["p-value"].loc[[code]])[0]
            if tmp == None:
    else : 
        show_count = None

    #Construct the plot 
    res = upset.UpSet(
        plot_data,
        element_size=None,
        orientation=orientation,
        sort_by=sort_by,
        show_counts=show_count,
        intersection_plot_elements=20,
        subset_size="auto"
        )
        
    if color_p_val : 
        p_val = []
        for code in list(df.index) :
            tmp = list(df["p_value"].loc[[code]])[0]
            if tmp == None : 
                p_val.append(0)
            elif tmp > 0:
                p_val.append(-round(log(float(tmp), 10), 2))
            elif tmp == 0:
                p_val.append(-round(log(float(1 * 10 ** -320), 10), 2))
            else:  # tmp == nan or na
                p_val.append(0)
        vmax = max(p_val) + 1 / 3 * max(p_val) if vmax is None else vmax
        for i, val in enumerate(list(df.index)):
            col = color_map_color(p_val[i], vmax=vmax_input if vmax is None else vmax)
            pres, abs = stest.tl.decode(val, name)
            # res.style_subsets(present = pres, absent=abs, facecolor=col) # label=f"-log10(p_value) = {p_val[i]}")

        # fig.show()
    # cbar states
    tickscolorbar = np.arange(0, vmax, (vmax - 0) / 4).astype(int)
    scalarmap, colorList = get_specific_color_gradient('YlOrRd', np.array(tickscolorbar), vmin=0, vmax=vmax)
    colorList = scalarmap.to_rgba(p_val)
    for i, val in enumerate(list(df.index)):
        col = color_map_color(p_val[i], vmax=vmax)
        pres, abs = stest.tl.decode(val, name)
        res.style_subsets(present=pres, absent=abs, facecolor=colorList[i] if p_val[
                                                                                  i] != 0 else '#CACACA')  # label=f"-log10(p_value) = {p_val[i]}")

    with plt.style.context(background_color):
        res.plot(fig)

    ax = plt.subplot(2, 12, 1)
    ax.legend().set_visible(False)
    plt.axis('off')
    cbar = fig.colorbar(scalarmap, orientation="vertical", format="%.0f", ticks=tickscolorbar, ax=ax, anchor=(1.0, 1.0))
    cbar.set_label(r'$-log_{10}(P_{val})$', fontsize=12)
        vmax = max(p_val)+1/3*max(p_val)
        for i,val in enumerate (list(df.index)) : 
            col = color_map_color(p_val[i],vmax=vmax)
            pres,abs = stest.tl.decode(val,name)
            res.style_subsets(present = pres, absent=abs,facecolor=col,label=f"-log10(p_value) = {p_val[i]}")
            res.style_subsets()
    with plt.style.context(background_color) :
        res.plot(fig = fig)
    plt.suptitle("test")
    plt.show()
    return df
