from anndata import AnnData
import upsetplot as upset
import pandas as pd 
from matplotlib.colors import Normalize,rgb2hex
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import superexacttestpy as stest
from math import log

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

    """
    if value == 0 : 
        return '#CACACA'
    norm = Normalize(vmin=vmin, vmax=vmax)
    cmap = get_cmap(cmap_name)  # PiYG
    rgb = cmap(norm(abs(value)))[:3]  # will return rgba, we take only first 3 so we get rgb
    color = rgb2hex(rgb)
    return color

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
    background_color:str="dark_background"
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
        Sort the results by "degree" or "p_val"
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
        Default = "dark_background" other possibility : 
            see `style.available` for list of available styles
    
    Returns
    -------
    None 

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
    if show_count : 
        show_count = '%d'
    else : 
        show_count = None

    #Construct the plot 
    if background_color != "dark_background" : 
        res = upset.UpSet(plot_data,orientation=orientation,sort_by=sort_by,show_counts=show_count,intersection_plot_elements=20,subset_size="auto")
    else : 
        res = upset.UpSet(plot_data,orientation=orientation,sort_by=sort_by,show_counts=show_count,intersection_plot_elements=20,subset_size="auto",facecolor='white')
    if color_p_val : 
        p_val = []
        for code in list(df.index) :
            tmp = list(df["p_value"].loc[[code]])[0]
            if tmp == None : 
                p_val.append(0)
            elif tmp > 0 : 
                p_val.append(-round(log(float(tmp),10),2))
            elif tmp == 0 : 
                p_val.append(-round(log(float(1*10**-320),10),2))
            else : #tmp == nan or na 
                p_val.append(0)
        vmax = max(p_val)+1/3*max(p_val)
        for i,val in enumerate (list(df.index)) : 
            col = color_map_color(p_val[i],vmax=vmax)
            pres,abs = stest.tl.decode(val,name)
            res.style_subsets(present = pres, absent=abs,facecolor=col,label=f"-log10(p_value) = {p_val[i]}")
            res.style_subsets()
    with plt.style.context(background_color) :
        res.plot(fig)