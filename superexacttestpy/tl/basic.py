from anndata import AnnData
from math import log,exp
from copy import deepcopy
from random import choices
from collections import Counter
import pandas as pd 
from re import finditer
from itertools import product


def basic_tool(adata: AnnData) -> int:
    """Run a tool on the AnnData object."""
    print("Implement a tool to run on the AnnData object.")
    return 0

def len_data(data): 
    """ Function to get the len of all thesub list """
    data_l = []
    for values in data : 
        data_l.append(len(values))
    return data_l

def log_choose(n,k): 
    """ Approximation of log binomial coefficient nCk """
    m = min(k,n-k)
    res = 0.0
    if m == 0 : 
        return res 

    for i in range(1,m+1):
        res += log((n-i+1))-log(i)

    return res 


def dhyper(x,w,b,n,logP): 
    """
    Probability of getting x white balls out of n draws from an urn with w white balls and b black balls
    """
    if x > w or x > n or b+x < n : 
        res = 0
        if logP : 
            return log(res)

    res = log_choose(w,x)+log_choose(b,n-x)
    res = res-log_choose(w+b,n)

    if not logP : 
        return exp(res)

    return res

def dmv_hyper(x,nl,set_len,p,n,logp): 
    """Function to calculate the probability density of x elems in all of the subsets"""
    a_size = max(set_len) - x + 1
    min_l = min(set_len)
    f1 = [0 for i in range(a_size)]
    i0 = 0

    if(nl == 2):
        p=dhyper(x,set_len[0],n-set_len[0],set_len[1],logp)
        return p
    for i in range(a_size) : 
        f1.append(0)
    p = 0
    #from inner-most to outer-most
    for i in range(1,nl+1): 
        if(i == 1) : 
            l = x
            f1[0] = dhyper(x,l,n - l,set_len[nl - 1],i0)
            for l in range(x+1,min(min_l,n+x-set_len[nl-1])): 
                f1[l - x] = f1[l - x -1] * ((n - l+1-set_len[nl - 1] + x)/(l - x))  * (l/(n -l+1))
        f0 = deepcopy(f1)
        if(nl - i >= 2):
            for k in range(x,min_l+1) : 
                f1[k - x] = 0
                l = max(x,k+set_len[nl - i] - n)
                temp = dhyper(l,set_len[nl - i],n - set_len[nl - i],k,i0)
                f1[k - x] += temp * f0[l - x]
                for l in range(l+1,k+1): 
                    temp = temp * ((set_len[nl - i]-l+1)/l) * ((k-l+1) /(n - set_len[nl - i]-k+l))
                    f1[k - x] += temp * f0[l - x]
            
        #final integration
        j = max(x,set_len[1]+set_len[0] - n)
        temp = dhyper(j,set_len[1],n - set_len[1],set_len[0],i0)
        p += temp * f1[j - x]
        for j in range(j+1,min_l+1): 
            temp = temp * ((set_len[1]-j+1)/j) * ((set_len[0]-j+1) /(n - set_len[1]-set_len[0]+j))
            p += temp * f1[j - x]

    if (p > 1 and not logp) : 
        p = 1.0
    if ( p <= 0 ): 
        p = 2.2*10**308
    if(logp) and p > 0 :
        p = log(p)
    if (logp and p<= 0) :  
        p = 0
    return p 

def dp_sets(x,data,n,logp=False): 
    """ Function to calculate the probability density of the intersection
    x (int): size of gene intersection. 
    data (list): the list of the data to compute 
    n (int): total number of gene 
    logP (bool): is the p value to return is log or not """
    l_data = len_data(data)
    mini = min(l_data)
    nl = len(data) # len of the dataset 
    for nb in l_data : 
        if nb > n or nb < x : 
            if not logp : 
                return 0 
            print("Invalid input")
            return False 
    if mini < 0 :
        if not logp : 
            return 0 
        print("Invalid input")
        return False 
    
    if nl < 2 :
        print("data should have at least 2 entries")
        return False 
    
    res = dmv_hyper(x=x,nl=len(data),l=len_data(data),n=n,p=0.0,logp=logp)
    return res

def intersect(data):
    """Performs set intersection on multiple input vector"""
    res = [list(set(data[0]).intersection(data[1]))]
    for i in range(1,len(data)): # 1 bc we compute the first value in the initialisation
        res.append(list(set(res[-1]).intersection(data[i])))
    return res[-1]

def log_choose_logval(n,k,logval):
    """log binomial coefficient nCk""" 
    m = min(k,n-k)
    result = 0.0
    if m == 0 : 
        return result
    for i in range(m) : 
        result = result+logval[n-i-1]-logval[i]
    return result

def dhyper_logval(x,w,b,n,logval,logp): 
    """probability of getting x white balls out of n draws from an urn with w white balls and b black balls"""
    if x > w or x > n or b+n < n : 
        result = 0 
        if logp : 
            result = log(result) ## Strange bc log(0) is invalid 
    else : 
        result = log_choose_logval(w,x,logval)+log_choose_logval(b,n-x,logval)
        result -= log_choose_logval(w+b,n,logval)
        if not logp : 
            result = exp(result)
    return result

def dmv_hyper_logval(x:int,nl:int,set_len:list,n:int,p:float,logval:list,logp:bool=False): 
    """Function to get the distribution of multiset intersection test"""
    a_size = max (set_len)-x+1
    min_l = min(set_len)
    f1 = [0 for i in range(a_size)]

    if nl == 2 : 
        p = dhyper_logval(x,set_len[0],n-set_len[0],set_len[1],logval,logp) 
        return p 
    
    p = 0
    for i in range(1,nl) : 
        if i == 1  : 
            l = x 
            f1[0] = dhyper_logval(x,l,n-l,set_len[nl-1],logval,logp)
            for ll in range(x+1,min(min_l,n+x-set_len[nl-1])+1) : 
                f1[ll-x] = f1[ll - x - 1] * ((n-ll+1-set_len[nl-1]+x)/(ll-x)) * (l / (n-ll+1))
        f0 = deepcopy(f1)
        if nl-i >= 2 : 
            for k in range(x,min_l+1) : 
                f1[k-x] = 0
                l = max(x,k+set_len[nl-i]-n)
                temp = dhyper_logval(l,set_len[nl-i],n-set_len[nl-i],k,logval,logp)
                f1[k-x] += temp*f0[l-x]
                for ll in range(l+1,k+1): 
                    temp = temp*((set_len[nl-i]-ll+1)/ll)*((k-ll+1)/(n-set_len[nl-i]-k+ll))
                    f1[k-x] += temp*f0[ll-x]
        j = max(x,set_len[1]+set_len[0]-n)
        temp = dhyper_logval(j,set_len[1],n-set_len[1],set_len[0],logval,logp)
        p += temp*f1[j-x]
        for jj in range(j+1,min_l+1): 
            temp = temp*((set_len[1]-jj+1)/jj)*((set_len[0]-jj+1)/(n-set_len[1]-set_len[0]+jj))
            p += temp*f1[jj-x]

        if p>1 : 
            p = 1.0
        if p < 0 : 
            p = 2.2*10**308
        if logp : 
            p = log(p)
        return p 

def pmvhyper(x:int,nl:int,set_len:list,n,p,lower_tail:bool=False,logp:bool=False): 
    """Function to get the probability of multiset intersection"""
    tiny = 1.0*10**308
    i0=0
    p0=0.0
    min_l = min(set_len)
    pp = [0 for i in range(min_l+1)]
    logVal = [log(i) for i in range(1,n+1)]
    
    if x == 0 : 
        p = dmv_hyper_logval(x,nl,set_len,n,0,logVal,logp)
        if not lower_tail : 
            p = 1.0 - p 
        if p > 1 : 
            p = 1.0
        if p < 0 : 
            p = 2.2*10**308
        if logp : 
            p = log(p)
        return p
    Xmean = n 
    for i in range(nl): 
        Xmean = Xmean *set_len[i] / n 

    for i in range(min_l+1): 
        pp[i] = 0
    p = 0 
    if x > Xmean : 
        i = x+1 
        while i <= min_l : 
            p0=dmv_hyper_logval(i,nl,set_len,n,0.0,logVal,logp)
            pp[i] = p0 
            if i > x+1 : 
                break 
            if p0 <= tiny : 
                break 
            if i > (x+1) and p0/pp[i-1] <= 0.01 : 
                break 
            if i > x+1 : 
                break 
 
        p = pp[x+1]#Not in the C code but it's works 
    if not lower_tail : 
        p = 1-p
    if p > 1 : 
        p == 1.0
    if p < 0 : 
        p = 2.2*10**308
    if logp:
        p = log(p)
    return p 

def cpsets_sim(x:int,set_len:list,n:int,lower_tail:bool=True,logp:bool=False,number_sim:int=10000): 
    """Function to simulate the intersection probability on multisets"""
    nl = len (set_len)
    count = 0
    for i in range(1,number_sim): 
        for length in set_len:
            tmp = choices(range(1,n),k=length)
            tmp = Counter(tmp).most_common() # We count the occurence of eatch value 
            for tuple in tmp : 
                if tuple[0] ==  nl : 
                    if tuple[1] <= x : 
                        count += 1
                    break
                elif tuple[0] > nl : 
                    break # break the loop bc if the nl value is 
    p = count/number_sim
    if not lower_tail : 
        return 1-p 
    if logp : 
        return log(p)
    return p 

def cpsets(x:int,data:list,n:int,lower_tail:bool=True,logp:bool=False,p_val_sim:bool=False,number_sim:int=100000):
    """Function to compute the distribution function of multiset intersection test. It's possible to simulate a p-value"""
    nl = len(data)
    l_data = len_data(data)
    if nl < 2 : 
        return False
    if x < 0 or n < 1 : 
        return False 
    if p_val_sim : 
        return cpsets_sim(x,l_data,n,lower_tail,logp,p_val_sim,number_sim)
    res = pmvhyper(x,nl,l_data,n,0,lower_tail,logp)
    return res 

def prod(data,total): 
    """Function to calculate the multiplication of each item in the list divided by the total """ 
    res = 1 # bc if it's 0 0*anything = 0 
    for i in range(len(data)):
        res *= len(data[i])/total
    return total*res 

def mset (x:int,n:int,lower_tail:bool=True,logp:bool=False) :
    """A wrapper to call cpsets
    x (int): list of sets 
    n (int): background size
    lower_tail (bool): if TRUE (default), probability is P[overlap < m], otherwise, P[overlap >= m], where m is the number of elements shared by all sets.
    logp (bool): if TRUE, probabilities p are given as log(p).
    """
    l_data = len_data(x)
    L = [i/n for i in l_data]
    if n < 1 :
        print("invalid input")
        return False
    for nb in l_data : 
        if nb>n or nb==0 : 
            print("invalid input")
            return False 

    intersects = intersect(x)
    Obs = len(intersects)
    Exp = prod(x,n)
    if Obs == 0 : 
        if lower_tail: 
            p = 0
        else : 
            p = 1
     
    else : 
        p = cpsets(Obs-1,x,n,lower_tail,logp)
    
    return {"Intersection":intersects, "FE" : round(Obs/Exp,1), "p-value":p}

def mk_barcode_degree(n,degree): 
    res = []

    if degree == None : 
        degree = [i for i in range(1,n+1)]

    if type(degree) == list : 
        for comb in product(["0","1"],repeat=n): 
            tmp = list(comb)
            if tmp.count('1') in degree:
                res.append("".join(tmp)) 
        return res

    if type(degree) == int: 
        if degree > n : 
            print("Impossible")
            return False
        for comb in product(["0","1"],repeat= n) : 
            tmp = list(comb)
            if tmp.count('1')==degree:
                res.append("".join(tmp)) 
        return res

    else : 
        print("Degree should be a list of int of an int")
        return False 

def inc_intersect(x,degree:int) : 
    if type(x) != list : 
        print("Input data must be list")
        return False 
    nl = len(x)
    if nl < 2 : 
        print("Input data should have at least two entries")
        return False 

    barcodes = mk_barcode_degree(nl,degree)
    otab = [0 for i in range(len(barcodes))]
    d = {}
    for i in range(len(barcodes)) : 
        i1 = [match.start() for match in finditer('1',barcodes[i])]
        if len(i1) == 1 : 
            otab[i] = len(x[i1[0]])
        else : 
            tmp = []
            for idx in i1 : 
                tmp.append(x[idx])
            otab[i] = len(intersect(tmp))
        d[barcodes[i]] = [otab[i]]
    return pd.DataFrame(d)

def intersect_elements(x): 
    if type(x) != list : 
        print("Input data must be list")
        return False 
    nl = len(x)
    if(nl < 2) : 
        print('Input x should have at least two entries')
    
    x_all = []
    for i in x : 
        x_all += i
    allE = list(set(x_all)) # get just the unique value

    barcode = []
    for elem in allE : 
        tmp = []
        for i in range(nl) : 
            if elem in x[i] : 
                tmp.append("1")
            else : 
                tmp.append("0")
        barcode.append("".join(tmp))
    return pd.DataFrame({"Entry":allE, "barcode":barcode})

def exclusive_intersect0(x): 
    intersects = intersect_elements(x)
    nl = len(x)
    barcode = mk_barcode_degree(nl,None)
    otab = {}
    for code in barcode : 
        nb = intersects[intersects.barcode == code].shape[0]
        otab[code] = [nb]
    return pd.DataFrame(otab)

def exc2inc_intersect(df:pd.DataFrame): 
    barcode = list(df.columns)
    d = {}
    for code in barcode : 
        i1 = [match.start() for match in finditer('1',code)]
        nb = 0

        for code1 in barcode : 
            i2 = [match.start() for match in finditer('1',code1)]
            if all(i in i2 for i in i1) : 
                nb += int(df[code1])
        
        d[code] = [nb]
    return pd.DataFrame(d)

def enumerate_intersect_sizes(x,degree:int=-1) : 
    if degree > 0 : 
        return inc_intersect(x,degree)
    otab = exclusive_intersect0(x)
    return exc2inc_intersect(otab)

def formating(barcode,names,collapse=' & '): 
    res = []
    barcode = list(barcode)
    for i in range(len(barcode)): 
        if barcode[i] == "1" : 
            res.append(names[i])
    return collapse.join(res)

def supertest(data,n:int,names:list=[],degree:int=-1,lower_tail=False): 
    if type(data) != list : 
        print("Input must be a list")
        return False
    if len(names) != len(data) :
        for i in range(len(data)): 
            names.append(f"Set{i+1}")
   
    x = data 
    size = [len(list(set(x[i]))) for i in range(len(x))]
    df_overlap_size = enumerate_intersect_sizes(x,degree=degree)
    
    barcode = list(df_overlap_size.columns)
    overlap_size = [int(df_overlap_size[code]) for code in barcode]

    if n > 0 : 
        for siz in size : 
            if siz > n : 
                print("Background population size should not be smaller than set size")
                return False

        overlap_expected = [None for i in df_overlap_size]
        for i in range(df_overlap_size.size): 
            nb = [match.start() for match in finditer('1',barcode[i])]
            if len(nb) > 1 : 
                overlap_expected[i] = n
                for size_nb in nb : 
                    overlap_expected[i] *= size[size_nb]/n
            else : 
                overlap_expected[i] = None

        p_val = [0 for i in range(df_overlap_size.size)]
        for idx,code in enumerate(barcode): 
            i1 = [match.start() for match in finditer('1',code)]
            if len(i1) == 1 : 
                p_val[idx] = None
            elif overlap_size[idx] == 0 :
                p_val[idx] = 1
            else : 
                L = [x[i] for i in i1]
                # print(cpsets(list(df_overlap_size.loc[0])[idx]-1,L,n))
                p_val[idx] = cpsets(list(df_overlap_size.loc[0])[idx]-1,L,n,lower_tail=lower_tail)

    nl = len(x)
    if type(degree) == int : 
        if degree <= 0 : 
            degree = [i for i in range(nl)]

    elif type(degree) == list : 
        for d in degree : 
            if d < 1 or d > nl : 
                print("Invalid degree value")
                return False
    else : 
        print("Degree should be a list of int or an int")
        return False 
    odegree = [i.count('1') for i in barcode]
    otab = exc2inc_intersect(exclusive_intersect0(x))
    otab = [int(otab[code]) for code in barcode]
    if len(otab) == 0 : 
        print("No data for output")
        return None
    etab = overlap_expected

    el = intersect_elements(x)
    elements = []
    for i in barcode : 
        tmp = []
        ncode = len(barcode[0])
        regex = ""
        for range_code in range(nl) :
            if i[range_code] == "1" : 
                regex += "1"
            else: 
                regex += "[01]"
        for elem in el[el["barcode"].str.contains(regex)]["Entry"] :# Search if the barcode match with the regex and get the entry name. 
            tmp.append(elem)
        tmp = sorted(tmp)
        elements.append(", ".join(tmp))
        
    decode = []
    for code in barcode : 
        decode.append(formating(code,names))

    FE = [otab[i]/etab[i] if etab[i]!= None and otab[i]!=None else None for i in range(len(etab))]

    res = pd.DataFrame({"Intersection":decode,"degree":odegree,"Observed_overlap":otab,"Expected_overlap":etab,"FE":FE,"p-value":p_val,"Elements":elements})
    res.index = barcode
    return res

def decode(barcode,name): 
    res_pres = []
    res_abs= []
    for i in range(len(list(barcode))) :
        if barcode[i] == "1" : 
            res_pres.append(name[i])
        else : 
            res_abs.append(name[i])
    return res_pres, res_abs