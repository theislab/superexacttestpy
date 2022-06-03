from anndata import AnnData
from math import log,exp


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

def logChoose(n,k): 
    """ Approximation of log binomial coefficient nCk """
    m = min(k,n-k)
    res=0.0
    if m == 0 : 
        return res 

    for i in range(1,m+1):
        res+=log((n-i+1))-log(i)

    return res 


def dhyper(x,w,b,n,logP): 
    """
    Probability of getting x white balls out of n draws from an urn with w white balls and b black balls
    """
    if x>w or x>n or b+x < n : 
        res=0
        if logP : 
            return log(res)

    res = logChoose(w,x)+logChoose(b,n-x)
    res=res-logChoose(w+b,n)

    if not logP : 
        return exp(res)

    return res

def dmvHyper2(x,nL,L,p,n,logp): 
    """Function to calculate the probability density of x elems in all of the subsets"""
    aSize=max(L) - x + 1
    minL=min(L)
    f1=[]
    f0=[]
    i0=0

    if(nL == 2):
        p=dhyper(x,L[0],n-L[0],L[1],logp)
        return p
    for i in range(aSize) : 
        f1.append(0)
        f0.append(1) #To not multiply by 0 at the 69 lines 
    p=0
    #from inner-most to outer-most
    for i in range(1,nL): 
        if(i==1) : 
            l = x
            f1[0]=dhyper(x,l,n - l,L[nL - 1],i0)
            for l in range(x+1,min(minL,n+x-L[nL-1])): 
                f1[l - x] = f1[l - x -1] * ((n - l+1-L[nL - 1] + x)/(l - x))  * (l/(n -l+1))
        if(nL - i>=2):
            for k in range(x,minL) : 
                f1[k - x]=0
                l=max(x,k+L[nL - i] - n)
                temp = dhyper(l,L[nL - i],n - L[nL - i],k,i0)
                f1[k - x] += temp * f0[l - x]
                for l in range(l+1,k): 
                    temp = temp * ((L[nL - i]-l+1)/l) * ((k-l+1) /(n - L[nL - i]-k+l))
                    f1[k - x] += temp * f0[l - x]
                    #print(k-x)
            
        #final integration
        j=max(x,L[1]+L[0] - n)
        temp=dhyper(j,L[1],n - L[1],L[0],i0)
        p += temp * f1[j - x]
        for j in range(j+1,minL): 
            temp=temp * ((L[1]-j+1)/j) * ((L[0]-j+1) /(n - L[1]-L[0]+j))
            p += temp * f1[j - x]

    if (p > 1 and not logp) : 
        p = 1.0
    if ( p < 0 ): 
        p = 2.2*10**308
    if(logp) and p > 0 :
        p = log(p)
    if (logp and p<= 0) :  
        p = 0
    print("update")
    return p 

def DpSets(x,data,n,logP=False): 
    """ Function to calculate the probability density of the intersection
    x (int): size of gene intersection. 
    data (list): the list of the data to compute 
    n (int): total number of gene 
    logP (bool): is the p value to return is log or not """
    l_data = len_data(data)
    mini = min(l_data) 
    nL = len(data) # len of the dataset 
    if mini < 0 or any(l_data)>n or any (l_data)<x  : 
        if not logP : 
            return 0 
        print("invalid input")
        return False 
    
    if nL < 2 :
        print("data should have at least 2 entries")
        return  False 
    
    res = dmvHyper2(x=x,nL=len(data),L=len_data(data),n=n,p=0.0,logp=logP)
    return res