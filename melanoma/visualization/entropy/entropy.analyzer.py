def entropyCalculator(v,ground,sky):

    '''
    This function computes the entropy of a vector under defined range of expression
    '''

    # f.1. calculate the probability distribution
    step=(sky-ground)/len(cellIDs)
    k=numpy.arange(ground,sky+step,step)
    n,bins=numpy.histogram(v,bins=k)
    y=[]
    y=numpy.array(n)
    y=y/float(sum(y))

    # f.2. calculate entropy
    s=scipy.stats.entropy(y,base=2)

    return s


def histogrammer(theData):

    '''
    This function creates a histogram.
    '''    

    x=[]; y=[]
    
    n,bins=numpy.histogram(theData,bins=int(numpy.sqrt(len(theData))))

    halfBin=(bins[1]-bins[0])/2.
    for bin in bins:
        center=bin+halfBin
        x.append(center)
    x.pop()
    y=numpy.array(n)
    y=list(y/float(sum(y)))

    return x,y
