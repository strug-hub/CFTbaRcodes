#------------------------------------------------------------------
# Original source:
#    Natsuhiko Kumasakaa, Ritei Shibata
#    High-dimensional data visualisation: The textile plot
#    https://doi.org/10.1016/j.csda.2007.11.016
#    Computational Statistics and Data Analysis
#
# Python implementation:
#    Scott Mastromatteo
#------------------------------------------------------------------

import numpy as np
import pandas as pd
import scipy.linalg
import networkx as nx
import matplotlib.pyplot as plt

def _flatten_list(lst):
    ''' Internal function '''

    n=lst[0].shape[0]
    cols = []
    for x in lst:
        if x.ndim == 1: 
            cols.append(x.reshape(n))
        else:
            for i in range(x.shape[1]):
                cols.append(x[:,i].reshape(n))
    return np.stack(cols, axis=1)

def _textile_transform_with_missing(Xj, wj, qj, eigen_choice=1, zscale=0):
    ''' Internal function '''

    n=len(Xj[0]) ; p=len(Xj)

    W = np.stack(wj, axis=1)
    d = np.sum(W, axis=1)
    if sum(d==0) > 0:
        raise ValueError("One or more records is completely missing data.")
    
    Wpad = np.stack([wj[i] for i in qj], axis=1)
    X=_flatten_list(Xj)
    
    # A11 ------------------------
    d11 = np.ones(n).transpose().dot(W)
    A11 = -W.transpose().dot(W/d.reshape(n, 1)) + np.diag(d11)
    
    # A12 ------------------------
    A12= Wpad.transpose().dot(Wpad*X/d.reshape(n, 1))
    d12 = np.diag(Wpad.transpose().dot(X))
    d12 = np.diag(d12)
    
    to_remove = [] ; s=set()
    for j,q_ in enumerate(qj):
        if q_ in s:
            to_remove.append(j)
            d12[q_] += d12[j]
        s.add(q_)
    
    A12 = np.delete(A12, to_remove, axis=0) - np.delete(d12, to_remove, axis=0)
    
    # A22 ------------------------
    d22_full = X.transpose().dot(Wpad)*Wpad.transpose().dot(X)/np.ones(n).reshape(1,n).dot(Wpad)
    d22 = np.zeros((d22_full.shape))
    skip = 0       
    for k,j in enumerate(qj):
        j=j+skip
        if j != k: skip+=1
        while j <= k:
            d22[j,k] = d22_full[j,k] ; d22[k,j] = d22_full[k,j]
            j+=1
    
    A22 = -(Wpad*X).transpose().dot(Wpad*X/d.reshape(n, 1)) + d22
    
    # B ------------------------
    B_full = np.dot(X.transpose(), Wpad*X) - d22_full
    B = np.zeros((B_full.shape))
    
    skip = 0       
    for k,j in enumerate(qj):
        j=j+skip
        if j != k: skip+=1
        while j <= k:
            B[j,k] = B_full[j,k] ; B[k,j] = B_full[k,j]
            j+=1
    
    # A ------------------------
    A11inv = np.linalg.pinv(A11, rcond=1e-10)
    
    A=(np.dot(A12.transpose(), A11inv).dot(A12) - A22)
    
    eigenvectors, eigenvalues, _ = scipy.linalg.lapack.dsygv(A, B, uplo="U")
    # in ascending order, last col has biggest eigenvalue
    beta = eigenvectors[:,-eigen_choice]

    # arbitray constant
    if zscale == 0:
        alpha=np.dot(A11inv, A12).dot(beta)
    else:
        z=np.ones(p) * zscale 
        alpha=np.dot(A11inv, A12).dot(beta) + (np.ones(p) - np.dot(A11inv, A11)).dot(z)
    
    return(A, B, alpha, beta)

def _textile_transform_no_missing(Xj, qj, eigen_choice=1, a0=0):
    ''' Internal function '''

    p=len(Xj) ; n=len(Xj[0])
    
    X = _flatten_list(Xj)
    
    XTX = X.transpose().dot(X)
    ones = np.ones(n).reshape(n, 1)
    XT11TX = X.transpose().dot(ones).dot(ones.transpose()).dot(X)
    A = (XTX - XT11TX/n)/p
    
    B_full = XTX - XT11TX/n
    B = np.zeros((B_full.shape))
    
    skip = 0
    for k,j in enumerate(qj):
        j=j+skip
        if j != k: skip+=1
        while j <= k:
            B[j,k] = B_full[j,k] ; B[k,j] = B_full[k,j]
            j+=1
        
    eigenvectors, eigenvalues, _ = scipy.linalg.lapack.dsygv(A, B, uplo="U")    
    # in ascending order, last col has biggest eigenvalue
    beta = eigenvectors[:,-eigen_choice]
    
    alpha = np.zeros(len(Xj))
    for j,xj in enumerate(Xj):
        k=[i for i,q in enumerate(qj) if q == j]
        alpha[j]=a0-(np.ones(n).reshape(1,n).dot(xj)/n).dot(beta[k])

    return(A, B, alpha, beta)

def textile_transform(X, is_categorical=None, eigen_choice=1, category_options=dict(), get_location_scale=False):
    '''  
    Textile plot transformation

    X : [p x n numpy matrix] Use None for missing data. Note: uses faster method if no data is missing
    
    is_categorical : [list] of length p. True if x[i] is categorical variable, False if quantitative.
                    Use is_categorical=None if no variables are categorical.
                    Use is_categorical="all" if ALL variables are categorical.
    
    category_options: [dict[int]:list] Key=variable index, value=all categorical classes available to that variable.
                      If category options are not provided, the categories are inferred and placed in sorted order.
    
    eigen_choice: [int] eigenvector to use for transformation. Default=1 (largest eigenvalue). Must be <= p
    
    get_location_scale: [boolean] If true, returns the location and scale vectors as a tuple along with the transformed matrix.
                        Intended for internal use for genetic application to calculate height for each allele.
    Returns p x n numpy matrix transformed data.
    '''
    
    Xj=[]; wj=[] ; qj=[] ; j=0
    
    x = [col.reshape(col.shape[0]) for col in np.hsplit(X, X.data.shape[1])]
    if is_categorical is None:
        is_categorical = [False for x_ in x]
    
    if is_categorical == "all":
        is_categorical = [True for x_ in x]

    p=len(x) ; n=len(x[0])

    def is_none(item):
        return item is None or np.isnan(item) or str(item) == "NA"
    def most_frequent(itr):
        lst = [x for x in itr if not is_none(x)]
        return max(set(lst), key = lst.count)
    
    for xi, cat in zip(x, is_categorical):
        wi = np.array([0 if is_none(x_) else 1 for x_ in xi])
        wj.append(wi)
        if cat:
            common = most_frequent(xi)
            xi=[common if w_ == 0 else x_ for w_,x_ in zip(wi,xi)]
            
            if j in category_options:
                cats = category_options[j]
            else:
                cats = sorted(list(dict.fromkeys(xi)))
            
            cat_dict = {cat_: l for l,cat_ in enumerate(cats)}
            Z = np.zeros((len(xi), len(cats)))
            for k,x_ in enumerate(xi):
                Z[k][cat_dict[x_]] = 1
                
            C=pd.get_dummies(cats, drop_first=True).values
            X_=Z.dot(C)
            Xj.append(X_)
            for col in range(X_.shape[1]):
                qj.append(j)
        else:
            xi=[0 if w_ == 0 else x_ for w_,x_ in zip(wi,xi)]
            Xj.append(np.array(xi).reshape(len(x[0]),1))
            qj.append(j)
        j+=1
        
    if sum([sum(1-w_) for w_ in wj]) > 0:
        A, B, alpha, beta = \
            _textile_transform_with_missing(Xj, wj, qj, eigen_choice=eigen_choice)
    else:
        A, B, alpha, beta = \
            _textile_transform_no_missing(Xj, qj, eigen_choice=1)
    
    alpha = alpha * (n*p)
    beta = beta * (n*p)

    yj=[]
    for j,xj in enumerate(Xj):
        k=[i for i,q in enumerate(qj) if q == j]
        yj.append(xj.dot(beta[k]) + alpha[j])
    
    Y=np.stack(yj, axis=1)
    #y0=(W*Y).sum(axis=1)/Y.shape[1]
    #yj = [y0] + yj
    #Y=np.stack(yj, axis=1)

    #print("A:\n", A.round(4))
    #print("B:\n", B.round(4))
    #print("alpha:\n", alpha.round(4))
    #print("beta:\n", beta.round(4)) 
    #print("Y:\n", Y.round(4))
    if get_location_scale:
        return(Y, alpha, beta)
    
    return Y

def genotype_textile(vcf_records, plot=True, haplotype=False, set_missing_to_ref=False):
    '''  
    Textile plot transformation for genetic data

    vcf_records : [list] A list of p VCF records (generated from PyVCF) to transform.
    haplotype : [boolean] Whether to transform genotype (GT) or phased haplotypes from records
                Defaults to genotype.
    set_missing_to_ref: [boolean] Set all GT values with missing data to 0/0.
    plot: [boolean] Whether to plot the results.

    Returns a list of length p corresponding to y-position of each variant in textile plot.
    If plot=True, draws a networkx plot of the variants.
    '''
    
    if not haplotype:
        allele_matrix=[]
        for record in vcf_records:
            if len(record.alleles) != 2:
                raise ValueError("Only bialleleic variants are allowed when haplotype=False : " + str(record))        
            samples=[]
            for sample in record.samples:
                gt=sample.gt_alleles
                if None in gt: 
                    if set_missing_to_ref:
                        samples.append(0)
                    else:
                        samples.append(None) 
                    continue
                samples.append(int(gt[0]) + int(gt[1]))
    
            if not set_missing_to_ref and len([x for x in samples if x is not None]) < 1:
                raise ValueError("All data for variant is missing: " + str(record))        
            allele_matrix.append(samples)
        options = {j : [0,1,2] for j in range(len(allele_matrix)) }

    else:
        allele_matrix=[]
        for record in vcf_records:
            samples=[]
            for sample in record.samples:
                a,b = sample.gt_alleles
                if a is None:
                    if set_missing_to_ref: a = 0
                else: a = int(a)
                if b is None:
                    if set_missing_to_ref: b = 0
                else: b = int(b)
                samples.append(a)
                samples.append(b)
            if not set_missing_to_ref and len([x for x in samples if x is not None]) < 1:
                raise ValueError("All data for variant is missing: " + str(record))        
            allele_matrix.append(samples)
        options = dict()
        for j,v in enumerate(allele_matrix):
            opt = set(v) ; opt.discard(None)
            options[j] = sorted(list(opt))
    
    M = np.stack(allele_matrix, axis=1)
    Y, alpha, beta = textile_transform(M, is_categorical="all", category_options=options, get_location_scale=True)
    
    y=[] ; k=0
    for j in range(M.shape[1]):
        y.append( {options[j][0] : alpha[j]} )
        for cat in options[j][1:]:
            y[j][cat] = alpha[j] + beta[k]
            k+=1

    if plot: _plot_textile(allele_matrix, y)
    
    return y

def _plot_textile(allele_matrix, y):
    ''' Internal function for plotting genetic data '''

    # create a list of all the edges (allele x -> allele y) in the population
    nvar=len(allele_matrix[:-1])
    edges=[]
    for i in range(nvar-1):
        edge_list = [(allele_matrix[i][j], allele_matrix[i+1][j]) for j in range(len(allele_matrix[i])) \
                      if allele_matrix[i][j] is not None and allele_matrix[i+1][j] is not None]
        edges.append(edge_list)
    
    # variant1 : {allele1 : nodeid, allele2: nodeid, ...}, variant2: ...
    node_dict=dict() ; nodeid=0
    graph_edges=[]
    graph_weights=[]
    graph_pos=dict()
    
    # build data structures for first variant
    node_dict[0] = dict()
    alleles = { e[0] for e in edges[0] }
    for i,allele in enumerate(alleles):
        node_dict[0][allele] = nodeid
        graph_pos[nodeid] = (0, y[0][allele])
        nodeid+=1

    # build data structures for the rest of the variants
    for varid,edge in enumerate(edges):
        to_alleles= { e[1] for e in edge }
        node_dict[varid+1] = dict() 
        for i,allele in enumerate(to_alleles):
            node_dict[varid+1][allele] = nodeid
            graph_pos[nodeid] = (varid+1, y[varid+1][allele])
            nodeid+=1
        
        edgecount = dict()
        for e in edge:
            if e not in edgecount: edgecount[e] = 0
            edgecount[e] += 1

        graph_edge = [(node_dict[varid][gt[0]], node_dict[varid+1][gt[1]]) \
                       for gt in edgecount]
        graph_weight = [edgecount[gt] for gt in edgecount]
        graph_edges.extend(graph_edge)
        graph_weights.extend(graph_weight)

    # build graph
    G = nx.DiGraph()
    for w,e in zip(graph_weights, graph_edges):
        G.add_edge(e[0], e[1], weight=w)
    
    labeldict = dict()
    for var in node_dict:
        for allele in node_dict[var]:
            labeldict[node_dict[var][allele]] = allele
    
    edges = G.edges()
    weights = np.log10([G[u][v]["weight"] for u,v in edges])+0.1
        
    nx.draw_networkx(G, graph_pos, node_size=300, font_size=10, 
                     width=weights, arrowstyle="-", font_color="white",
                     labels=labeldict, with_labels = True)
        
    plt.xlim([-1, len(allele_matrix)])
    
    # Set margins for the axes so that nodes aren't clipped
    ax = plt.gca()
    ax.margins(0.20)
    plt.show()