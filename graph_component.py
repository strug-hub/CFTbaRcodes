from dash import dcc, dash_table
from plotly import graph_objs as go
import vcf
import textile_plot as textile
import numpy as np
import pandas as pd
import cvxpy as cp
import random
from pairing import pair
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

EDGE_COLOR="rgba(0,0,0,0.5)"
EDGE_COLOR_HIGHLIGHT="rgba(255,0,0,0.95)"
EDGE_SCALE_FACTOR=10

MISSING_TO_REF=True
VARIANT_INDEL="INDEL"
VARIANT_SNP="SNP"
VARIANT_ETC="ETC"

SYMBOL=["♈","♉","♊","♋","♌","♍","♎","♏","♐","♑","♒","♓","⛎"]
EMPTY_SYMBOL="◻️"

class Node:
    def __init__(self, variant, gt, seq):
        self.variant = variant
        self.gt = gt
        self.bin = id
        self.seq = seq
        self.x = None
        self.y = None
        self.dy = None
        self.shape = "circle" ; self.color = "blue"
        if variant.get_class() == VARIANT_INDEL:
            self.shape = "diamond" ; self.color = "orange"

    def coords(self):
        return (self.x, self.y)
    def dcoords(self):
        return (self.x, self.dy)
    def id(self):
        return( pair(self.variant.id, self.gt) )

class Variant:
    def __init__(self, id, record):
        self.record = record
        self.id = id
        self.alleles=[]
        for gt,allele in enumerate(record.alleles):
            new_allele = Node(self, gt, allele)
            self.alleles.append(new_allele)
        
            self.gt=[]
            for sample in record.samples:
                self.gt.extend(sample.gt_alleles)
                none_value = 0 if MISSING_TO_REF else None
                self.gt = [ none_value if a is None else int(a) for a in self.gt ]

    def set_coordinates(self, y, dy):
        for gt in y:
            self.alleles[gt].y = y[gt]
            self.alleles[gt].dy = dy[gt]
            self.alleles[gt].x = self.id
    def get_sample(self, i):
        return self.alleles[self.gt[i]]
    def get_class(self):
        if self.record.is_snp: return VARIANT_SNP
        if self.record.is_indel: return VARIANT_INDEL
        return VARIANT_ETC
    
class Edge:
    def __init__(self, node1, node2, nsamp):
        self.node1 = node1
        self.node2 = node2
        self.frequency = 1
        self.total = nsamp
        self.haps = []

    def observed(self):
        self.frequency += 1

    def weight(self):
        return self.frequency/self.total
    
    def get_x(self, pad=False):
        x = [self.node1.coords()[0], self.node2.coords()[0]]
        if pad: x.append(None)
        return x
    
    def get_y(self, pad=False):
        y= [self.node1.coords()[1], self.node2.coords()[1]]
        if pad: y.append(None)
        return y
    def get_dy(self, pad=False):
        dy= [self.node1.dcoords()[1], self.node2.dcoords()[1]]
        if pad: dy.append(None)
        return dy

    def id(self):
        return( pair(self.node1.id(), self.node2.id()) )

class Haplotype:
    def __init__(self, id):
        self.id = id
        self.hap = []
    
    def add(self, edge):
        self.hap.append(edge)
    def to_vector(self):
        return [edge.id() for edge in self.hap]
    def nodes(self):
        return [self.hap[0].node1] + [edge.node2 for edge in self.hap]


def define_haps(vcf_records=None, missing_to_ref=True):
    
    if vcf_records is None or len(vcf_records) < 1:        
        VCF="cftr.vcf.gz"
        vcf_reader = vcf.Reader(filename=VCF)
        vcf_records = [record for record in vcf_reader]
        selected=[5,17]
        
    # Create variant objects
    variants = [] ; nodes = []
    for i,record in enumerate(vcf_records):
        variant = Variant(i, record)
        variants.append(variant)
        for node in variant.alleles:
            nodes.append(node)
            
    # Calculate textile heights
    height_dict = textile.genotype_textile(vcf_records, plot=False, haplotype=True, set_missing_to_ref=MISSING_TO_REF)
    #height_dict = optimize_spread(height_dict, 15)
    height_dict_spread = optimize_spread(height_dict)

    for i,variant in enumerate(variants):
        variant.set_coordinates( height_dict[i], height_dict_spread[i] )        
    
    nodes = [node for node in nodes if node.y is not None]       
    
    # Create haplotype objects
    nsamp = len(variants[0].gt) ; nvar = len(variants)
    haplotypes=[] ; edge_dict = dict() ; edges = []
    for i in range(nsamp):
        haplotype = Haplotype(i)

        for j in range(nvar-1):
            allele_from = variants[j].get_sample(i)
            allele_to = variants[j+1].get_sample(i)
            edge_key = (allele_from, allele_to)
            if edge_key in edge_dict:
                edge = edge_dict[edge_key]
                edge.observed()
            else:
                edge = Edge(allele_from, allele_to, nsamp)
                edge_dict[edge_key] = edge
                edges.append(edge)
        
            haplotype.add(edge)
        haplotypes.append(haplotype)
    
    
    groups = dict()
    for i,haplotype in enumerate(haplotypes):
        hap_nodes = haplotype.nodes()
        group = tuple(hap_nodes[s].gt for s in selected)
        if not group in groups: groups[group] = []
        groups[group].append(i)
    
    def fake_rsid():
        return "rs"+ str(round(random.random()*1e6))
    rsids= [ fake_rsid() for s in selected ]
        
    total = sum([len(groups[g]) for g in groups])
    rows = []
    for group in groups:
        size = len(groups[group])
        r = {"freq" : str(round(size*100/total, 2)) + "%" }
        for i,gt in enumerate(group): r[SYMBOL[i]] = gt
        rows.append(r)
    
    
    tooltip = { SYMBOL[i]: {"value": rsid, "use_with": "both"} for i,rsid in enumerate(rsids) }
    
    return dash_table.DataTable(rows, id="haplotype-table", 
                                row_selectable="single",
                                tooltip=tooltip,
                                tooltip_delay=0)






# optimization problem to spread out nodes
def optimize_spread(y, height_fraction=10):    
    def spread(points, min_distance, max_left, max_right):
        points = np.array(sorted(points))
        #if len(y) > 4: print(y)
        p = cp.Variable(points.shape[0])
          
        # constraints
        dist_constr = [p[i] - p[i-1] >= min_distance for i in range(1, points.shape[0])]
        left_constr = [p[0] >= max_left]
        right_constr = [p[-1] <= max_right]
        constrs = dist_constr + left_constr + right_constr
        
        # setup problem + solve
        obj = cp.Minimize(sum(cp.abs(p - points)))
        problem = cp.Problem(obj, constrs)
        problem.solve()
        
        return p.value
    
    y_all = [var[allele] for var in y for allele in var]
    y_min, y_max = min(y_all), max(y_all)
    gap = (abs(y_max) + abs(y_min))/height_fraction
    
    y_spread = []
    for var in y:
        alleles = [] ; y_original = []
        for allele in var:
            alleles.append(allele)
            y_original.append(var[allele])
            
        y_order = np.argsort(y_original)
        y_rank = np.argsort(y_order)
        y_new = spread(y_original, gap, y_min, y_max)
        y_spread.append({ allele: y_new[i] for i,allele in zip(y_rank,alleles) })
    
    return y_spread

def plot_graph(vcf_records=None):
    
    if vcf_records is None or len(vcf_records) < 2:        
        VCF="cftr.vcf.gz"
        vcf_reader = vcf.Reader(filename=VCF)
        vcf_records = [record for record in vcf_reader]
    
    # Create variant objects
    variants = [] ; nodes = []
    for i,record in enumerate(vcf_records):
        variant = Variant(i, record)
        variants.append(variant)
        for node in variant.alleles:
            nodes.append(node)
    
    # Calculate textile heights
    height_dict = textile.genotype_textile(vcf_records, plot=False, haplotype=True, set_missing_to_ref=MISSING_TO_REF)
    #height_dict = optimize_spread(height_dict, 15)
    height_dict_spread = optimize_spread(height_dict)

    for i,variant in enumerate(variants):
        variant.set_coordinates( height_dict[i], height_dict_spread[i] )        
    
    nodes = [node for node in nodes if node.y is not None]       
    
    # Create haplotype objects
    nsamp = len(variants[0].gt) ; nvar = len(variants)
    haplotypes=[] ; edge_dict = dict() ; edges = []
    for i in range(nsamp):
        haplotype = Haplotype(i)

        for j in range(nvar-1):
            allele_from = variants[j].get_sample(i)
            allele_to = variants[j+1].get_sample(i)
            edge_key = (allele_from, allele_to)
            if edge_key in edge_dict:
                edge = edge_dict[edge_key]
                edge.observed()
            else:
                edge = Edge(allele_from, allele_to, nsamp)
                edge_dict[edge_key] = edge
                edges.append(edge)
        
            haplotype.add(edge)
        haplotypes.append(haplotype)

    
    #collapse_haps(haplotypes)
    
    node_traces = get_node_traces(nodes)
    edge_traces = get_edge_traces(edges)
    haplotype_traces = get_haplotype_traces(haplotypes)
    
    def fake_rsid():
        return "rs"+ str(round(random.random()*1e6))

    y = [node.y for node in nodes]
    y_min, y_max = min(y), max(y)
    yspan = y_max - y_min
    yann = y_min - yspan*0.3
    
    label_traces =[
        go.Scatter(
            x=[x for x,variant in enumerate(variants)],
            y=[y_min - yspan*0.15 for x,variant in enumerate(variants)],
            text=[EMPTY_SYMBOL for x,variant in enumerate(variants)],
            mode="text",
            hoverinfo="name",
            showlegend=True,
            customdata=[ "label", [], EMPTY_SYMBOL] 
            )
        ]
    
    fig = go.Figure(data=label_traces + edge_traces + node_traces,
             layout=go.Layout(
                plot_bgcolor="#F6F6F6",
                hovermode="x",
                hoverdistance = 1,
                margin=dict(b=20,l=5,r=5,t=40),
                transition={"duration": 300, "easing": "elastic-in"},
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                )
             )

    for x,variant in enumerate(variants):
        fig.add_annotation(x=x, y=yann,
                text=fake_rsid(), #variant.ID,
                showarrow=False,
                align="left",
                textangle=90
                )

    return dcc.Graph(figure=fig, id="graph-plot")

def flatten(t):
    return [item for sublist in t for item in sublist]

def get_node_traces(nodes):
    
    node_bins = dict();
    for node in nodes:
        if node.bin not in node_bins:
            node_bins[node.bin] = []
        node_bins[node.bin].append(node)

    node_traces = []
    for i in node_bins:
        
        node_traces.append(
                go.Scatter(
                    x=[node.x for node in node_bins[i]],
                    y=[node.y for node in node_bins[i]],
                    text=[node.seq for node in node_bins[i]],
                    mode="markers+text",
                    hoverinfo="text",
                    textfont_color="white",
                    showlegend=True,
                    marker=dict(
                        color=[node.color for node in node_bins[i]],
                        symbol=[node.shape for node in node_bins[i]],
                        size=20, line=dict(color="#888888", width=2)),
                    customdata=[ "node", 
                                [node.y for node in node_bins[i]],
                                [node.dy for node in node_bins[i]]]
                    )
                )
    return node_traces

def get_edge_traces(edges):
    
    weights = [edge.weight() for edge in edges]
    bin_values = [b/100 for b in range(0, 100, 10)][1:]
    bins = np.digitize(weights, bin_values, right=True)
    bin_values.append(1)

    edge_bins = dict();
    for b,edge in zip(bins,edges):
        if b not in edge_bins:
            edge_bins[b] = []
        edge_bins[b].append(edge)

    edge_traces = []
    for i in edge_bins:

        edge_traces.append(
                go.Scatter(
                    x=flatten([edge.get_x(pad=True) for edge in edge_bins[i]]),
                    y=flatten([edge.get_y(pad=True) for edge in edge_bins[i]]),
                    mode="lines",
                    hoverinfo="none",
                    line_color=EDGE_COLOR,
                    line=dict(width=bin_values[i]*EDGE_SCALE_FACTOR),
                    showlegend=False,
                    customdata=[ "edge", 
                                flatten([edge.get_y(pad=True) for edge in edge_bins[i]]),
                                flatten([edge.get_dy(pad=True) for edge in edge_bins[i]])]
                    )
                )
    return edge_traces

def get_haplotype_traces(haplotypes):
    
    haplotype_traces = []

    for i,hap in enumerate(haplotypes):
        haplotype_traces.append(
                go.Scatter(
                    x=flatten([edge.get_x(pad=True) for edge in hap.hap]),
                    y=flatten([edge.get_y(pad=True) for edge in hap.hap]),
                    mode="lines",
                    hoverinfo="none",
                    line_color=EDGE_COLOR_HIGHLIGHT,
                    line=dict(width=EDGE_SCALE_FACTOR*1.2),
                    visible=False,
                    showlegend=False,
                    customdata=[ "haplotype", 
                                flatten([edge.get_y(pad=True) for edge in hap.hap]),
                                flatten([edge.get_dy(pad=True) for edge in hap.hap]),
                                i]
                    )
                )
    return haplotype_traces
    
CLIENTSIDE_UPDATE= \
    """
    function(hover_data, click_data, selected_ids, figure) {
        SYMBOL="""+str(SYMBOL)+""";

        console.log(figure);
        var triggered = dash_clientside.callback_context.triggered;
        if (triggered.length > 0 ){
        
            // update edge color with selected_ids haplotype ------------------------
            if(triggered[0]["prop_id"]=="haplotype-table.selected_row_ids"){
    
                var idx = selected_ids[0];
    
                for (let i = 0; i < figure["data"].length; i++) {
                    if ("customdata" in figure["data"][i] && figure["data"][i]["customdata"][0] == "haplotype"){
                        figure["data"][i]["visible"] = false;
    
                        if (idx == figure["data"][i]["customdata"][3]){
                            figure["data"][i]["visible"] = true;
                        }
                    }
                }
                
                return { "data": figure["data"], "layout": figure["layout"] }
            }
            
            
            // handle click ------------------------
            if(triggered[0]["prop_id"]=="graph-plot.clickData"){
                console.log(click_data["points"]);

                idx = click_data["points"][0]["x"];
                for (let i = 0; i < figure["data"].length; i++) {
                    data=figure["data"][i];
                    if ("customdata" in data && data["customdata"][0] == "label"){
                        
                        if (idx in data["customdata"][1]){
                            console.log('remove it');
                        }
                        else{
                            data["customdata"][1].push(idx);
                        }
                        data["customdata"][1].sort(function(a, b){return a - b});

                        new_text=[]; k=0;
                        for (let j = 0; j < data["text"].length; j++) {
                            new_text[j] = data["customdata"][2];
                            if (j in data["customdata"][1]){
                                new_text[j] = SYMBOL[k];
                                k=k+1;
                            }
                        }
                        data["text"] = new_text;
                        return { "data": figure["data"], "layout": figure["layout"] }

                    }
                }

            
            
            
            }
            
            // hover on plot ------------------------
            if(triggered[0]["prop_id"]=="graph-plot.hoverData"){
                
                if (typeof hover_data == "undefined"){
                    return figure;
                } else{
                    var idx = hover_data["points"][0]["x"]
                }
        
                if ("hover_id" in figure["data"][0]){
                    if (figure["data"][0]["hover_id"] == idx){
                        return figure;
                    }
                }
                figure["data"][0]["hover_id"] = idx;    
                for (let i = 0; i < figure["data"].length; i++) {
                    if ("customdata" in figure["data"][i]){
                        var trace = figure["data"][i]["customdata"][0];
                        
                        if (trace == "node" || trace == "edge" || trace == "haplotype"){

                            var x = figure["data"][i]["x"];
                            var init_y = figure["data"][i]["customdata"][1];
                            var spread_y = figure["data"][i]["customdata"][2];
                        
                            var update_y = [];
                            for (let j = 0; j < x.length; j++) {
                                update_y[j] = init_y[j]
                                if (x[j] == idx){
                                    update_y[j] = spread_y[j];
                                }
                            }
                            figure["data"][i]["y"] = update_y;
                        }
                    }
                }
                
                return { "data": figure["data"], "layout": figure["layout"] }
            }
        }

        return figure
    }
    """









'''
def collapse_haps(vcf_records=None):
    
    if vcf_records is None or len(vcf_records) < 2:        
        VCF="cftr.vcf.gz"
        vcf_reader = vcf.Reader(filename=VCF)
        vcf_records = [record for record in vcf_reader]
    
    # Create variant objects
    variants = [] ; nodes = []
    for i,record in enumerate(vcf_records):
        variant = Variant(i, record)
        variants.append(variant)
        for node in variant.alleles:
            nodes.append(node)
            
    # Calculate textile heights
    height_dict = textile.genotype_textile(vcf_records, plot=False, haplotype=True, set_missing_to_ref=MISSING_TO_REF)
    #height_dict = optimize_spread(height_dict, 15)
    height_dict_spread = optimize_spread(height_dict)

    for i,variant in enumerate(variants):
        variant.set_coordinates( height_dict[i], height_dict_spread[i] )        
    
    nodes = [node for node in nodes if node.y is not None]       
    
    # Create haplotype objects
    nsamp = len(variants[0].gt) ; nvar = len(variants)
    haplotypes=[] ; edge_dict = dict() ; edges = []
    for i in range(nsamp):
        haplotype = Haplotype(i)

        for j in range(nvar-1):
            allele_from = variants[j].get_sample(i)
            allele_to = variants[j+1].get_sample(i)
            edge_key = (allele_from, allele_to)
            if edge_key in edge_dict:
                edge = edge_dict[edge_key]
                edge.observed()
            else:
                edge = Edge(allele_from, allele_to, nsamp)
                edge_dict[edge_key] = edge
                edges.append(edge)
        
            haplotype.add(edge)
        haplotypes.append(haplotype)
    
    import plotly.figure_factory as ff
    np.random.seed(1)
    
    #haps = [haplotype.to_vector() for haplotype in haplotypes]
    #X = np.matrix(haps)
    
    X = np.random.rand(15, 12)
    fig = ff.create_dendrogram(X)
    
    return dcc.Graph(figure=fig, id="graph-dendogram")
    #plot_bgcolor="#F6F6F6",
    
    
    haps = [haplotype.to_vector() for haplotype in haplotypes]
    X = np.matrix(haps)
    dendrogram = sch.dendrogram(sch.linkage(X, method="ward"), p=15, truncate_mode="lastp")
    plt.xlabel("Haplotypes") ; plt.ylabel("Euclidean distances")
    plt.show()
    
    from sklearn.cluster import AgglomerativeClustering

    cluster = AgglomerativeClustering(n_clusters=2, affinity='euclidean', linkage='ward')

    def recusive_cluster(X, ids=None):

        if X.shape[0] < 3: return None

        if ids is None:
            ids = [i for i in range(X.shape[0])]
        
        pred = cluster.fit_predict(X)
        
        idx0 = [i for i,p in enumerate(pred) if p==0]
        ids0 = [ids[i] for i in idx0]
        X0 = np.take(X, idx0, 0)
        idx1 = [i for i,p in enumerate(pred) if p==1]
        ids1 = [ids[i] for i in idx1]
        X1 = np.take(X, idx1, 0)

        res0 = recusive_cluster(X0, ids0)
        res1 = recusive_cluster(X1, ids1)

        if res0 is None or res
        
        return res0 + res1

    recusive_cluster(X)

    
    
    
    
    
    cluster = AgglomerativeClustering(n_clusters=X.shape[0], affinity='euclidean', linkage='ward')
    fit = cluster.fit(X)
    
    ward = ward_tree(X, return_distance=True)
'''
