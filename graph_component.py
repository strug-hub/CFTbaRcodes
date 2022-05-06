from dash import dcc, dash_table
from plotly import graph_objs as go
import vcf
import textile_plot as textile
import numpy as np
import cvxpy as cp
import random
from pairing import pair
import json
from dash.dependencies import Input, Output, State
from cryptography.fernet import Fernet

KEY = Fernet.generate_key()
def encrypt(message):
    return Fernet(KEY).encrypt(message)
def decrypt(token):
    return Fernet(KEY).decrypt(token)

def get_haps(figure):
    for data in figure["data"]:
        if data["customdata"][0] == "haplotype":
            hap_token=data["customdata"][1]
            string = decrypt(hap_token.encode()).decode();
            haps = string.split(",")
            haps = [h.split(" ") for h in haps]
            return haps

EDGE_COLOR="rgba(0,0,0,0.5)"
EDGE_COLOR_HIGHLIGHT="rgba(255,0,0,0.7)"
EDGE_SCALE_FACTOR=10

MISSING_TO_REF=True
VARIANT_INDEL="INDEL"
VARIANT_SNP="SNP"
VARIANT_ETC="ETC"

SYMBOL=["♈","♉","♊","♋","♌","♍","♎","♏","♐","♑","♒","♓","⛎"]
EMPTY_SYMBOL="◻️"

def json_Node(j):
    serialized = json.loads(j)
    node = Node(None, serialized["gt"], serialized["seq"])
    node.x = serialized["x"] ; node.y = serialized["y"]
    node.shape = serialized["shape"] ; node.color = serialized["color"]
    return node

class Node:
    def __init__(self, variant, gt, seq):
        self.variant = variant
        self.gt = gt
        self.seq = seq
        self.x = None
        self.y = None
        self.shape = "circle" ; self.color = "blue"
        if variant.get_class() == VARIANT_INDEL:
            self.shape = "diamond" ; self.color = "orange"

    def coords(self):
        return (self.x, self.y)
    def id(self):
        return( pair(self.variant.id, self.gt) )
    def to_json(self):
        serialized = dict(gt=self.gt, seq=self.seq,
                         x=self.x, y=self.y,
                         shape=self.shape, color=self.color)
        return json.dumps(serialized)
    
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

    def set_coordinates(self, y):
        for gt in y:
            self.alleles[gt].y = y[gt]
            self.alleles[gt].x = self.id
    def get_sample(self, i):
        return self.alleles[self.gt[i]]
    def get_class(self):
        if self.record.is_snp: return VARIANT_SNP
        if self.record.is_indel: return VARIANT_INDEL
        return VARIANT_ETC

    def to_json(self):
        serialized = dict(gt=self.gt, seq=self.seq,
                         x=self.x, y=self.y,
                         shape=self.shape, color=self.color)
        return json.dumps(serialized)
    
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
    
    def gts(self):
        return ((self.node1.gt, self.node2.gt))

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
    
    def stringify(self):
        return " ".join([str(node.gt) for node in self.nodes()])
    
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
    height_dict = optimize_spread(height_dict, 12)

    for i,variant in enumerate(variants):
        variant.set_coordinates( height_dict[i] )        
    
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
        
    string_haps = ",".join([hap.stringify() for hap in haplotypes])
    
    haplotype_trace =[
        go.Scatter(
            x=[], y=[], showlegend=False,
            customdata=[ "haplotype", encrypt(string_haps.encode()).decode() ]
            )
        ]

    
    node_traces = get_node_traces(nodes)
    edge_traces = get_edge_traces(edges)
    
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
            showlegend=False,
            customdata=[ "label" ] 
            )
        ]
    
    fig = go.Figure(data=haplotype_trace + label_traces + edge_traces + node_traces,
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
        if node.gt not in node_bins:
            node_bins[node.gt] = []
        node_bins[node.gt].append(node)

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
                    showlegend=False,
                    marker=dict(
                        color=[node.color for node in node_bins[i]],
                        symbol=[node.shape for node in node_bins[i]],
                        size=20, line=dict(color="#888888", width=2)),
                    customdata=[ "node", [node.gt for node in node_bins[i]] ] 
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
                    customdata=[ "edge", [edge.gts() for edge in edge_bins[i]] ]
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
                    customdata=[ "haplotype", i ]
                    )
                )
    return haplotype_traces

def unflatten(lst):
    v = [] ; tmp = []
    for item in lst:
        if item is None: 
            v.append(tuple(tmp)) ; tmp =[]
        else: 
            tmp.append(item)
    if len(tmp) > 0: v.append(tuple(tmp))
    return v
    
def draw_haplotype(figure, haplogroup):
    
    # remove existing traces
    to_rm = []
    for i,data in enumerate(figure["data"]):
        if "customdata" in data and data["customdata"][0] == "haplogroup":
            to_rm.append(i)
    for i in sorted(to_rm, reverse=True):
        del figure["data"][i]
    
    increment = 1/len(haplogroup)
    edges = [(i, i+1) for i,_ in enumerate(haplogroup[0][:-1])]
    edge_values = [dict() for _ in haplogroup[0][:-1]]
    for haplotype in haplogroup:
        for i,hi in enumerate(haplotype[:-1]):
            value = (hi, haplotype[i+1])
            if value not in edge_values[i]: 
                edge_values[i][value] = 0
                
            edge_values[i][value] += increment

    node_values = [dict() for _ in haplogroup[0]]
    for haplotype in haplogroup:
        for i,v in enumerate(haplotype):
            if v not in node_values[i]: 
                node_values[i][v] = 0
                
            node_values[i][v] += increment

    # place new traces under node and over edges
    placement = -1
    for i,(d1,d2) in enumerate(zip(figure["data"][:-1], figure["data"][1:])):
        if "customdata" in d1 and "customdata" in d2:
            if d1["customdata"][0] == "edge" and d2["customdata"][0] == "node":
                placement = i
                break
        
    traces_to_add=[]
    # add edges
    for data in figure["data"]:
        if "customdata" in data and data["customdata"][0] == "edge":

            edge_xs = unflatten(data["x"])
            edge_ys = unflatten(data["y"])
            
            for i,edge_x in enumerate(edge_xs):
                match_idx = edges.index(edge_x)
                
                value = tuple(str(v) for v in data["customdata"][1][i])
                
                if value in edge_values[match_idx]:
                    proportion = edge_values[match_idx][value]
                    trace = {"customdata": ["haplogroup"], 
                             "hoverinfo": "none", 
                             "line": {"color": EDGE_COLOR_HIGHLIGHT,
                                      "width": proportion*EDGE_SCALE_FACTOR*1.2},
                             "mode": "lines", 
                             "showlegend": False, 
                             "x": [x for x in edge_x], 
                             "y": [y for y in edge_ys[i]],
                             "type": "scatter"}
                    
                    traces_to_add.append(trace)
                    
    # add nodes
    for data in figure["data"]:
        if "customdata" in data and data["customdata"][0] == "node":

            for i,x in enumerate(data["x"]):
                value = str(data["customdata"][1][i])
                
                if value in node_values[x]:
                    proportion = node_values[x][value]
                    
                    trace = {"customdata": ["haplogroup"], 
                             "hoverinfo": "text", 
                             "marker": {"color": EDGE_COLOR_HIGHLIGHT,
                                        "symbol": "circle",
                                        "size": 22 + 10*proportion},
                             "text": str(round(proportion*100,1)) + "%",
                             "showlegend": False, 
                             "x": [x], 
                             "y": [data["y"][i]],
                             "type": "scatter"}
                    
                    traces_to_add.append(trace)

    for trace in traces_to_add:
        figure["data"].insert(placement, trace)


    return figure

def update_symbols(figure, symbols):
    
    for data in figure["data"]:
        if data["customdata"][0] == "label":
            data["text"] = symbols
    return figure


"""
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
    height_dict = optimize_spread(height_dict, 15)

    for i,variant in enumerate(variants):
        variant.set_coordinates( height_dict[i] )        
    
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

"""