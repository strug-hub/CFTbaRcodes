from dash import Dash, html, dash_table
import pandas as pd
import plotly.graph_objects as go
import urllib, json
from dash import Dash, html, dcc, dash_table
import random
from dash.dependencies import Input, Output, State
import dash
import vcf
import matplotlib.pyplot as plt
import networkx as nx

app = Dash(__name__)

url = 'https://raw.githubusercontent.com/plotly/plotly.js/master/test/image/mocks/sankey_energy.json'
response = urllib.request.urlopen(url)
data = json.loads(response.read())

fig_ = go.Figure(data=[go.Sankey(
    valueformat = ".0f",
    valuesuffix = "TWh",
    node = dict(
      pad = 15,
      thickness = 15,
      line = dict(color = "black", width = 0.5),
      label =  data['data'][0]['node']['label'],
      color =  data['data'][0]['node']['color']
    ),
    link = dict(
      source =  data['data'][0]['link']['source'],
      target =  data['data'][0]['link']['target'],
      value =  data['data'][0]['link']['value'],
      label =  data['data'][0]['link']['label']
  ))])

fig_.update_layout(
    hovermode = 'x',
    font=dict(size = 10, color = 'white'),
    plot_bgcolor='black',
    paper_bgcolor='black'
)

SANKEY_ = dcc.Graph(figure=fig_, id="sankey-plot_")



# --------------------------------------
# READ VARIANT DATA
# --------------------------------------
VCF="cftr.vcf.gz"
vcf_reader = vcf.Reader(filename=VCF)
vcf_records = [record for record in vcf_reader]

def fake_rsid():
    return "rs"+ str(round(random.random()*1e6))

VARIANT_INS="INS"
VARIANT_DEL="DEL"
VARIANT_SNP="SNP"
VARIANT_ETC="ETC"
def get_variant_type(record):
    if record.is_snp: return VARIANT_SNP
    elif record.is_deletion: return VARIANT_DEL
    elif record.is_indel: return VARIANT_INS
    return VARIANT_ETC    

variant_dict = {
    "id": [i for i in range(len(vcf_records))],
    "CHROM": [v.CHROM for v in vcf_records],
    "POS": [v.POS for v in vcf_records],
    "REF": [v.REF for v in vcf_records],
    "ALT": [v.POS for v in vcf_records],

    "rsid": [fake_rsid() for v in vcf_records],
    "type": [get_variant_type(v) for v in vcf_records],
    "af": [random.random()/2 for v in vcf_records],
    "called": [v.num_called for v in vcf_records],
    "nsamp": [len(v.samples) for v in vcf_records],
    }
variant_data = pd.DataFrame(data=variant_dict)

def get_hap(record, sample, hap):
    gt = record.samples[sample].gt_alleles[hap]
    if gt is None: return 0
    return int(gt)

nvar = len(vcf_records) ; nsamp = len(vcf_records[0].samples)
haps1 = [ [get_hap(v,i,0) for v in vcf_records] for i in range(nsamp) ]
haps2 = [ [get_hap(v,i,1) for v in vcf_records] for i in range(nsamp) ]
haps=haps1 + haps2

edges=[]
for i in range(nvar-1):
    edge_list = [(hap[i], hap[i+1]) for hap in haps]
    edges.append(edge_list)



# variant1 : {allele1 : nodeid, allele2: nodeid, ...}, variant2: ...
node_dict=dict() ; nodeid=0
graph_edges=[]
graph_weights=[]
graph_pos=dict()

alleles = { e[0] for e in edges[0] }
node_dict[0] = dict() 
for i,allele in enumerate(alleles):
    node_dict[0][allele] = nodeid
    graph_pos[nodeid] = (0,i)
    nodeid+=1

for varid,edge in enumerate(edges):
    
    #from_alleles = { e[0] for e in edge }
    to_alleles= { e[1] for e in edge }
    
    node_dict[varid+1] = dict() 
    for i,allele in enumerate(to_alleles):
        node_dict[varid+1][allele] = nodeid
        graph_pos[nodeid] = (varid+1,i)
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


G = nx.DiGraph()
for w,e in zip(graph_weights, graph_edges):
    G.add_edge(e[0], e[1], weight=w)

labels = nx.get_edge_attributes(G,'weight')

plt.figure(figsize=(18,6)) 
nx.draw_networkx(G, graph_pos, node_size=20, font_size=8)

# Set margins for the axes so that nodes aren't clipped
ax = plt.gca()
ax.margins(0.20)
plt.show()




'''
nvar=10 ; nhaps=100
afs = [random.random() for i in range(nvar)]
def generate_hap(n, af):
    return [int(random.random()>af[i]) for i in range(n)]
haps = [generate_hap(nvar, af=afs) for i in range(nhaps)]
'''

x = [] ; y = [] ; label = []
source = [] ; target = [] ; value = []
link_colors = [] ; lblue="rgba(135,206,235,0.6)"
for i in range(nvar):
    x.extend([i/10,i/10])
    label.extend(["X0", "X1"])
    if i < (nvar-1):      
        source.extend([(i*2),(i*2),(i*2)+1,(i*2)+1])
        target.extend([(i*2)+2,(i*2)+3,(i*2)+2,(i*2)+3])
        link_colors.extend([lblue,lblue,lblue,lblue])
        x00 = [1 if (hap[i] == 0 and hap[i+1] == 0) else 0 for hap in haps]
        x01 = [1 if (hap[i] == 0 and hap[i+1] == 1) else 0 for hap in haps]
        x10 = [1 if (hap[i] == 1 and hap[i+1] == 0) else 0 for hap in haps]
        x11 = [1 if (hap[i] == 1 and hap[i+1] == 1) else 0 for hap in haps]
    
        value.extend([sum(x00),sum(x01),sum(x10),sum(x11)])
    

fig = go.Figure(data=[go.Sankey(
    arrangement = "snap",
    valueformat = ".0f",
    valuesuffix = " haps",
    node = dict(
        label=label, x=x, pad=10,
        color="rgb(240,145,145)", 
        line = dict(color="black", width=1),
        thickness=40
        ),
    link = dict(
        source=source, target=target,
        color=link_colors,
        value=value))])

fig.update_layout(
    hovermode = 'x',
    title="CFTR2 Haplotypes",
    font=dict(size = 10, color = 'white'),
    plot_bgcolor='black',
    paper_bgcolor='black'
)

SANKEY = dcc.Graph(figure=fig, id="sankey-plot")

@app.callback(
    Output("sankey-plot", "figure"),
    Input("btn", 'n_clicks'),
    State("sankey-plot", "figure"))
def update_component_selection(clicks, figure):
    if clicks == 0 : return dash.no_update
    
    lblue="rgba(135,206,235,0.6)"
    red="rgba(255,0,0,0.8)"

    hap = haps[clicks % nhaps]
    
    figure["data"][0]["link"]["color"] = [lblue for _ in figure["data"][0]["link"]["color"]]
    i=0; j=0
    while True:
        if i > len(hap)-2: break
        i00 = j ; i01 = j+1
        i10 = j+2 ; i11 = j+3
        gt = hap[i] ; gt_ = hap[i+1]
        if gt == 0:
            if gt_ == 0: figure["data"][0]["link"]["color"][i00] = red 
            if gt_ == 1: figure["data"][0]["link"]["color"][i01] = red
        elif gt == 1:
            if gt_ == 0: figure["data"][0]["link"]["color"][i10] = red 
            if gt_ == 1: figure["data"][0]["link"]["color"][i11] = red 
        i=i+1 ; j=j+4
    
    return figure

BUTTON = html.Div([
    html.Button("Button", id="btn", n_clicks=0),
    ])

app.layout = html.Div( [SANKEY_, BUTTON, SANKEY, ] )

if __name__ == '__main__':
    app.run_server(debug=True)

