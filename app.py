import dash
from dash import Dash, html, dcc, dash_table
import dash_bio as dashbio
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd
import dash_bootstrap_components as dbc
import vcf
import random
import plotly.graph_objects as go

import dash_tabulator

# 3rd party js to export as xlsx
external_scripts = ['https://oss.sheetjs.com/sheetjs/xlsx.full.min.js']
# bootstrap css
external_stylesheets = ['https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css']

CFTR_START=117479025 ; CFTR_END=117669665 ; FLANK=10000
CFTR_VIEW="chr7:" + str(CFTR_START-FLANK) + "-" + str(CFTR_END+FLANK)

VCF="quick_cftr_merge_chr7_117287120-117715971.vcf"
vcf_reader = vcf.Reader(filename=VCF)
vcf_records = [record for record in vcf_reader]

LOGO = "assets/logo.svg"

# --------------------------------------
# READ ANNOTATION DATA
# --------------------------------------

annotation_data = pd.read_csv("cftr_exons.txt", sep="\t")

INTRON_SCALE=100
atype=[] ; x0=[] ; x1=[]
xpos = FLANK/INTRON_SCALE ; lastpos=None
introns = [] ; prevtype=None
for i, r in annotation_data.iterrows():
    if "exon" in r["id"]: a = "exon" 
    if "UTR" in r["id"]: a = "UTR" 
    atype.append(a)
    span = r["end"]-r["start"]
    
    if lastpos is None:
        x0.append(xpos)
    else:
        intronspan= r["start"] - lastpos
        if "exon" in prevtype and "exon" in r["id"] and intronspan > 0:
            introns.append({ "start": lastpos+1, 
                             "end":  lastpos+intronspan-1,
                             "type": "intron",
                             "id": "intron" + str(len(introns)+1),
                             "x0": xpos+1, 
                             "x1": xpos+(intronspan/INTRON_SCALE) })
        xpos = xpos+intronspan/INTRON_SCALE
        x0.append(xpos)
        

    xpos = xpos+span
    x1.append(xpos)
    lastpos = r["end"]
    prevtype = r["id"]

annotation_data["type"] = atype
annotation_data["x0"] = x0
annotation_data["x1"] = x1

intron_data = pd.DataFrame(introns)
annotation_data = pd.concat([annotation_data, intron_data])

# --------------------------------------
# READ VARIANT DATA
# --------------------------------------

def relative_positon(pos):
    prev_end = 0
    for i, r in annotation_data.iterrows():
        if pos >= r["start"] and pos <= r["end"]: # in exon
            relativepos = pos - r["start"]
            return r["x0"]+relativepos
        elif pos < r["start"]: # in intron
            return round(prev_end + (r["x0"] - prev_end)/2)
        prev_end=r["x1"]
    return round(prev_end + (FLANK/INTRON_SCALE)/2)

variant_dict = {
    "id": [i for i in range(len(vcf_records))],
    "CHROM": [v.CHROM for v in vcf_records],
    "POS": [v.POS for v in vcf_records],
    "REF": [v.REF for v in vcf_records],
    "ALT": [v.POS for v in vcf_records],
    "TYPE": [round(random.random()) for v  in vcf_records],
    "x": [relative_positon(v.POS) for v in vcf_records],
    "y": [1 for v in vcf_records]
    }
variant_data = pd.DataFrame(data=variant_dict)

# --------------------------------------
# CREATE APP
# --------------------------------------

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

@app.callback(
    Output('cftr-plot-container', 'children'),
    Input('flip-cftr-plot-view', 'n_clicks')
)
def make_cftr_plot(n_clicks):
    
    LINE_COL="slateblue"
    EXON_COL="lightsteelblue"
    UTR_COL="sandybrown"
    UTR_LINE_COL="firebrick"
    TEXT_COL="dimgray"
    
    if n_clicks % 2 == 1:
        EXCLUDE_EXON_LABEL=["exon6", "exon17", "exon19", "exon26"]
        df = annotation_data[annotation_data["type"] == "exon"]
        fig = px.scatter(variant_data, x="POS", y="y")

        fig.add_shape(
                type="line",
                x0=min(df["start"]), y0=-0.5, x1=max(df["end"]), y1=-0.5,
                line=dict(color=LINE_COL, width=4)
                )
        for i, r in df.iterrows():
            fig.add_shape(type="rect",
               xref="x", yref="y",
               x0=r["start"], y0=-1, x1=r["end"], y1=0,
               line=dict( color=LINE_COL, width=1 ),
               fillcolor=EXON_COL 
               )
            if r["id"] not in EXCLUDE_EXON_LABEL:
                fig.add_annotation(
                        x=(r["start"]+r["end"])/2, y=-1.5,
                        text=r["id"],
                        showarrow=False,
                        textangle=-90,
                        font=dict(
                            family="Courier New, monospace",
                            size=12,
                            color=TEXT_COL
                            )
                        )
        
        fig.update_xaxes(range=[CFTR_START-FLANK, CFTR_END+FLANK])
        fig.update_yaxes(range=[-2, 4])
    
    else:
        ROTATE_LABEL={"exon14", "3’ UTR"}
        EXCLUDE_EXON_LABEL=["exon16"]

        df = annotation_data
        vdf = variant_data ; vdf["y"] = vdf["y"]*random.random()*4
        fig = px.scatter(vdf, x="x", y="y")

        for x in set(variant_data["x"]):
            fig.add_shape(
                type="line",
                x0=x, y0=-0.5, x1=x, y1=3,
                line=dict(color=TEXT_COL, width=1)
                )
 
        fig.add_shape(
                type="line",
                x0=min(df["x0"]), y0=-0.5, x1=max(df["x1"]), y1=-0.5,
                line=dict(color=LINE_COL, width=4)
                )
        
        for i, r in df.iterrows():

            col = UTR_COL if r["type"] == "UTR" else EXON_COL
            linecol = UTR_LINE_COL if r["type"] == "UTR" else LINE_COL

            fig.add_shape(type="rect",
               xref="x", yref="y",
               x0=r["x0"], y0=-1,
               x1=r["x1"], y1=0,
               line=dict( color=linecol, width=1 ),
               fillcolor=col
               )
            if r["id"] not in EXCLUDE_EXON_LABEL:
                fig.add_annotation(
                        x=(r["x0"]+r["x1"])/2, y=-0.5,
                        text=r["id"],
                        showarrow=False,
                        textangle=0 if r["id"] in ROTATE_LABEL else -90,
                        font=dict(
                            family="Courier New, monospace",
                            size=12,
                            color=TEXT_COL
                            )
                        )
    
        fig.update_xaxes(range=[0, max(df["x1"])+FLANK/INTRON_SCALE])
        fig.update_yaxes(range=[-1.5, 4])
        fig.update_xaxes(visible=False)

    fig.update_layout(plot_bgcolor="#F6F6F6", paper_bgcolor="LightSteelBlue",
            margin={"l": 0, "r": 0, "t": 0, "b": 0})
    fig.update_xaxes(fixedrange=True, automargin=False)
    fig.update_yaxes( fixedrange=True, visible=False, automargin=True)

    
    plot = dcc.Graph(figure=fig, id="cftr-plot")
    return plot

HEADER = html.Div(children=[
    html.Img(src=LOGO, id="logo-image"),
    html.Div(children="A web application for understanding CFTR haplotypes.")
    ])

PLOT_BUTTON = dbc.Button('Switch View', id='flip-cftr-plot-view', n_clicks=0)
PLOT = html.Div([
    dcc.Loading(id='cftr-plot-container')
    ])


CURVENUMBER_DICT=dict()
FIGURE_OBJ=None
@app.callback(
    Output("cftr-mini-plot", "children"),
    Input("flip-cftr-plot-view", "n_clicks")
)
def make_cftr_miniplot(n_clicks):
    
    LINE_COL="black"
    EXON_COL="lightskyblue"
    INTRON_COL="coral"
    UTR_COL="plum"
    TEXT_COL="dimgray"
    
    ROTATE_LABEL={"exon14", "3’ UTR"}
    EXCLUDE_EXON_LABEL=["exon16"]

    df = annotation_data
    
    fig = go.Figure()

    for i, r in df.iterrows():
        if r["type"] in ["exon", "UTR"]:
            fig.add_trace(
                go.Scatter(
                    x=[r["x0"], r["x0"], r["x1"], r["x1"], r["x0"]],
                    y=[0, 2, 2, 0, 0],
                    fill='toself', fillcolor=UTR_COL if r["type"] == "UTR" else EXON_COL,
                    hoveron = 'fills', # select where hover is active
                    line = dict(width=1),
                    line_color=LINE_COL,
                    marker={"size": 0.1},
                    text=r["id"],
                    customdata=[r["id"]]*5,
                    hoverinfo = "text+x+y",
                    showlegend=False))
            if r["id"] not in EXCLUDE_EXON_LABEL:
                fig.add_annotation(
                        x=(r["x0"]+r["x1"])/2, y=1,
                        text=r["id"],
                        showarrow=False,
                        textangle=0 if r["id"] in ROTATE_LABEL else -90,
                        font=dict(
                            family="Courier New, monospace",
                            size=11,
                            color=TEXT_COL
                            ))

    fig.add_trace(
        go.Scatter(
            x=[min(df["x0"]), max(df["x1"])],
            y=[1, 1],
            fill='toself', fillcolor="black",
            marker={"size": 0.1},
            showlegend=False,
            hoveron = "fills", # select where hover is active
            line_color="black")
        )


    fig.data = fig.data[::-1]

    for i, r in df.iterrows():        
        if r["type"] =="intron":
            fig.add_trace(go.Scatter(mode="markers", x=[(r["x0"]+r["x1"])/2], y=[2.25], marker_symbol="triangle-down",
                           marker_line_color=LINE_COL, marker_color=INTRON_COL,
                           marker_line_width=1, marker_size=20,
                           text=r["id"],
                           customdata=[r["id"]],
                           hoverinfo = "text",
                           showlegend=False))
            fig.add_annotation(
                    x=(r["x0"]+r["x1"])/2, y=3.2,
                    text=r["id"],
                    showarrow=False,
                    textangle=0 if r["id"] in ROTATE_LABEL else -90,
                    font=dict(
                        family="Courier New, monospace",
                        size=11,
                        color=TEXT_COL
                        ))
            
    newdict=dict()
    for i,x in enumerate(fig.data):
        newdict[i] = x

    global CURVENUMBER_DICT
    CURVENUMBER_DICT = newdict
    global FIGURE_OBJ
    FIGURE_OBJ = fig
    
    fig.update_xaxes(range=[0, max(df["x1"])+FLANK/INTRON_SCALE])
    fig.update_yaxes(range=[-0.1, 4])
    fig.update_xaxes(visible=False)

    fig.update_layout(plot_bgcolor="#F6F6F6", paper_bgcolor="LightSteelBlue",
            margin={"l": 0, "r": 0, "t": 0, "b": 0})
    fig.update_xaxes(fixedrange=True)
    fig.update_yaxes(fixedrange=True, visible=False)
    
    plot = dcc.Graph(figure=fig, id="cftr-plot", style={"height": "150px"})
    return plot


MINIPLOT = html.Div([
    dcc.Loading(id="cftr-mini-plot")
    ])


@app.callback(
    Output("selection-text", "children"),
    Input("cftr-plot", "hoverData"),
    Input("cftr-plot", "clickData"))
def update_selection_text(hoverData, clickData):
    #print("click" , clickData)
    #print("hover ", hoverData)

    TEXT_COL="dimgray"

    for i in range(len(FIGURE_OBJ.data)):
        curveId = -1 if clickData is None else clickData["points"][0]["curveNumber"]

        if i == curveId: 
            FIGURE_OBJ.data[i].fillcolor=TEXT_COL
    
    if clickData is not None:
        curve = CURVENUMBER_DICT[clickData["points"][0]["curveNumber"]]
        return html.P("click: " + curve["text"])
    curve = CURVENUMBER_DICT[hoverData["points"][0]["curveNumber"]]
    return html.P("hover: " + curve["text"])


LABEL = html.Div([
    html.Pre(id="selection-text")
    ])

@app.callback(
    Output('tbl_out', 'children'), 
    Input('tbl', 'active_cell'))
def update_graphs(active_cell):
    return str(active_cell) if active_cell else "Click the table"

@app.callback(
    Output("tbl", "style_data_conditional"),
    Input("tbl", "derived_viewport_selected_row_ids"),
)
def style_selected_rows(selRows):
    if selRows is None:
        return dash.no_update
    return [
        {"if": 
         {"filter_query": "{{id}} ={}".format(i)},
         "backgroundColor": "#fff5be",
         }
        for i in selRows
    ]

TABLE = dbc.Container([
    dbc.Label('Click a cell in the table:'),
    dash_table.DataTable(
        data=variant_data.to_dict('records'),
        columns=[{"name": c, "id": c} for c in variant_data.columns],
        id='tbl', 
        row_selectable='multi',
        page_size=10,
        fixed_rows={'headers': True},
        style_table={'height': 300}  # default is 500
    ),
    dbc.Alert(id='tbl_out'),
])

downloadButtonType = {"css": "btn btn-primary", "text":"Export", "type":"xlsx"}
clearFilterButtonType = {"css": "btn btn-outline-dark", "text":"Clear Filters"}

# Setup some columns 
# This is the same as if you were using tabulator directly in js 
# Notice the column with "editor": "input" - these cells can be edited
# See tabulator editor for options http://tabulator.info/docs/4.8/edit
columns = [
    { "title": "CHROM", "field": "chrom"},
    { "title": "POS", "field": "pos", "headerFilter":True, "hozAlign": "center" },
    { "title": "REF", "field": "ref",  },
    { "title": "ALT", "field": "alt", "hozAlign": "center" },
    ]

# Setup some data
data = []
for i, r in variant_data.iterrows():
    item={"id":i, "chrom": r["CHROM"], "pos": r["POS"],
          "ref": r["REF"], "alt": r["ALT"]}
    data.append(item)
    
    if len(data) > 100: break

TABULAR=html.Div([
    dash_tabulator.DashTabulator(
        id='tabulator',
        theme='tabulator_simple',  #optional
        columns=columns,
        data=data,
        options = { "selectable":"true", "pagination":"true", "paginationSize":5,
                   "paginationSizeSelector":"[10, 25, 50, 100, true]", "height":"300px"},
        downloadButtonType=downloadButtonType,
        clearFilterButtonType=clearFilterButtonType
    ),
    html.Div(id='output'),
    dcc.Interval(
                id='interval-component-iu',
                interval=1*10, # in milliseconds
                n_intervals=0,
                max_intervals=0
            )

])

@app.callback([ Output('tabulator', 'columns'), 
                Output('tabulator', 'data')],
                [Input('interval-component-iu', 'n_intervals')]) 
def initialize(val):
    return columns, data

@app.callback(Output('output', 'children'), 
    [Input('tabulator', 'rowClicked'),
    Input('tabulator', 'cellEdited'),
    Input('tabulator', 'dataChanged'), 
    Input('tabulator', 'dataFiltering'),
    Input('tabulator', 'dataFiltered')])
def display_output(row, cell, dataChanged, filters, dataFiltered):
    print(row)
    print(cell)
    print(dataChanged)
    print(filters)
    print(dataFiltered)
    return 'You have clicked row {} ; cell {}'.format(row, cell)

app.layout = html.Div(style={"padding": 20}, children=[
    HEADER,
    html.Hr(),
    dbc.Card([
                dbc.CardHeader([ html.H5("CFTR Plot") ]),
                dbc.CardBody([ MINIPLOT , PLOT_BUTTON ])
                ], style={"box-shadow": "0 4px 8px 0 rgba(0,0,0,0.2)"},
        id="plot-card"),   
    html.Hr(),
    LABEL,
    dbc.Card(
            id="table-card",
            children=[
                    dbc.CardHeader([ html.H5("CFTR Variants") ]),
                    dbc.CardBody([ TABLE ])],
                    className="w-75"
                    ),
    
    ])


if __name__ == '__main__':
    app.run_server(debug=True)
