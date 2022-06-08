import dash
from dash import Dash, html, dcc, dash_table
from dash.dependencies import Input, Output, State
import pandas as pd
import dash_bootstrap_components as dbc
import vcf
import random
import plotly.graph_objects as go
import plotly.colors
from PIL import ImageColor
import graph_component

SYMBOL=["🟥","🟧","🟨","🟩","🟦","🟪","🟫","⬛",
        "❌","⭕","❗","❓",
        "♈","♉","♊","♋","♌","♍","♎","♏","♐","♑","♒","♓","⛎"]
EMPTY_SYMBOL="◻️"

# bootstrap css
external_stylesheets = ['https://stackpath.bootstrapcdn.com/bootstrap/4.1.3/css/bootstrap.min.css']

CFTR_START=117479025 ; CFTR_END=117669665 ; FLANK=10000
CFTR_VIEW="chr7:" + str(CFTR_START-FLANK) + "-" + str(CFTR_END+FLANK)

VCF="quick_cftr_merge_chr7_117287120-117715971.vcf"
vcf_reader = vcf.Reader(filename=VCF)
vcf_records = [record for record in vcf_reader if record.num_called > 10]


LOGO = "assets/logo.svg"

LINE_COL="black"
EXON_COL_SELECTED="#afd4da"
EXON_COL="lightskyblue"
INTRON_COL_SELECTED="#dac8af"
INTRON_COL="coral"
UTR_COL_SELECTED="#caafda"
UTR_COL="plum"
TEXT_COL="dimgray"
PLOT_BG="#F6F6F6"
UNSELECTED_MARKER_COL="white"
DEFAULT_IDX=2
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
annotation_data = pd.concat([annotation_data, intron_data], sort=False)

# --------------------------------------
# READ VARIANT DATA
# --------------------------------------

def relative_positon(pos):
    prev_x = 0  ; prev_end = 0
    for i, r in annotation_data.iterrows():
        if r["type"] == "intron": continue
        if pos >= r["start"] and pos <= r["end"]: # in exon
            relativepos = pos - r["start"]
            return r["x0"]+relativepos
        elif pos < r["start"]: # in intron
            break
        prev_x=r["x1"] ; prev_end=r["end"]
    return prev_x + (pos - prev_end)/INTRON_SCALE

def get_group(pos):
    if pos < min(annotation_data["start"]):
        return("5' Flanking")
    if pos > max(annotation_data["end"]):
        return("3' Flanking")
    for i, r in annotation_data.iterrows():
        if pos >= r["start"] and pos <= r["end"]: # in exon
            return r["id"]

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
    "x": [relative_positon(v.POS) for v in vcf_records],
    "y": [1 for v in vcf_records],
    "rsid": [fake_rsid() for v in vcf_records],
    "type": [get_variant_type(v) for v in vcf_records],
    "af": [random.random()/2 for v in vcf_records],
    "called": [v.num_called for v in vcf_records],
    "nsamp": [len(v.samples) for v in vcf_records],
    "group":[get_group(v.POS) for v in vcf_records]
    }
variant_data = pd.DataFrame(data=variant_dict)

# --------------------------------------
# CREATE APP
# --------------------------------------

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "CFTbaRcodes"

HEADER = html.Div(children=[
    html.Img(src=LOGO, id="logo-image"),
    html.Div(children="A web application for understanding CFTR haplotypes.")
    ])

def get_empty_div(id, text="Nothing to show", height="400px"):
    return html.Div([text], id=id,
             style={"display": "flex",
                    "align-items": "center",
                    "justify-content":"center",
                    "color": TEXT_COL, 
                    "background":PLOT_BG,
                    "height":height})
def hide_mode():
    return {"displaylogo": False, "modeBarButtonsToRemove": ["select2d", "lasso2d", "toImage"]}


# --------------------------------------
# DRAW CFTR FIGURE
# --------------------------------------

def make_cftr_plot():
    
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
                    customdata=[r["id"]],
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
    
    fig.update_xaxes(range=[0, max(df["x1"])+FLANK/INTRON_SCALE])
    fig.update_yaxes(range=[-0.1, 4])

    fig.update_layout(plot_bgcolor=PLOT_BG,
            margin={"l": 0, "r": 0, "t": 0, "b": 0})
    fig.update_xaxes(fixedrange=True, visible=False)
    fig.update_yaxes(fixedrange=True, visible=False)
    
    return dcc.Graph(figure=fig, id="cftr-plot", style={"height": "150px"},
                     config=hide_mode())

def make_barcode_plot():
        
    fig = go.Figure()
    x=[r["x"] for _, r in variant_data.iterrows()]
    y=[0 for _, r in variant_data.iterrows()]
    t=[r["rsid"] for _, r in variant_data.iterrows()]
    fig.add_trace(go.Scatter(mode="markers", x=x, y=y, marker_symbol="line-ns",
       marker_line_color=LINE_COL, marker_color=LINE_COL,
       marker_line_width=0.1, marker_size=10, text=t, hoverinfo = "text",
       showlegend=False))
    
    fig.update_xaxes(range=[0, max(annotation_data["x1"])+FLANK/INTRON_SCALE])
    fig.update_yaxes(range=[-0.5, 0.5])
    fig.update_layout(plot_bgcolor=PLOT_BG,
            margin={"l": 0, "r": 0, "t": 0, "b": 0})
    fig.update_xaxes(fixedrange=True, visible=False)
    fig.update_yaxes(fixedrange=True, visible=False)
    return dcc.Graph(figure=fig, id="barcode-plot", style={"height": "50px"},
                      config=hide_mode())


CFTRPLOT = html.Div([ make_cftr_plot() ])
BARCODEPLOT = html.Div([ make_barcode_plot() ])

# --------------------------------------
# CLICK RESPONSE FUNCTIONS
# --------------------------------------

# Show selected element
@app.callback(
    Output("selection-text", "children"),
    Input("cftr-plot", "clickData"),
    Input("cftr-plot", "figure"))
def update_selection_text(click_data, figure):
    
    if click_data is not None:
        idx = click_data["points"][0]["curveNumber"]
        text=figure["data"][idx]["text"]
    else:
        text = "exon"
        
    color=UTR_COL
    if "exon" in text: 
        text="all exons"
        color=EXON_COL
    elif "intron" in text:
        color=INTRON_COL
    
    return html.Div(style={"display": "flex", "align-items": "center"}, 
                    children=[html.Div("Selected: "),
                              html.Div(text,
                                       style={"background": color,
                                              "padding": "4px",
                                              "border-radius": "4px",
                                              "border": "2px solid "+LINE_COL,
                                              "color": TEXT_COL,
                                              "font-family": "Courier New, monospace"
                                              })])

# Change outline for selected element(s) in the gene box
@app.callback(
    Output("cftr-plot", "figure"),
    Input("cftr-plot", "clickData"),
    Input("cftr-plot", "figure"))
def update_cftr_plot_selection(click_data, figure):
    
    if click_data is not None:
        idx = click_data["points"][0]["curveNumber"]
    else:
        idx=DEFAULT_IDX
    
    trace = figure["data"][idx]
    
    for i in range(len(figure["data"])):
        if "customdata" in figure["data"][i]:
            w=0.1

            clickId=trace["text"]
            if figure["data"][i]["customdata"][0] == clickId or \
                ("exon" in clickId and "exon" in figure["data"][i]["customdata"][0]):
                    w=1.5
                
            if "line" in figure["data"][i]:
                figure["data"][i]["line"]["width"] = w
            if "marker" in figure["data"][i]:
                if "line" not in figure["data"][i]["marker"]:
                    figure["data"][i]["marker"]["line"] = dict()
                figure["data"][i]["marker"]["line"]["width"] = w
        
        
    return figure

# --------------------------------------
# VARIANT TABLE
# --------------------------------------

TABLE_COLUMNS = [
    { "id": "rsid", "name": "rsid" },
    { "id": "CHROM", "name": "chrom" },
    { "id": "POS", "name": "pos" },
    { "id": "REF", "name": "ref" },
    { "id": "ALT", "name": "alt" },
    { "id": "type", "name": "type"},
    { "id": "af", "name": "allele freq"},
    { "id": "group", "name": "location"}
    ]

def empty_variant_table():
    empty_div = get_empty_div("variant-table-empty", text="No variants to display", height="250px")
    empty_div.children.append(dash_table.DataTable(id="variant-table"))
    return empty_div

# Update table with variants from selected elements
@app.callback(
    Output("variant-table-container", "children"),
    Input("cftr-plot", "clickData"),
    Input("cftr-plot", "figure"),
    Input("selected-variants", "data"))
def update_table_target(click_data, figure, selected_data):
    
    if selected_data is None: selected_data = {"selected":[]}
    selected_ids = selected_data["selected"]
        
    if click_data is not None: 
        idx = click_data["points"][0]["curveNumber"]
    else:
        idx = DEFAULT_IDX
        
    trace = figure["data"][idx]
    target_group=trace["text"]
    
    if "exon" in target_group:
        vdf = variant_data[variant_data["group"].str.contains("exon")]
    else:
        vdf = variant_data[variant_data["group"] == target_group]
    
    if len(vdf) == 0 : return empty_variant_table()
    
    # convert "id" to row index
    i = 0 ; selected_rows = []
    for _, r in vdf.iterrows():
        if r["id"] in selected_ids : selected_rows.append(i)
        i+=1

    variant_table = dash_table.DataTable(
        data=vdf.to_dict("records"),
        id="variant-table",
        columns=TABLE_COLUMNS,
        row_selectable='multi',
        selected_rows=selected_rows,
        page_size=50,
        fixed_rows={'headers': True},
        style_table={'height': 250}
        )
    
    condition = []
    for i in selected_ids:
        condition.append({"if": {"filter_query": "{{id}} ={}".format(i)},
                          "backgroundColor": "#fff5be" })
        
    variant_table.style_data_conditional = condition

    return variant_table


TABLE = dbc.Container([
    html.Pre(id="selection-text"),
    html.Div(id="variant-table-container",
             children=[empty_variant_table()])], 
    style={"width": "48vw"})


# --------------------------------------
# COMPONENT PLOT
# --------------------------------------

def empty_component_plot():
    empty_div = get_empty_div("component-plot-empty", text="No variants to display", height="300px")
    empty_div.children.append(dcc.Graph(id="component-plot", style={"display": "none"}))
    return empty_div

@app.callback(
    Output("component-plot-container", "children"),
    Input("cftr-plot", "clickData"),
    Input("cftr-plot", "figure"))
def make_cftr_component_plot(click_data, figure):
    
    COLLAPSE_PC=0.025

    if click_data is not None: 
        idx = click_data["points"][0]["curveNumber"]
    else:
        idx = DEFAULT_IDX
        
    trace = figure["data"][idx]
    target_group=trace["text"]

    if "exon" in target_group:
        vdf = variant_data[variant_data["group"].str.contains("exon")]
    else:
        vdf = variant_data[variant_data["group"] == target_group]

    if len(vdf) == 0 : return empty_component_plot()

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=[min(annotation_data["start"]), max(annotation_data["end"])],
            y=[0.15, 0.15],
            fill='toself', fillcolor="black",
            marker={"size": 0.1},
            showlegend=False,
            hoveron = "fills",
            line_color="black")
        )

    for i, r in annotation_data.iterrows():
        if not "intron" in r["type"]:
            col = UTR_COL if r["type"] == "UTR" else EXON_COL
            #adding an extra marker to assist in hover detection
            fig.add_trace(go.Scatter(mode="markers", x=[(r["start"]+ r["end"])/2],
                                     y=[0.15], marker_symbol="line-ns",
                                     marker_line_color=col, marker_color=col,
                                     marker_line_width=1, text=r["id"], 
                                     hoverinfo = "text", showlegend=False))

            fig.add_trace(go.Scatter(
                    x=[r["start"], r["start"], r["end"], r["end"], r["start"]],
                    y=[0.1, 0.2, 0.2, 0.1, 0.1],
                    fill='toself',
                    hoveron = 'fills',
                    line = dict(width=1),
                    line_color=col,
                    fillcolor=col,
                    marker={"size": 0.1},
                    text=r["id"],
                    customdata=[r["id"]],
                    hoverinfo = "text+x+y",
                    showlegend=False))

    #fig.update_layout(hovermode="x")
    collapse_range=(max(vdf["POS"])-min(vdf["POS"]))*COLLAPSE_PC
    collapse=[] ; collapse_pos = -collapse_range*2 ; c = 0
    for i, r in vdf.iterrows():
        if r["POS"] - collapse_pos > collapse_range:
            collapse_pos = r["POS"]
            c+=1
        collapse.append(c)
        
    vdf["collapse"] = collapse
    x = [r["POS"] for i, r in vdf.iterrows()]
    y = [0.3 for i, r in vdf.iterrows()]
    t = [r["rsid"] for i, r in vdf.iterrows()]
    fig.add_trace(go.Scatter(mode="markers", x=x, y=y, marker_symbol="line-ns",
                   marker_line_color="black", marker_color="black",
                   marker_line_width=1, marker_size=10, text=t, hoverinfo = "text",
                   showlegend=False))

    shapedict = {VARIANT_DEL: "diamond", VARIANT_INS: "diamond",
                 VARIANT_SNP: "circle", VARIANT_ETC: "circle"}
    colorscaledict = {VARIANT_DEL: "oranges", VARIANT_INS: "oranges",
                 VARIANT_SNP: "blues", VARIANT_ETC: "blues"}

    for i in set(vdf["collapse"]):
        col = LINE_COL  
        y=0.4
        for j, r in vdf[vdf["collapse"] == i].iterrows():
            
            shape = shapedict[r["type"]]
            colorscale = colorscaledict[r["type"]]

            fig.add_trace(go.Scatter(
                mode="markers", x=[r["POS"]], y=[y], marker_symbol=shape,
                marker_line_color=LINE_COL, marker_color=get_color(colorscale, r["af"]),
                marker_line_width=0.5, marker_size=8, customdata=[r["id"]],
                showlegend=False, text=[r["rsid"]], hoverinfo = "text"))
            y+=0.08

    
    xmin = min(annotation_data[annotation_data["id"] == target_group]["start"])
    xmax = min(annotation_data[annotation_data["id"] == target_group]["end"])
    if "exon" in target_group:
        xmin=min(annotation_data["start"])
        xmax=max(annotation_data["end"])
    fig.update_yaxes(range=[0, 2])
    fig.update_xaxes(range=[xmin-(xmax-xmin)*0.05, xmax+(xmax-xmin)*0.05 ])

    fig.update_layout(plot_bgcolor="#F6F6F6",
            margin={"l": 0, "r": 0, "t": 0, "b": 0},
            height=300)
    fig.update_xaxes(fixedrange=True, title_text="GRCh38 Coordinate")
    fig.update_yaxes(fixedrange=True, visible=False)
    return dcc.Graph(figure=fig, id="component-plot")


COMPONENT_PLOT = dbc.Container([
        html.Div(id="component-plot-container", 
                 children=[empty_component_plot()])
        ], style={"width": "48vw"})

@app.callback(
    Output("component-plot", "figure"),
    Input("selected-variants", "data"),
    State("component-plot", "figure"))
def update_component_selection(selected_data, figure):

    if figure is None: return dash.no_update

    if selected_data is None: selected_data = {"selected":[]}
    selected = selected_data["selected"]

    for curve in figure["data"]:
        if "marker" in curve and "customdata" in curve:
            id=curve["customdata"][0]
            if not "line" in curve["marker"]:
                curve["marker"]["line"]=dict()
            curve["marker"]["line"]["width"] = 2 if id in selected else 0.1

    return figure

# --------------------------------------
# VARIANT SELECTION
# --------------------------------------
def empty_variant_table_selection():
    empty_div = get_empty_div("variant-table-selection-empty", text="No variants to display", height="250px")
    empty_div.children.append(dash_table.DataTable(id="variant-table-selection"))
    return empty_div

@app.callback(
    Output("variant-table-selection-container", "children"),
    Input("selected-variants", "data"))
def add_selected_variant(selected_data):
    
    if selected_data is None: selected_data = {"selected":[]}
    selected_ids = selected_data["selected"]

    if selected_ids is None or len(selected_ids) == 0: return empty_variant_table_selection()
    vdf = variant_data[variant_data["id"].isin(selected_ids)]
    
    vdf["symbol"] = [EMPTY_SYMBOL for _ in selected_ids]

    dropdown={"symbol": { "options": [ {"label": EMPTY_SYMBOL, "value": EMPTY_SYMBOL} ] + [
                         {"label": m, "value": m} for m in SYMBOL ] } }
    
    selection_table = dash_table.DataTable(
        data=vdf.to_dict("records"),
        id="variant-table-selection",
        columns=TABLE_COLUMNS + [{"id": "symbol", "name": "symbol", "presentation": "dropdown"}],
        editable=True,
        row_deletable=True,
        dropdown=dropdown,
        page_size=50,
        fixed_rows={'headers': True},
        style_table={'height': 250})

    return selection_table


@app.callback(
    Output("variant-table-selection", "data"),
    Input("variant-table-selection", "data"))
def update_table_dropdown(selection_table_data):
    if selection_table_data is None: return dash.no_update
    for x in selection_table_data:
        if x["symbol"] is None: x["symbol"] = EMPTY_SYMBOL
    return selection_table_data

@app.callback(
    Output("selected-variants", "data"),
    Input("variant-table", "selected_row_ids"),
    Input("variant-table", "data"),
    Input("variant-table-selection", "data"),
    Input("component-plot", "clickData"),
    Input("predefined-dropdown", "value"),
    State("selected-variants", "data"))
def variant_row_change(selected_ids, table_data, selection_table_data, 
                       variant_plot_click, predefined_dropdown_value, store_data):
        
    store_data = store_data or {"selected": []}
    selected = [id for id in store_data["selected"]]
    to_rm = set()

    ctx = dash.callback_context

    # handle predefined dropdown selection
    if ctx.triggered and "predefined-dropdown" in ctx.triggered[0]["prop_id"]:
        print("Loading predefined: " + str(predefined_dropdown_value))
            
        if predefined_dropdown_value == PREDEFINED0:
            return dash.no_update            
        if predefined_dropdown_value == PREDEFINED1:
            store_data = {"selected": [172, 173, 174, 175, 178, 183]}
            
        return store_data

    # handle variant clicked on component plot
    if ctx.triggered and "component-plot" in ctx.triggered[0]["prop_id"]:
        if variant_plot_click is not None:
            for clicked in variant_plot_click["points"]:
                if "customdata" in clicked:
                    clicked_id = clicked["customdata"]
                    if clicked_id in selected and len(variant_plot_click["points"])==1:
                        selected.remove(clicked_id)
                    else:
                        selected.append(clicked_id)
        store_data["selected"] = selected
        return store_data

    
    # handle variant removed from selection table
    if selection_table_data is not None:
        table_ids = [x["id"] for x in selection_table_data]
        for id in selected:
            if id not in table_ids: to_rm.add(id)

    # handle variant [un]selected from variant table
    if selected_ids is not None and table_data is not None:
        table_ids = [x["id"] for x in table_data]
        selected = store_data["selected"] + selected_ids
        selected = list(set(selected))
        for id in selected:
            if id in table_ids and id not in selected_ids:
                to_rm.add(id)

    for id in to_rm: 
        selected.remove(id)

    if store_data["selected"] == selected:
        return dash.no_update

    print("selected: ", selected)
    store_data["selected"] = selected
    return store_data

TABLE_SELECTION = dbc.Container(
    id="variant-table-selection-container",
    children=[ empty_variant_table_selection() ])

PREDEFINED0="Select Predefined Variants"
PREDEFINED1="Predefined 1"
PREDEFINED2="Predefined 2"
PREDEFINED3="Predefined 3"

PREDEFINED_SELECTIONS = html.Div([
    dcc.Dropdown([PREDEFINED0, PREDEFINED1, PREDEFINED2, PREDEFINED3], PREDEFINED0, id="predefined-dropdown"),
    html.Div(id="predefined-dropdown-container")
])

# --------------------------------------
# GRAPH PLOT
# --------------------------------------

@app.callback(
    Output("graph-plot-container", "children"), 
    Input("selected-variants", "data"))
def graph_replot(store_data):
    records = [vcf_records[i] for i in store_data["selected"]]
    plot = graph_component.plot_graph(records)
    return [dcc.Graph(figure=plot, id="graph-plot", style={'width': '100%'})]

@app.callback(
    Output("graph-plot", "figure"),
    Input("variant-table-selection", "data"),
    Input("haplotype-table", "selected_rows"),
    State("haplotype-table", "data"),
    State("graph-plot", "figure"),
    State("selected-variants", "data"))
def update_graph_figure(selection_data, haplotype_selection, haplotype_data, figure, store_data):
    
    ctx = dash.callback_context

    if ctx.triggered and "variant-table-selection" in ctx.triggered[0]["prop_id"]:
        symbols = [d["symbol"] for d in selection_data] 
        figure = graph_component.update_symbols(figure, symbols)
    
    if ctx.triggered and "haplotype-table" in ctx.triggered[0]["prop_id"]:

        if haplotype_selection is None: return dash.no_update
        if haplotype_selection[0] is None: return dash.no_update
        
        target_variants = dict()
        target = haplotype_data[haplotype_selection[0]]
        print(target)
        for i,row in enumerate(selection_data):
            if row["symbol"] in target:
                target_variants[i] = target[row["symbol"]]
        
        haplotypes = graph_component.get_haps(figure)
        haplogroup = []
        for haplotype in haplotypes:
            
            check = [ target_variants[vid] != haplotype[vid] for vid in target_variants ]
            if sum(check) == 0:
                haplogroup.append(haplotype)

        figure = graph_component.draw_haplotype(figure, haplogroup)
    return figure 

@app.callback(
    Output("haplotype-container", "children"),
    Input("variant-table-selection", "data"),
    State("graph-plot", "figure"))
def update_haplotype_table(selection_data, figure):
    haplotypes = graph_component.get_haps(figure)
    
    selection = [i for i,s in enumerate(selection_data) if s["symbol"] != EMPTY_SYMBOL]
    symbol = [s["symbol"] for s in selection_data if s["symbol"] != EMPTY_SYMBOL]

    groups = dict()
    for i,haplotype in enumerate(haplotypes):
        group = tuple(haplotype[s] for s in selection)
        if not group in groups: groups[group] = []
        groups[group].append(i)
    
    def fake_rsid():
        return "rs"+ str(round(random.random()*1e6))
    rsids= [ fake_rsid() for s in selection ]
    
    total = sum([len(groups[g]) for g in groups])
    rows = []
    for group in groups:
        size = len(groups[group])
        r = {"freq" : str(round(size*100/total, 2)) + "%" , "n" : size }
        for i,gt in enumerate(group): r[symbol[i]] = gt
        rows.append(r)

    tooltip = { SYMBOL[i]: {"value": rsid, "use_with": "both"} for i,rsid in enumerate(rsids) }
    #columns=[ {"name": i, "id": i, "selectable": False} for i in range(len(rows))]
    n="#Haps"
    df = pd.DataFrame(dict( 
        [(n, [row["n"] for row in rows])] + \
        [("Frequency", [row["freq"] for row in rows])] + \
        [(s, [row[s] for row in rows]) for s in symbol]

    ))
    
    df.sort_values(by=[n], ascending=False, inplace=True)
    before_len = len(df)
    df = df[df[n] > 3]
    after_len = len(df)

    return [dash_table.DataTable(df.to_dict("records"), id="haplotype-table",
                                row_selectable="single",
                                tooltip=tooltip,
                                tooltip_delay=0),
            html.Div(str(before_len-after_len)+" rare haplotypes hidden")
            ]

HAPLOTYPE_TABLE = dbc.Container(
        id="haplotype-container",
        children=[get_empty_div("haplotype-table", height=250)],
        style={"width": "30vw"})

GRAPH_PLOT = dbc.Container(
        id="graph-plot-container",
        children=[],
        style={"max-width": "100%"})


# --------------------------------------
# APP LAYOUT
# --------------------------------------

cardstyle={"box-shadow": "0 4px 8px 0 rgba(0,0,0,0.2)"}   

app.layout = html.Div(style={"padding": 20}, children=[
    dcc.Store(id="selected-variants", storage_type="memory"),
    HEADER, PREDEFINED_SELECTIONS,
    html.Hr(),
    dbc.Card(
        children = [ dbc.CardHeader([ html.H5("CFTR Plot") ]),
                    dbc.CardBody([ BARCODEPLOT, CFTRPLOT ])],
        style=cardstyle,
        id="cftr-card"),   
    html.Br(),
    dbc.Card(
        children=[dbc.CardHeader([ html.H5("CFTR Variants") ]),
                  dbc.CardBody(html.Div(style={"display": "flex", "align-content":"stretch"},
                                        children=[ TABLE, COMPONENT_PLOT ]))],
        style=cardstyle,
        id="variant-card"),   
    html.Br(),
    dbc.Card(
        children=[dbc.CardHeader([ html.H5("Selected Variants") ]),
                  dbc.CardBody(html.Div(style={"display": "flex", "align-content":"stretch"},
                                        children=[ TABLE_SELECTION, HAPLOTYPE_TABLE ]))],
        style=cardstyle,
        id="selected-card"),
    html.Br(),
    dbc.Card(
        children=[dbc.CardHeader([ html.H5("Haplotype Graph") ]),
                  dbc.CardBody(children=[ GRAPH_PLOT ])],
        style=cardstyle,
        id="graph-card"),

    ])


# Interpolating from a plotly color scale

def get_color(colorscale_name, loc):
    from _plotly_utils.basevalidators import ColorscaleValidator
    # first parameter: Name of the property being validated
    # second parameter: a string, doesn't really matter in our use case
    cv = ColorscaleValidator("colorscale", "")
    # colorscale will be a list of lists: [[loc1, "rgb1"], [loc2, "rgb2"], ...] 
    colorscale = cv.validate_coerce(colorscale_name)
    
    if hasattr(loc, "__iter__"):
        return [get_continuous_color(colorscale, x) for x in loc]
    return get_continuous_color(colorscale, loc)
        

def get_continuous_color(colorscale, intermed):

    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    hex_to_rgb = lambda c: "rgb" + str(ImageColor.getcolor(c, "RGB"))

    if intermed <= 0 or len(colorscale) == 1:
        c = colorscale[0][1]
        return c if c[0] != "#" else hex_to_rgb(c)
    if intermed >= 1:
        c = colorscale[-1][1]
        return c if c[0] != "#" else hex_to_rgb(c)

    for cutoff, color in colorscale:
        if intermed > cutoff:
            low_cutoff, low_color = cutoff, color
        else:
            high_cutoff, high_color = cutoff, color
            break

    if (low_color[0] == "#") or (high_color[0] == "#"):
        # some color scale names (such as cividis) returns:
        # [[loc1, "hex1"], [loc2, "hex2"], ...]
        low_color = hex_to_rgb(low_color)
        high_color = hex_to_rgb(high_color)

    return plotly.colors.find_intermediate_color(
        lowcolor=low_color,
        highcolor=high_color,
        intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
        colortype="rgb",
    )

if __name__ == '__main__':
    app.run_server(debug=True)