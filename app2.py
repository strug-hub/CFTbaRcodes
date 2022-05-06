from dash import Dash, html
import dash
from dash.dependencies import Input, Output, State
import vcf
import dash_bootstrap_components as dbc
from dash import dash_table
import pandas as pd
import random

import graph_component
SYMBOL=["🟥","🟧","🟨","🟩","🟦","🟪","🟫","⬛","❌","⭕","❗","❓","♈","♉","♊","♋","♌","♍","♎","♏","♐","♑","♒","♓","⛎"]
EMPTY_SYMBOL="◻️"


app = Dash(__name__)

VCF="quick_cftr_merge_chr7_117287120-117715971.vcf"
vcf_reader = vcf.Reader(filename=VCF)
vcf_records = [record for record in vcf_reader if record.num_called > 10]

# --------------------------------------
# GRAPH PLOT
# --------------------------------------

@app.callback(
    Output("selection-table", "data"),
    Input("selection-table", "data"))
def update_table_dropdown(data):
    for x in data:
        if x["symbol"] is None: x["symbol"] = EMPTY_SYMBOL
    return data

@app.callback(
    Output("graph-plot", "figure"),
    Input("selection-table", "data"),
    Input("haplotype-table", "selected_rows"),
    State("haplotype-table", "data"),
    State("graph-plot", "figure"))
def update_graph_figure(selection_data, haplotype_selection, haplotype_data, figure):
    
    ctx = dash.callback_context

    if ctx.triggered and "selection-table" in ctx.triggered[0]["prop_id"]:
        symbols = [d["symbol"] for d in selection_data] 
        figure = graph_component.update_symbols(figure, symbols)
    
    if ctx.triggered and "haplotype-table" in ctx.triggered[0]["prop_id"]:
        if haplotype_selection is None: return dash.no_update
        if haplotype_selection[0] is None: return dash.no_update

        target_variants = dict()
        target = haplotype_data[haplotype_selection[0]]
        for row in selection_data:
            if row["symbol"] in target:
                target_variants[row["idx"]] = target[row["symbol"]]
        
        haplotypes = graph_component.get_haps(figure)
        
        haplogroup = []
        for haplotype in haplotypes:
            
            check = [ target_variants[vid] != haplotype[vid] for vid in target_variants ]
            if sum(check) == 0:
                haplogroup.append(haplotype)
        
        figure = graph_component.draw_haplotype(figure, haplogroup)
        
        
    return figure 

def _selected_variant(vcf_records=None):
    
    if vcf_records is None or len(vcf_records) < 1:        
        VCF="cftr.vcf.gz"
        vcf_reader = vcf.Reader(filename=VCF)
        vcf_records = [record for record in vcf_reader]
        
    df = pd.DataFrame(dict([
        ("symbol", [EMPTY_SYMBOL for _ in vcf_records]),
        ("idx", [i for i,_ in enumerate(vcf_records)])
    ]))
    columns=[
        {"id": "symbol", "name": "Symbol", "presentation": "dropdown"},
        {"id": "idx", "name": "ID"}   ]
    dropdown={"symbol": { "options": [ {"label": EMPTY_SYMBOL, "value": EMPTY_SYMBOL} ] + [
                         {"label": m, "value": m} for m in SYMBOL ] } }

    selection_table = dash_table.DataTable(
        data=df.to_dict("records"),
        id="selection-table",
        columns=columns,
        dropdown=dropdown,
        editable=True,
        row_deletable=True,
        fixed_rows={"headers": True},
        style_table={"height": 400})

    return selection_table

@app.callback(
    Output("haplotype-container", "children"),
    Input("selection-table", "data"),
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
        r = {"freq" : str(round(size*100/total, 2)) + "%" }
        for i,gt in enumerate(group): r[symbol[i]] = gt
        rows.append(r)

    tooltip = { SYMBOL[i]: {"value": rsid, "use_with": "both"} for i,rsid in enumerate(rsids) }
    #columns=[ {"name": i, "id": i, "selectable": False} for i in range(len(rows))]

    df = pd.DataFrame(dict( 
        [("freq", [row["freq"] for row in rows])] + \
        [(s, [row[s] for row in rows]) for s in symbol]
    ))
    return dash_table.DataTable(df.to_dict("records"), id="haplotype-table",
                                row_selectable="single",
                                tooltip=tooltip,
                                tooltip_delay=0)

TEXT_COL="dimgray"
PLOT_BG="#F6F6F6"
def get_empty_div(id, text="Nothing to show", height="400px"):
    return html.Div([text], id=id,
             style={"display": "flex",
                    "align-items": "center",
                    "justify-content":"center",
                    "color": TEXT_COL, 
                    "background":PLOT_BG,
                    "height":height})


SELECTION_TABLE = dbc.Container(
        id="selection-container",
        children=[_selected_variant()],
        style={"width": "10vw"})

GRAPH_PLOT = dbc.Container(
        id="graph-plot-container",
        children=[graph_component.plot_graph()],
        style={"width": "80vw"})

HAPLOTYPE_TABLE = dbc.Container(
        id="haplotype-container",
        children=[get_empty_div("haplotype-table")],
        style={"width": "20vw"})

cardstyle={"box-shadow": "0 4px 8px 0 rgba(0,0,0,0.2)"}   

app.layout = html.Div(    
    dbc.Card(
        children=[dbc.CardHeader([ html.H5("Graph") ]),
                  dbc.CardBody(html.Div(style={"display": "flex", "align-content":"stretch"},
                                        children=[SELECTION_TABLE, GRAPH_PLOT, HAPLOTYPE_TABLE]))],
        style=cardstyle,
        id="graph-card"))


if __name__ == "__main__":
    app.run_server(debug=True, port=8069)

