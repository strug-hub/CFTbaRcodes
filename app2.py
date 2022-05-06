from dash import Dash, html
from dash.dependencies import Input, Output, State
import vcf
import dash_bootstrap_components as dbc
from dash import dash_table
import pandas as pd

import graph_component
SYMBOL=["♈","♉","♊","♋","♌","♍","♎","♏","♐","♑","♒","♓","⛎"]
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
    print(data)
    return data

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
        row_selectable="multi",
        fixed_rows={"headers": True},
        style_table={"height": 400})

    return selection_table

SELECTION_TABLE = dbc.Container(
        id="selection-container",
        children=[_selected_variant()],
        style={"width": "10vw"})

GRAPH_PLOT = dbc.Container(
        id="graph-plot-container",
        children=[graph_component.plot_graph()],
        style={"width": "80vw"})

GRAPH_TABLE = graph_component.define_haps()


app.layout = html.Div([SELECTION_TABLE, GRAPH_PLOT, GRAPH_TABLE])



app.clientside_callback(
    graph_component.CLIENTSIDE_UPDATE,
    Output("graph-plot", "figure"), 
    Input("graph-plot", "hoverData"),
    Input("graph-plot", "clickData"),
    Input("haplotype-table", "selected_row_ids"),
    State("graph-plot", "figure")
)


if __name__ == "__main__":
    app.run_server(debug=True, port=8069)

