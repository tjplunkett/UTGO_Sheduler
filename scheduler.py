"""
scheduler.py 

Author: Thomas Plunkett

Date: 24/03/24

Purpose: Creates an interactive dashboard to assist in observations scheduling. 
"""
# Import necessary packages 
from dash import Dash, html, dcc, callback, Output, Input, dash_table, State
from dash.exceptions import PreventUpdate
import plotly.express as px
from datetime import datetime, timedelta, date
from scanner import *

# Initialize the location/observatory
bisdee_tier = EarthLocation.from_geodetic(lon = 147.3*u.deg, lat = -42.43*u.deg, height = 646*u.m)
GHO = Observer(location=bisdee_tier, name='Greenhill Observatory', timezone='Australia/Tasmania')

# Read in the table of events from Google sheet
ms_df = read_targets()

app = Dash(__name__)

app.layout = html.Div(id='main', children = [
    #Title
    html.H1(children='Greenhill Observatory Target Scheduler', style={'textAlign':'center'}), html.Br(),
    
    # Event table
    dash_table.DataTable(id='TargetTable', data=ms_df.to_dict('records'), columns = [{"name": i, "id": i} for i in ms_df.columns], editable=True, style_data={'whiteSpace': 'normal','height': 'auto'}, style_header={'fontWeight': 'bold', 'color':'#FFFFFF', 'backgroundColor':'#008B8B'}),
    
    html.Br(),
    
    html.Button('Add Target', id='target_ad', n_clicks=0, style={'width':'200 px', 'height':'200 px','color':'#FFFFFF', 'backgroundColor':'#008B8B'}),
    
    #html.Button('Read CSV', id='csv_add', n_clicks=0, style={'width':'200 px', 'height':'200 px','color':'#FFFFFF', 'backgroundColor':'#008B8B'}),
    
    html.Button('Generate schedule', id='sched_button', n_clicks=0, style={'width':'200 px', 'height':'200 px','color':'#FFFFFF', 'backgroundColor':'#008B8B'}),   
    
    html.Br(), 
    
    html.Div(id = 'sched_div'),
    
    html.Br(),
    
    html.Div(id='skychart'),
])

@callback(
    Output('TargetTable', 'data'),
    Input('target_ad', 'n_clicks'),
    State('TargetTable', 'data'),
    State('TargetTable', 'columns'))
def add_row(n_clicks, rows, columns):
    if n_clicks > 0:
        rows.append({c['id']: '' for c in columns})
    return rows

@app.callback(
    Output('sched_div', 'children'),
    Input('sched_button', 'n_clicks'),
    State('TargetTable', 'data'),
    State('TargetTable', 'columns'))
def generate_schedule(n_clicks, rows, columns):
    if n_clicks > 0:
        master_df = pd.DataFrame(rows, columns=[c['id'] for c in columns])
        master_df = master_df.astype({'PRIORITY': 'float64', 'EXP. TIME [s]': 'float64', 'IMAGES': 'int64', 'REPEATS': 'int64'})
        sched = construct_schedule(master_df)
        return html.Br(), dash_table.DataTable(id='Schedule', data=sched.round(4).to_dict('records'), style_data={'whiteSpace': 'normal','height': 'auto'}, style_header={'fontWeight': 'bold', 'color':'#FFFFFF', 'backgroundColor':'#008B8B'}), html.Button('Send to Telescope', id='tele_button', n_clicks=0, style={'width':'200 px', 'height':'200 px','color':'#FFFFFF', 'backgroundColor':'#008B8B'})
    
@app.callback(
    Output('skychart', 'children'),
    Input('sched_button', 'n_clicks'),
    State('TargetTable', 'data'),
    State('TargetTable', 'columns'))
def generate_sky(n_clicks, rows, columns):
    if n_clicks > 0:
        master_df = pd.DataFrame(rows, columns=[c['id'] for c in columns])
        master_df = master_df.astype({'PRIORITY': 'float64', 'EXP. TIME [s]': 'float64', 'IMAGES': 'int64', 'REPEATS': 'int64'})
        sched = construct_schedule(master_df)
        sky_polar = construct_sky(sched)
        
        return dcc.Graph(id='polar', figure=sky_polar, style={'width':'80vw', 'margin-left': '10vw', 'margin-right': '10vw'}), html.Br()
        
if __name__ == '__main__':
    app.run_server(debug=False, port = 8080)