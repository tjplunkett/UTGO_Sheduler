"""
scheduler.py 

Author: Thomas Plunkett

Date: 29/06/23

Purpose: Creates an interactive dashboard to assist in observations scheduling 
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

# Create table for microlensing events
ms_df = construct_final()

app = Dash(__name__)

app.layout = html.Div(id='main', children = [
    #Title
    html.H1(children='Greenhill Observatory Target Scheduler', style={'textAlign':'center'}), html.Br(),
    
    # Microlensing event table
    dash_table.DataTable(id='ulensTable', data=ms_df.to_dict('records'), columns = [{"name": i, "id": i} for i in ms_df.columns], editable=True, style_data={'whiteSpace': 'normal','height': 'auto'}, style_header={'fontWeight': 'bold', 'color':'#FFFFFF', 'backgroundColor':'#008B8B'}),
    
    html.Br(),
    
    html.Button('Add Target', id='target_ad', n_clicks=0, style={'width':'200 px', 'height':'200 px','color':'#FFFFFF', 'backgroundColor':'#008B8B'}),
    
    html.Button('Generate schedule', id='sched_button', n_clicks=0, style={'width':'200 px', 'height':'200 px','color':'#FFFFFF', 'backgroundColor':'#008B8B'}),   
    
    html.Br(), 
    
    html.Div(id = 'sched_div'),
])

@callback(
    Output('ulensTable', 'data'),
    Input('target_ad', 'n_clicks'),
    State('ulensTable', 'data'),
    State('ulensTable', 'columns'))
def add_row(n_clicks, rows, columns):
    if n_clicks > 0:
        rows.append({c['id']: '' for c in columns})
    return rows

@app.callback(
    Output('sched_div', 'children'),
    Input('sched_button', 'n_clicks'),
    State('ulensTable', 'data'),
    State('ulensTable', 'columns'))
def generate_schedule(n_clicks, rows, columns):
    if n_clicks > 0:
        master_df = pd.DataFrame(rows, columns=[c['id'] for c in columns])
        sched = construct_schedule(master_df)
        return dash_table.DataTable(id='Schedule', data=sched.to_dict('records'), style_data={'whiteSpace': 'normal','height': 'auto'}, style_header={'fontWeight': 'bold', 'color':'#FFFFFF', 'backgroundColor':'#008B8B'})
        

if __name__ == '__main__':
    app.run_server(debug=True)