import dash
from dash import dcc, html
import plotly.graph_objs as go
import pandas as pd
import colorlover as cl
import argparse
import colorsys

# Initialize the parser
parser = argparse.ArgumentParser(description="Interactive BLAST Results Visualization")

# Add command-line arguments
parser.add_argument("blast_results_file", help="Path to the BLAST results file")
parser.add_argument("--bannotat", help="Path to the BLAST annotations file")
parser.add_argument("--sc", dest="show_containing", action="store_true", help="Include containing hits in the visualization")

# Parse the command-line arguments
args = parser.parse_args()

blast_results_file = args.blast_results_file
blast_annotations_file = args.bannotat
show_containing = args.show_containing if hasattr(args, 'show_containing') and args.show_containing else False

# Read tabular BLAST results into a DataFrame
blast_df = pd.read_csv(blast_results_file, sep='\t', header=None)
blast_df.columns = ["query_name", "subject_name", "percent_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"]

# Initialize the annotations DataFrame if the file path is provided
annotations_df = pd.DataFrame
if blast_annotations_file:
    annotations_df = pd.read_csv(blast_annotations_file, sep='\t', header=None)
    annotations_df.columns = ["query_name", "subject_name", "percent_identity", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"]


# Get unique queries and subjects
queries = blast_df['query_name'].unique()
subjects = blast_df['subject_name'].unique()

# Define the number of colors you want
num_colors = len(queries)

# Get a color scale from colorlover, such as "Set3"
color_scale = cl.scales['9']['qual']['Set1'][:num_colors]

# Create a dictionary to map queries to colors
query_colors = {query: color for query, color in zip(queries, color_scale)}

# Initialize the Dash app
app = dash.Dash(__name__)

# Define the layout of the web application
app.layout = html.Div([
    dcc.Tabs([
        dcc.Tab(label=subject, children=[
            dcc.Graph(
                id=f'blast-visualization-{subject}',
                config={'displayModeBar': True},  # Show the mode bar for zooming
                figure={'data': [], 'layout': {'xaxis.autorange': True}},  # Initialize with an empty figure
            )
        ]) for subject in subjects
    ]),
])

# Create a callback function to update the graph based on relayoutData
@app.callback(
    [dash.dependencies.Output(f'blast-visualization-{subject}', 'figure') for subject in subjects],
    [dash.dependencies.Input(f'blast-visualization-{subject}', 'relayoutData') for subject in subjects]
)
def update_blast_visualization(*args):
    figures = []

    for subject in subjects:
        figure = go.Figure()
         # Keep track of encountered annotation names
        encountered_names = set()

     
        # Create a list to store legend entries dynamicall
        for index, row in blast_df[blast_df['subject_name'] == subject].iterrows():
            query_name = row['query_name']
            hsp_start, hsp_end = row['s_start'], row['s_end']
            query_start, query_end = row['q_start'], row['q_end']
            
            color = query_colors[query_name]

            showLegend = True
            if query_name in encountered_names:
                 showLegend = False
            
            # Check for containment with the previous trace
            
            y_shift = 0
            isContaining = any((trace.x[0] < hsp_end < trace.x[2]) or (trace.x[1] < hsp_start < trace.x[3]) for trace in figure.data)
            trancparency = 0.3 if show_containing and isContaining else 0.5

            # Inside the loop where you set the fillcolor
            """hsl_color = color[4:-1].replace('%', '').split(",")  # Remove % from HSL values
            rgb_color = colorsys.hls_to_rgb(float(hsl_color[0]) / 360, float(hsl_color[1]) / 100, float(hsl_color[2]) / 100)
            rgba_color = f'rgba({int(rgb_color[0] * 255)}, {int(rgb_color[1] * 255)}, {int(rgb_color[2] * 255)}, {trancparency})'"""

            if show_containing and isContaining:
                 y_shift = -1
            
            y_hight = 1 if show_containing and isContaining else 1
          
            if isContaining and not show_containing: 
                    pass
            else: 
                figure.add_trace(go.Scatter(
                        x=[hsp_start, hsp_end, hsp_end, hsp_start, hsp_start],
                        y=[y_shift, y_shift, y_shift + y_hight, y_shift + y_hight, y_shift],
                        mode='lines',
                        line=dict(color='black', width=1),
                        fill='toself',
                        fillcolor=f'rgba({color[4:-1]}, {trancparency})',
                        text=f"{query_name}<br>Start: {hsp_start}<br>End: {hsp_end}<br>queryStart: {query_start}<br>queryEnd: {query_end}",
                        hoverinfo='text',
                        name=query_name, 
                        showlegend=showLegend, 
                    ))
                encountered_names.add(query_name)
                
        # Check if the annotations DataFrame is not empty
        if not annotations_df.empty:

           
            # Define a color scale for annotations
            sorted_annotation_df = annotations_df[annotations_df['subject_name'] == subject].sort_values(by='s_start')
            annotation_color_scale = cl.interp(cl.scales['9']['seq']['Reds'], len(sorted_annotation_df))
            # Iterate over the annotations for the current subject
            for index, row in sorted_annotation_df.iterrows():
                annotation_start, annotation_end = row['s_start'], row['s_end']
                annotation_name = row['query_name']
                annotation_color = annotation_color_scale[index]

                showLegend = True
                if annotation_name in encountered_names:
                    showLegend = False

                # Add a trace for each annotation
                figure.add_trace(go.Scatter(
                    x=[annotation_start, annotation_end, annotation_end, annotation_start, annotation_start],
                    y=[-2, -2, -1, -1, -2],
                    mode='lines',
                    line=dict(color='red', width=1),
                    fill='toself',
                    fillcolor=annotation_color,
                    text=f"{annotation_name}<br>Start: {annotation_start}<br>End: {annotation_end}",
                    hoverinfo='text',
                    name=annotation_name,
                    showlegend=showLegend  
                ))
                encountered_names.add(annotation_name)

        figure.update_layout(yaxis=dict(visible=False))
        figures.append(figure)


    return figures


if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_hot_reload=False)
