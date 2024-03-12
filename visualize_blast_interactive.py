import dash
from dash import dcc, html
import plotly.graph_objs as go
import pandas as pd
import colorlover as cl
import argparse

# Function to parse command-line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Interactive BLAST Results Visualization")
    parser.add_argument("blast_results_file", help="Path to the BLAST results file")
    parser.add_argument("--annotat", help="Path to the BLAST annotations file")
    parser.add_argument("--sc", dest="show_containing", action="store_true", help="Include containing hits in the visualization")

    # Parse the command-line arguments
    args = parser.parse_args()
    
    # Return parsed arguments
    return args.blast_results_file, args.annotat, args.show_containing if hasattr(args, 'show_containing') and args.show_containing else False

# Function to read BLAST results from a file
def read_blast_file(file_path):
    if not file_path:
        return pd.DataFrame()
    df = pd.read_csv(file_path, sep='\t', header=None)
    df.columns = ["query_name", "subject_name", "percent_identity", "alignment_length", "mismatches", "gap_opens",
                  "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score"]
    return df


# Function to create scatter trace 
def create_trace(x, y, color, lineColor, text, name, show_legend):
    return go.Scatter(
        x=x,
        y=y,
        mode='lines',
        line=dict(color=lineColor, width=1),
        fill='toself',
        fillcolor=color,
        text=text,
        hoverinfo='text',
        name=name,
        showlegend=show_legend
    )

# Function to update BLAST visualization based on user interaction
def update_blast_visualization(*args):
    figures = []

    for subject in subjects:
        # Create a new figure for each subject
        figure = go.Figure()
        # Keep track of encountered annotation names
        encountered_names = set()




        # Iterate over BLAST results for the current subject
        for index, row in blast_df[blast_df['subject_name'] == subject].iterrows():

            hsp_start, hsp_end = row['s_start'], row['s_end']
            y_shift = 0
            is_containing = any((trace.x[0] < hsp_end < trace.x[2]) or (trace.x[1] < hsp_start < trace.x[3]) for trace in figure.data)
            transparency = 0.3 if show_containing and is_containing else 0.5
            y_shift = -1.2 if show_containing and is_containing else 0

            if is_containing and not show_containing:
                pass
            else:
                # Add scatter trace for each BLAST result
                query_name = row['query_name']
                query_start, query_end = row['q_start'], row['q_end']
                y_height = 1 if show_containing and is_containing else 1
                show_legend = query_name not in encountered_names
                color = query_colors[query_name]
                colorRGB = f'rgba({color[4:-1]},{transparency})'
                figure.add_trace(create_trace(
                    x=[hsp_start, hsp_end, hsp_end, hsp_start, hsp_start],
                    y=[y_shift, y_shift, y_shift + y_height, y_shift + y_height, y_shift],
                    color=colorRGB,
                    lineColor='black',
                    text=f"{query_name}<br>Start: {hsp_start}<br>End: {hsp_end}<br>queryStart: {query_start}<br>queryEnd: {query_end}",
                    name=query_name,
                    show_legend=show_legend
                ))
                encountered_names.add(query_name)

        # Check if the annotations DataFrame is not empty
        if not annotations_df.empty:
            # Sort annotations by start position
            sorted_annotation_df = annotations_df[annotations_df['subject_name'] == subject].sort_values(by='s_start')
            # Define a color scale for annotations
            annotation_color_scale = cl.interp(cl.scales['9']['seq']['Reds'], len(sorted_annotation_df))

            # Iterate over annotations for the current subject
            for index, row in sorted_annotation_df.iterrows():
                annotation_start, annotation_end = row['s_start'], row['s_end']
                annotation_name = row['query_name']
                annotation_color = annotation_color_scale[index]

                show_legend = annotation_name not in encountered_names

                # Add annotation trace for each annotation
                figure.add_trace(create_trace(
                    x=[annotation_start, annotation_end, annotation_end, annotation_start, annotation_start],
                    y=[-2, -2, -1, -1, -2],
                    color=annotation_color,
                    lineColor='red',
                    text=f"{annotation_name}<br>Start: {annotation_start}<br>End: {annotation_end}",
                    name=annotation_name,
                    show_legend=show_legend  
                ))
                encountered_names.add(annotation_name)

        # Update layout and append the figure to the list of figures
        figure.update_layout(yaxis=dict(visible=False))
        figures.append(figure)

    return figures


# Call the function to get command-line arguments
blast_results_file, blast_annotations_file, show_containing = parse_arguments()

# Read BLAST results and annotations
blast_df = read_blast_file(blast_results_file)
annotations_df = read_blast_file(blast_annotations_file)

# Get unique queries and subjects
queries = blast_df['query_name'].unique()
subjects = blast_df['subject_name'].unique()

# Define the number of colors and color scale
num_colors = len(queries)
color_scale = cl.scales['9']['qual']['Set1'][:num_colors]

# Create a dictionary to map queries to colors
query_colors = {query: color for query, color in zip(queries, color_scale)}

# Initialize the Dash app
app = dash.Dash(__name__)


app.layout = html.Div([
    dcc.Tabs([
        dcc.Tab(label=subject, children=[
            dcc.Graph(
                id=f'blast-visualization-{subject}',
                config={'displayModeBar': True, 'scrollZoom': True, 'doubleClick': 'reset'},
                figure={'data': [], 'layout': {'xaxis.autorange': True, 'zoom': 'reset'}},
            )
        ]) for subject in subjects
    ]),
])

app.callback(
    [dash.dependencies.Output(f'blast-visualization-{subject}', 'figure') for subject in subjects],
    [dash.dependencies.Input(f'blast-visualization-{subject}', 'relayoutData') for subject in subjects]
)(update_blast_visualization)

if __name__ == '__main__':
    app.run_server(debug=True, dev_tools_hot_reload=False)
