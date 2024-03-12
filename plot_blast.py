
import plotly.express as px
import pandas as pd
import plotly.figure_factory as ff
# Input BLAST results file in outfmt 6 format
'''
blast_results_file = '/home/michele/git/4_in_house_seq/1_vacs_insert_sites/2_integration_site_fromChr13ToInsert/7_consensus_sequences/1_cons_seq_chr13_insert_flanks/3_blast_unknown_part_cr/blast_cons_unknown_cr_table.txt'

def parse_blast_results(blast_file):
    database_sequences = []
    alignment_scores = []
    alignment_positions = []

    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 12:
                continue
            database_seq = fields[1]
            alignment_score = float(fields[11])
            alignment_start = int(fields[8])
            alignment_end = int(fields[9])
            if alignment_start > alignment_end:
                alignment_start, alignment_end = alignment_end, alignment_start  # Swap positions for negative strand alignments
            database_sequences.append(database_seq)
            alignment_scores.append(alignment_score)
            alignment_positions.append((alignment_start, alignment_end))
    
    return database_sequences, alignment_scores, alignment_positions

# Parse BLAST results
database_sequences, alignment_scores, alignment_positions = parse_blast_results(blast_results_file)
'''
import pandas as pd
import plotly.figure_factory as ff

# Sample data (replace with your parsed BLAST results)
database_sequences = ["Seq1", "Seq2", "Seq2", "Seq3"]  # Unsorted sequence names
alignment_scores = [90, 80, 95, 70]  # Replace with your actual scores
alignment_positions = [(10, 50), (30, 80), (100, 120), (5, 60)]  # Replace with actual positions

# Pair sequences with their scores and positions
data = list(zip(database_sequences, alignment_scores, alignment_positions))

# Sort the data based on sequence names and alignment start position
data.sort(key=lambda x: (x[0], x[2][0]))  

# Create a DataFrame for Gantt chart
df = pd.DataFrame([
    dict(Task=seq, Start=start, Finish=end, Score=score)
    for seq, score, (start, end) in data
])

# Create the Gantt chart using ff.create_gantt
fig = ff.create_gantt(df, index_col='Task', colors=['blue'], title='Interactive BLAST Results',
                      show_colorbar=True, bar_width=1)  # Set bar_width to an integer value

fig.update_layout(xaxis_type='linear', autosize=False, width=800, height=400)
fig.update_traces(text=df['Task'] + ' (Score: ' + df['Score'].astype(str) + ')')

fig.show()

