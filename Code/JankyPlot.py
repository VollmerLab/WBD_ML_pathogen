from Python_Data_Wrangling import *

import plotly.graph_objs as go

# Define the source and target nodes for the Sankey plot
sources = fullSampledf[features].values.flatten()
targets = fullSampledf['health'].values.flatten()

print("Features:", features)
print("Sources:", sources)
# Create a dictionary to count the flow between source and target nodes
flow_counts = {}
for source, target in zip(sources, targets):
    flow_counts[(source, target)] = flow_counts.get((source, target), 0) + 1

source_indices = []
target_indices = []
values = []

# Convert the flow_counts dictionary into the required format for Plotly Sankey plot
for (source, target), value in flow_counts.items():
    source_idx = features.get_loc(source)  # Get the index of the source feature
    target_idx = len(features) + target.index(target)
    source_indices.append(source_idx)
    target_indices.append(target_idx)
    values.append(value)

# Create the Plotly Sankey plot
fig = go.Figure(go.Sankey(
    arrangement="snap",
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=features.tolist() + y[target].unique().tolist(),
    ),
    link=dict(
        source=source_indices,
        target=target_indices,
        value=values,
    )
))

# Show the plot
fig.show()