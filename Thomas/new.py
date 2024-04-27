import plotly.graph_objects as go

# Create a figure
fig = go.Figure()
x = [1,2,3,4,5,6]
y = [1,3,5,7,9,11]
# Add traces
fig.add_trace(go.Scatter(x=x, y=y))

# Enable zooming and panning
fig.update_layout(
    xaxis=dict(
        rangeslider=dict(
            visible=True
        ),
        type="linear"
    )
)

# Show the plot
fig.show()