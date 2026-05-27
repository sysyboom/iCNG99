import plotly.graph_objects as go

before_rev = 2514
before_irrev = 871
total = before_rev + before_irrev

after_rev = 583
after_irrev = 2802

rev_to_rev = after_rev                  # 583
rev_to_irrev = before_rev - after_rev   # 1931
irrev_to_irrev = before_irrev           # 871

# sanity check
assert rev_to_rev + rev_to_irrev == before_rev
assert rev_to_irrev + irrev_to_irrev == after_irrev
assert after_rev + after_irrev == total

labels = [
    f"Before correction<br>Reversible<br>{before_rev:,}<br>({before_rev / total * 100:.1f}%)",
    f"Before correction<br>Irreversible<br>{before_irrev:,}<br>({before_irrev / total * 100:.1f}%)",
    f"After correction<br>Reversible<br>{after_rev:,}<br>({after_rev / total * 100:.1f}%)",
    f"After correction<br>Irreversible<br>{after_irrev:,}<br>({after_irrev / total * 100:.1f}%)"
]

source = [0, 0, 1]
target = [2, 3, 3]
value = [rev_to_rev, rev_to_irrev, irrev_to_irrev]

link_labels = [
    f"{rev_to_rev:,} ({rev_to_rev / total * 100:.1f}%)",
    f"{rev_to_irrev:,} ({rev_to_irrev / total * 100:.1f}%)",
    f"{irrev_to_irrev:,} ({irrev_to_irrev / total * 100:.1f}%)"
]

orange_node = "rgba(255, 130, 20, 0.95)"
green_node = "rgba(120, 190, 80, 0.95)"

orange_link_light = "rgba(255, 130, 20, 0.35)"
orange_link_dark = "rgba(255, 130, 20, 0.55)"
green_link = "rgba(120, 190, 80, 0.35)"

fig = go.Figure(data=[go.Sankey(
    arrangement="fixed",

    node=dict(
        pad=35,
        thickness=40,
        line=dict(color="white", width=1),
        label=labels,
        color=[
            orange_node,
            green_node,
            orange_node,
            green_node
        ],
        x=[0.08, 0.08, 0.92, 0.92],
        y=[0.20, 0.72, 0.20, 0.72],
    ),

    link=dict(
        source=source,
        target=target,
        value=value,
        color=[
            orange_link_light,
            orange_link_dark,
            green_link
        ],
        label=link_labels
    )
)])

fig.update_layout(
    title=dict(
        text="<b>Correction of reaction reversibility</b>",
        x=0.50,
        y=0.95,
        xanchor="center",
        yanchor="top",
        font=dict(size=26)
    ),
    font=dict(
        size=18,
        family="Arial"
    ),
    width=1000,
    height=650,
    margin=dict(l=50, r=50, t=90, b=70)
)

fig.add_annotation(
    x=0.5,
    y=-0.08,
    xref="paper",
    yref="paper",
    text=(
        f"<b>Reversible → Irreversible:</b> "
        f"{rev_to_irrev:,} ({rev_to_irrev / total * 100:.1f}%)"
    ),
    showarrow=False,
    font=dict(size=20, color="rgb(255,130,20)"),
    bordercolor="rgba(80,140,200,0.7)",
    borderwidth=2,
    borderpad=10,
    bgcolor="rgba(255,255,255,0.8)"
)

fig.write_html("reaction_reversibility_sankey.html")

# Requires: pip install kaleido
fig.write_image("reaction_reversibility_sankey.pdf")
fig.write_image("reaction_reversibility_sankey.png", scale=3)

fig.show()