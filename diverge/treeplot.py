from Bio import Phylo
import plotly.graph_objects as go
import math
import streamlit as st
import os
from pathlib import Path
from streamlit_plotly_events import plotly_events
import re
import copy
state = st.session_state


# Some constants for the drawing calculations

def get_clade_name(clade_dict,position):
      return [k for k,v in clade_dict.items() if v == position]

def get_x_positions(tree):
        """Create a mapping of each clade to its horizontal position.

        Dict of {clade: x-coord}
        """
        depths = tree.depths()
        # If there are no branch lengths, assume unit branch lengths
        if not max(depths.values()):
            depths = tree.depths(unit_branch_lengths=True)
        return depths


def get_y_positions(tree):
        """Create a mapping of each clade to its vertical position.

        Dict of {clade: y-coord}.
        Coordinates are negative, and integers for tips.
        """
        maxheight = tree.count_terminals()
        # Rows are defined by the tips
        heights = {
            tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
        }

        # Internal nodes: place at midpoint of children
        def calc_row(clade):
            for subclade in clade:
                if subclade not in heights:
                    calc_row(subclade)
            # Closure over heights
            heights[clade] = (
                heights[clade.clades[0]] + heights[clade.clades[-1]]
            ) / 2.0

        if tree.root.clades:
            calc_row(tree.root)
        return heights
    
def calculate_font_size(tree):
    """Calculate font size based on the complexity of the tree."""
    num_terminals = tree.count_terminals()
    base_size = 12  # Base font size
    # Adjust the font size based on the number of terminals (leaves)
    size = max(base_size - math.log(num_terminals, 2), 8)  # Minimum size 6
    return size




def draw_clade_lines(
        clade_dict,fig,
        orientation="horizontal",
        y_here=0,
        x_start=0,
        x_here=0,
        y_bot=0,
        y_top=0,
        color="black",
        lw= 1.5,
    ):
    if orientation == "horizontal":
        fig.add_trace(go.Scatter(x=[x_start, x_here], y=[y_here, y_here],marker_size=1,hoverinfo='skip',line=dict(color=color, width=lw)))
        # fig.add_trace(type="line", x0=x_start, y0=y_here, x1=x_here, y1=y_here, line=dict(color=color, width=lw))
    elif orientation == "vertical":
        fig.add_trace(go.Scatter(x=[x_here, x_here], y=[y_bot, y_top],marker_size=1,hoverinfo='skip',line=dict(color=color, width=lw)))
        # fig.add_shape(type="line", x0=x_here, y0=y_bot, x1=x_here, y1=y_top, line=dict(color=color, width=lw))

def draw_clade(clade_dict,fig,clade, x_start ,color = 'black', lw = 1.5):
    x_here = clade_dict[clade][0]
    y_here = clade_dict[clade][1]
    draw_clade_lines(
            clade_dict,fig,
            orientation="horizontal",
            y_here=y_here,
            x_start=x_start,
            x_here=x_here,
            color=color,
            lw=lw,
        )
    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = clade_dict[clade.clades[0]][1]
        y_bot = clade_dict[clade.clades[-1]][1]
        # Only apply widths to horizontal lines, like Archaeopteryx
        draw_clade_lines(
            clade_dict,fig,
            orientation="vertical",
            x_here=x_here,
            y_bot=y_bot,
            y_top=y_top,
            color=color,
            lw=lw,
        )
        # Draw descendents
        for child in clade:
            draw_clade(clade_dict,fig,child, x_here, color, lw)


def draw_phylogenetic_tree(tree):    
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    x_coords = []
    y_coords = []
    node_names = []
    clade_dict = {}
# set name for nonterminal clades
    index = 0
    for clade in tree.find_clades():
        x_coords.append(x_posns[clade])
        y_coords.append(y_posns[clade])
        node_names.append(clade.name)
        clade_dict[clade] = [x_posns[clade],y_posns[clade]]
    font_size = calculate_font_size(tree)
    scatter = go.Scatter(
        x=x_coords,
        y=y_coords,
        mode='markers+text',
        text=node_names,
        textfont=dict(size=font_size,family='Arial'),
        textposition='middle right',
        hoverinfo='none',
        marker=dict(
            size=5,
            color='grey',
        )
    )

    layout = go.Layout(
        # title='Evolutionary Tree',
        xaxis=dict(
            title='',showgrid=False, zeroline=False, showticklabels=False
        ),
        yaxis=dict(
            title='',showgrid=False, zeroline=False, showticklabels=False
        ),
        showlegend=False,
        paper_bgcolor='white',
    )
    # 绘制进化树
    fig = go.Figure(data=[scatter], layout=layout)
    draw_clade(clade_dict,fig,tree.root, 0)
    clade_name_dict = {}
    for k,v in clade_dict.items():
        clade_name_dict[k.name] = v
    return fig,clade_name_dict,tree
    # selected_points = plotly_events(fig, click_event=True, hover_event=False,override_height=600,override_width="130%")
    # selected_clade = get_clade_name(clade_dict,[selected_points[0]['x'],selected_points[0]['y']])[0]
    # if selected_clade is not None:
    #     return tree,selected_clade
    # else:
    #     return tree,None

  
  

