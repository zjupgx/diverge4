from typing import List, Tuple, Dict, Any
import math
import os
from pathlib import Path

import plotly.graph_objects as go
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade

def get_clade_name(clade_dict: Dict[Any, Tuple[float, float]], position: Tuple[float, float]) -> List[Any]:
    """
    Get the name of a clade based on its position.

    Args:
        clade_dict (Dict[Any, Tuple[float, float]]): Dictionary mapping clades to their positions.
        position (Tuple[float, float]): Position to search for.

    Returns:
        List[Any]: List of clade names at the given position.
    """
    return [k for k, v in clade_dict.items() if v == position]

def get_x_positions(tree: Tree) -> Dict[Clade, float]:
    """
    Create a mapping of each clade to its horizontal position.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        Dict[Clade, float]: Dictionary of {clade: x-coord}.
    """
    depths = tree.depths()
    # If there are no branch lengths, assume unit branch lengths
    if not max(depths.values()):
        depths = tree.depths(unit_branch_lengths=True)
    return depths

def get_y_positions(tree: Tree) -> Dict[Clade, float]:
    """
    Create a mapping of each clade to its vertical position.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        Dict[Clade, float]: Dictionary of {clade: y-coord}.
    """
    maxheight = tree.count_terminals()
    # Rows are defined by the tips
    heights = {
        tip: maxheight - i for i, tip in enumerate(reversed(tree.get_terminals()))
    }

    # Internal nodes: place at midpoint of children
    def calc_row(clade: Clade) -> None:
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

def calculate_font_size(tree: Tree) -> float:
    """
    Calculate font size based on the complexity of the tree.

    Args:
        tree (Tree): The phylogenetic tree.

    Returns:
        float: Calculated font size.
    """
    num_terminals = tree.count_terminals()
    base_size = 12  # Base font size
    # Adjust the font size based on the number of terminals (leaves)
    size = max(base_size - math.log(num_terminals, 2), 8)  # Minimum size 8
    return size

def draw_clade_lines(
    clade_dict: Dict[Clade, Tuple[float, float]],
    fig: go.Figure,
    orientation: str = "horizontal",
    y_here: float = 0,
    x_start: float = 0,
    x_here: float = 0,
    y_bot: float = 0,
    y_top: float = 0,
    color: str = "black",
    lw: float = 1.5,
) -> None:
    """
    Draw clade lines on the figure.

    Args:
        clade_dict (Dict[Clade, Tuple[float, float]]): Dictionary mapping clades to their positions.
        fig (go.Figure): The plotly figure to draw on.
        orientation (str, optional): Orientation of the line. Defaults to "horizontal".
        y_here (float, optional): Y-coordinate for horizontal lines. Defaults to 0.
        x_start (float, optional): Starting X-coordinate. Defaults to 0.
        x_here (float, optional): Ending X-coordinate. Defaults to 0.
        y_bot (float, optional): Bottom Y-coordinate for vertical lines. Defaults to 0.
        y_top (float, optional): Top Y-coordinate for vertical lines. Defaults to 0.
        color (str, optional): Line color. Defaults to "black".
        lw (float, optional): Line width. Defaults to 1.5.
    """
    if orientation == "horizontal":
        fig.add_trace(go.Scatter(x=[x_start, x_here], y=[y_here, y_here], marker_size=1, hoverinfo='skip', line=dict(color=color, width=lw)))
    elif orientation == "vertical":
        fig.add_trace(go.Scatter(x=[x_here, x_here], y=[y_bot, y_top], marker_size=1, hoverinfo='skip', line=dict(color=color, width=lw)))

def draw_clade(clade_dict: Dict[Clade, Tuple[float, float]], fig: go.Figure, clade: Clade, x_start: float, color: str = 'black', lw: float = 1.5) -> None:
    """
    Recursively draw a clade and its descendants.

    Args:
        clade_dict (Dict[Clade, Tuple[float, float]]): Dictionary mapping clades to their positions.
        fig (go.Figure): The plotly figure to draw on.
        clade (Clade): The clade to draw.
        x_start (float): Starting X-coordinate.
        color (str, optional): Line color. Defaults to 'black'.
        lw (float, optional): Line width. Defaults to 1.5.
    """
    x_here, y_here = clade_dict[clade]
    draw_clade_lines(
        clade_dict, fig,
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
            clade_dict, fig,
            orientation="vertical",
            x_here=x_here,
            y_bot=y_bot,
            y_top=y_top,
            color=color,
            lw=lw,
        )
        # Draw descendents
        for child in clade:
            draw_clade(clade_dict, fig, child, x_here, color, lw)

def draw_phylogenetic_tree(tree: Tree) -> Tuple[go.Figure, Dict[str, Tuple[float, float]], Tree]:
    """
    Draw a phylogenetic tree using plotly.

    Args:
        tree (Tree): The phylogenetic tree to draw.

    Returns:
        Tuple[go.Figure, Dict[str, Tuple[float, float]], Tree]: 
        The plotly figure, a dictionary mapping clade names to positions, and the tree.
    """
    x_posns = get_x_positions(tree)
    y_posns = get_y_positions(tree)
    x_coords = []
    y_coords = []
    node_names = []
    clade_dict = {}
    # set name for nonterminal clades
    for clade in tree.find_clades():
        x_coords.append(x_posns[clade])
        y_coords.append(y_posns[clade])
        node_names.append(clade.name)
        clade_dict[clade] = [x_posns[clade], y_posns[clade]]
    font_size = calculate_font_size(tree)
    scatter = go.Scatter(
        x=x_coords,
        y=y_coords,
        mode='markers+text',
        text=node_names,
        textfont=dict(size=font_size, family='Arial'),
        textposition='middle right',
        hoverinfo='none',
        marker=dict(
            size=5,
            color='grey',
        )
    )

    layout = go.Layout(
        xaxis=dict(
            title='', showgrid=False, zeroline=False, showticklabels=False
        ),
        yaxis=dict(
            title='', showgrid=False, zeroline=False, showticklabels=False
        ),
        showlegend=False,
        paper_bgcolor='white',
    )
    fig = go.Figure(data=[scatter], layout=layout)
    draw_clade(clade_dict, fig, tree.root, 0)
    clade_name_dict = {k.name: v for k, v in clade_dict.items()}
    return fig, clade_name_dict, tree