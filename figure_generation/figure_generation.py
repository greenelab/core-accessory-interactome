# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.9.1+dev
#   kernelspec:
#     display_name: Python [conda env:core_acc] *
#     language: python
#     name: conda-env-core_acc-py
# ---

# # Figure generation

from IPython.display import Image, display, SVG
import svgutils.transform as sg
import numpy as np
from lxml import etree
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# Directory of output figures
# local_directory = "/home/alexandra/Documents/Data/Generic_expression_patterns/"
output_directory = "output/"
os.makedirs(output_directory, exist_ok=True)


# ## Function to plot

def make_figure_panel(filename, scale_x_input, scale_y_input, x_loc, y_loc):
    panel = sg.fromfile(filename)

    panel_size = (
        np.round(float(panel.root.attrib["width"][:-2]) * 1.33, 0),
        np.round(float(panel.root.attrib["height"][:-2]) * 1.33, 0),
    )

    scale_x = scale_x_input
    scale_y = scale_y_input

    print(f"original: {panel_size}")
    print(f"scaled:{(panel_size[0]*scale_x,panel_size[1]*scale_y)}")

    panel = panel.getroot()
    panel.scale_xy(x=scale_x, y=scale_y)
    panel.moveto(x_loc, y_loc)

    return panel


# ## Figure 1

# Create panels for figure 1
panel_1a = make_figure_panel(
    "../1_processing/compendia_media.svg",
    scale_x_input=0.5,
    scale_y_input=0.5,
    x_loc=30,
    y_loc=10,
)
panel_1b = make_figure_panel(
    "../1_processing/compendia_gene_function.svg",
    scale_x_input=0.5,
    scale_y_input=0.5,
    x_loc=600,
    y_loc=10,
)
panel_1c = make_figure_panel(
    "../1_processing/compendia_kegg.svg",
    scale_x_input=0.5,
    scale_y_input=0.5,
    x_loc=30,
    y_loc=300,
)

panel_1a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_1b_label = sg.TextElement(600, 20, "B", size=18, weight="bold", font="Verdana")
panel_1c_label = sg.TextElement(10, 300, "C", size=18, weight="bold", font="Verdana")

figure_1 = sg.SVGFigure("1200", "600")
figure_1.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_1a,
        panel_1b,
        panel_1c,
        panel_1a_label,
        panel_1b_label,
        panel_1c_label,
    ]
)
display(SVG(figure_1.to_str()))

# save generated SVG files
figure_1.save("output/figure_1.svg")

# ## Figure 2

# Create panels for figure 1
panel_2a = make_figure_panel(
    "fig2A_core_stability_workflow.svg",
    scale_x_input=3,
    scale_y_input=3,
    x_loc=30,
    y_loc=10,
)
panel_2b = make_figure_panel(
    "../3_core_core_analysis/pao1_similarity_scores_dist_spell.svg",
    scale_x_input=0.9,
    scale_y_input=0.9,
    x_loc=700,
    y_loc=30,
)
panel_2c = make_figure_panel(
    "fig2C_homolog_dist.svg",
    scale_x_input=3,
    scale_y_input=3,
    x_loc=30,
    y_loc=400,
)
panel_2d = make_figure_panel(
    "../3_core_core_analysis/stability_percent_match_homolog.svg",
    scale_x_input=0.7,
    scale_y_input=0.7,
    x_loc=700,
    y_loc=400,
)

panel_2a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_2b_label = sg.TextElement(700, 20, "B", size=18, weight="bold", font="Verdana")
panel_2c_label = sg.TextElement(10, 400, "B", size=18, weight="bold", font="Verdana")
panel_2d_label = sg.TextElement(700, 400, "C", size=18, weight="bold", font="Verdana")

figure_2 = sg.SVGFigure("1200", "800")
figure_2.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_2a,
        panel_2b,
        panel_2c,
        panel_2d,
        panel_2a_label,
        panel_2b_label,
        panel_2c_label,
        panel_2d_label,
    ]
)
display(SVG(figure_2.to_str()))

# save generated SVG files
figure_2.save("output/figure_2.svg")

# ## Figure 3

# Create panels for figure 2
panel_3a = make_figure_panel(
    "fig3A_core_acc_calc.svg", scale_x_input=4, scale_y_input=4, x_loc=30, y_loc=10
)
panel_3b_right = make_figure_panel(
    "../5_core_acc_analysis/PAO1_stablility_expression_relationships_operon_corrected_spell.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=30,
    y_loc=500,
)
panel_3b_left = make_figure_panel(
    "../5_core_acc_analysis/PA14_stability_expression_relationships_operon_corrected_spell.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=450,
    y_loc=500,
)
panel_3c = make_figure_panel(
    "../5_core_acc_analysis/core_genes_correlated_with_exo.svg",
    scale_x_input=1,
    scale_y_input=1,
    x_loc=30,
    y_loc=900,
)

panel_3a_label = sg.TextElement(10, 20, "A", size=18, weight="bold", font="Verdana")
panel_3b_label = sg.TextElement(10, 500, "B", size=18, weight="bold", font="Verdana")
panel_3c_label = sg.TextElement(10, 900, "C", size=18, weight="bold", font="Verdana")

figure_3 = sg.SVGFigure("1000", "1200")
figure_3.append(
    [
        etree.Element("rect", {"width": "100%", "height": "100%", "fill": "white"}),
        panel_3a,
        panel_3b_right,
        panel_3b_left,
        panel_3c,
        panel_3a_label,
        panel_3b_label,
        panel_3c_label,
    ]
)
display(SVG(figure_3.to_str()))

# save generated SVG files
figure_3.save("output/figure_3.svg")

# ## Output png version

# Exports low resolution png just for easy viewing,
# But will need to use inkscape to export high resolution images
# !inkscape --export-png=output/figure_1.png output/figure_1.svg
# !inkscape --export-png=output/figure_2.png output/figure_2.svg
# !inkscape --export-png=output/figure_3.png output/figure_3.svg
