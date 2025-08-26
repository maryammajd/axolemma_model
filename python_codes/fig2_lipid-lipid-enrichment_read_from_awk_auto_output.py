#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 16:56:47 2024

@author: maryamma
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

ratio_bulk = {'axl': {"pc-up": 0.3679, "pe-up": 0.1925, "sm-up": 0.0585, "gl-up": 0.3811,
              "pc-down": 0.2574, "pe-down": 0.4813, "pi-down": 0.1336, "ps-down": 0.1022, "sm-down": 0.0255},
              'node': {"pc-up": 0.38, "pe-up": 0.19, "sm-up": 0.07, "gl-up": 0.36,
                       "pc-down": 0.27, "pe-down": 0.47, "pi-down": 0.13, "ps-down": 0.10, "sm-down": 0.03},
              'large': {"pc-up": 0.37, "pe-up": 0.19, "sm-up": 0.06, "gl-up": 0.38,
                        "pc-down": 0.26, "pe-down": 0.48, "pi-down": 0.14, "ps-down": 0.10, "sm-down": 0.02}}

absolute_neighboring, relative_neighboring = {}, {}
colnames = ['time', 'number']

source_dir = "/upload_files"

directories = [source_dir + "/node-of-Ranvier/neighbor/neighbor-all/",
               source_dir + "/one-protein-axolemma/neighbor/",
               source_dir + "/axolemma/neighbor/"]
systems = ['node-of-Ranvier', 'one-protein-axolemma', 'axolemma']


comment_chars = ['@', '#']
for dir, system_chosen in zip(directories, systems):
    first_folder_list = [d for d in sorted(
        os.listdir(dir)) if os.path.isdir(os.path.join(dir, d))]
    absolute_neighboring[system_chosen], relative_neighboring[system_chosen] = {
    }, {}
    for first_folder in first_folder_list:
        second_folder_list = [d for d in sorted(os.listdir(dir + '/' + first_folder)) if os.path.isdir(
            os.path.join(dir + '/' + first_folder, d)) and not d.startswith('.')]
        for second_folder in second_folder_list:
            if second_folder.startswith('all'):
                new_dir = dir + first_folder + '/' + second_folder
                filename_list = sorted(os.listdir(new_dir))
                data_frame_all = pd.DataFrame()
                os.chdir(new_dir)
                df = pd.read_csv(new_dir + '/' + second_folder + '.xvg')
                time_ref = df['time'].astype(float).to_numpy()/1000000
                neighboring_value_ref = df['number'].astype(float).to_numpy()
            elif not second_folder.startswith('chol'):
                new_dir = dir + first_folder + '/' + second_folder
                filename_list = sorted(os.listdir(new_dir))
                data_frame_all = pd.DataFrame()
                os.chdir(new_dir)
                df = pd.read_csv(new_dir + '/' + second_folder + '.xvg')
                time = df['time'].astype(float).to_numpy()/1000000
                neighboring_value = df['number'].astype(float).to_numpy()
                relative = []
                for rel_num in range(0, len(neighboring_value)):
                    relative.append(
                        np.round(neighboring_value[rel_num]/neighboring_value_ref[rel_num], 4))
                absolute_neighboring[system_chosen][second_folder] = neighboring_value
                relative_neighboring[system_chosen][second_folder] = relative

mean_neigh, std_neigh, enrichment_factor = {}, {}, {}
for system in relative_neighboring.keys():
    for lipid_lipid, value_of_neigh in relative_neighboring[system].items():
        if lipid_lipid not in mean_neigh.keys():
            mean_neigh[lipid_lipid], enrichment_factor[lipid_lipid], std_neigh[lipid_lipid] = {
            }, {}, {}
        mean_neigh[lipid_lipid][system] = round(
            np.mean(value_of_neigh[-1000:]), 4)
        std_neigh[lipid_lipid][system] = round(
            np.std(value_of_neigh[-1000:]), 4)
        enrichment_factor[lipid_lipid][system] = round(
            mean_neigh[lipid_lipid][system]/ratio_bulk[system][lipid_lipid[:3]+lipid_lipid[6:]], 2)
for i in enrichment_factor.keys():
    print(i, enrichment_factor[i])

fig, axes = plt.subplots(
    nrows=2, ncols=1, figsize=(1.2*10, 10), dpi=600)

color = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3',
         '#999999', '#e41a1c', '#dede00']


for plot_num, leaflet in enumerate(['up', 'down']):
    if leaflet == 'up':
        leaf = 'Extracellular'
    if leaflet == 'down':
        leaf = 'Intracellular'

    x_values = list(key[:-len(leaflet)-1].upper()
                    for key in enrichment_factor.keys() if (key.endswith(leaflet)))
    y_values_node = [enrichment_factor[key.lower() + '-' + leaflet]['node']
                     for key in x_values]
    y_err_node = [std_neigh[key.lower() + '-' + leaflet]['node']
                  for key in x_values]
    y_values_large = [enrichment_factor[key.lower() + '-' + leaflet]['large']
                      for key in x_values]
    y_err_large = [std_neigh[key.lower() + '-' + leaflet]['large']
                   for key in x_values]
    y_values_axl = [enrichment_factor[key.lower() + '-' + leaflet]['axl']
                    for key in x_values]
    y_err_axl = [std_neigh[key.lower() + '-' + leaflet]['axl']
                 for key in x_values]

    bar_width = 0.25

    index = np.arange(len(x_values))


    axes[plot_num].bar(index, y_values_axl,
                       bar_width, label='axolemma', color=color[0], align='edge', zorder=5)
    axes[plot_num].errorbar(index + bar_width/2, y_values_axl,
                            yerr=y_err_axl, fmt=".", markersize=2, capsize=3, zorder=6, color=color[0])

    axes[plot_num].bar(index + bar_width + 0.02, y_values_large,
                       bar_width, label='one-protein-axolemma', color=color[1], align='edge', zorder=5)
    axes[plot_num].errorbar(index + 3/2*bar_width + 0.02, y_values_large,
                            yerr=y_err_large, fmt=".", markersize=2, capsize=3, zorder=6, color=color[1])

    axes[plot_num].bar(index + 2 * bar_width + 0.04, y_values_node, bar_width,
                       label='node-of-Ranvier', color=color[2], align='edge', zorder=5)
    axes[plot_num].errorbar(index + 5/2 * bar_width + 0.02, y_values_node,
                            yerr=y_err_node, fmt=".", markersize=2, capsize=3, zorder=6, color=color[2])

    axes[plot_num].grid(color="lightgrey",
                        linewidth="0.6", zorder=0)

    axes[plot_num].set_xticks(index + 3/2 * bar_width,
                              x_values, rotation=45, ha='right')
    axes[plot_num].set_ylabel('Enrichment factor', fontsize=15)
    axes[plot_num].text(-0.05,
                        (1.64 - 0.01 * plot_num),  leaf + ' leaflet', fontsize=12)
    axes[plot_num].set_ylim(0, 1.751)

    axes[plot_num].margins(x=0.01)  

plt.legend(fontsize=12, bbox_to_anchor=(0.81, -0.12), ncol=3)
plt.show()
