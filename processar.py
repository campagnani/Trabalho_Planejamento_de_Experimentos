import os
os.system("clear")

import pandas as pd
import matplotlib.pyplot as plt
import ast
import numpy as np

#import resultados_experimentos_10000_220_20_2 as r
#import resultados_experimentos_1000_220_20 as r

#import resultados_experimentos_100000_220_20 as r
#import resultados_experimentos_1000_100_10 as r
#import resultados_experimentos_10000_100_10 as r
#import resultados_experimentos_10000_220_20 as r

import resultados_experimentos_1200_220_20_2 as r


reat = []
for i in range(len(r.vetor_keff)):
        reat.append((r.vetor_keff[i]-r.vetor_keff[0])*10**5)

# Create DataFrame
# Columns: F1, F2, F3, F4, F5
df = pd.DataFrame(r.matriz_planejamento_2k, columns=['F1', 'F2', 'F3', 'F4', 'F5'])
df['Reatividade'] = r.vetor_keff#reat

# Define Factor Labels
factor_labels = ["Enriquecimento (%)", "Temp Comb (°C)", "Dens Comb (g/cm³)", "Temp Mod (°C)", "Dens Mod (g/cm³)"]
factor_cols = ['F1', 'F2', 'F3', 'F4', 'F5']

# --- Plot 1: Main Effects ---
fig1, axes1 = plt.subplots(1, 5, figsize=(20, 5), sharey=True)
global_mean = df['Reatividade'].mean()

# Determine global min and max for y-axis scaling to make effects comparable
# Calculate all means to find range
all_means = []
for col in factor_cols:
    means = df.groupby(col)['Reatividade'].mean()
    all_means.extend(means.values)
y_min, y_max = min(all_means), max(all_means)
margin = (y_max - y_min) * 0.1

for i, col in enumerate(factor_cols):
    ax = axes1[i]
    means = df.groupby(col)['Reatividade'].mean()
    ax.plot(means.index, means.values, marker='o', linestyle='-', color='blue')
    ax.axhline(global_mean, color='gray', linestyle='--', linewidth=0.8, label='Média Global')
    ax.set_xlabel(factor_labels[i])
    ax.set_xticks([-1, 0, 1])
    ax.grid(True, linestyle=':', alpha=0.6)
    if i == 0:
        ax.set_ylabel('Média de Reatividade')

fig1.suptitle('Gráfico de Efeitos Principais no Reatividade', fontsize=16)
plt.ylim(y_min - margin, y_max + margin)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('main_effects_plot.png')
plt.show()

# --- Plot 2: Interaction Matrix ---
# We will create a 5x5 matrix.
# Diagonal: Factor names
# Off-diagonal (row i, col j): Interaction plot where x-axis is Factor j and lines are Factor i
# Only use data where factors are -1 or 1 (exclude center points for interactions)

df_corners = df[(df['F1'] != 0)].copy() # Assuming if F1 is 0, all are 0 (center points)
# Actually, verify if all 0s are center points. Yes, rows 0-3 are 0,0,0,0,0.
# Rows 4-35 are factorial points.

fig2, axes2 = plt.subplots(5, 5, figsize=(20, 20), sharey=True)

# Y-axis limits for interactions might differ, but usually better to share to compare magnitude.
# Let's find the range for interaction means.
int_means_list = []
for i in range(5):
    for j in range(5):
        if i != j:
            fi = factor_cols[i]
            fj = factor_cols[j]
            # Group by both factors
            means = df_corners.groupby([fj, fi])['Reatividade'].mean()
            int_means_list.extend(means.values)
y_min_int, y_max_int = min(int_means_list), max(int_means_list)
margin_int = (y_max_int - y_min_int) * 0.1


for i in range(5): # Row (Trace Factor / Legend Factor)
    for j in range(5): # Col (X-axis Factor)
        ax = axes2[i, j]
        
        if i == j:
            # Diagonal: Just text
            ax.text(0.5, 0.5, factor_labels[i], horizontalalignment='center', verticalalignment='center', fontsize=12, wrap=True)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            # Interaction Plot
            # X-axis: Factor j (fj)
            # Lines: Factor i (fi)
            fi = factor_cols[i] # This varies the line style/color
            fj = factor_cols[j] # This is the x-axis
            
            # Calculate means
            means = df_corners.groupby([fj, fi])['Reatividade'].mean().unstack()
            
            # Plot lines
            # means.index is fj levels (-1, 1)
            # means.columns is fi levels (-1, 1)
            
            ax.plot(means.index, means[-1], marker='o', linestyle='-', color='blue', label=f'{factor_labels[i]} -1')
            ax.plot(means.index, means[1], marker='s', linestyle='--', color='red', label=f'{factor_labels[i]} +1')
            
            ax.set_xticks([-1, 1])
            ax.grid(True, linestyle=':', alpha=0.6)
            
            # Only show y labels on first column
            if j == 0:
                ax.set_ylabel('Reatividade')
            
            # Fixed Y limits
            ax.set_ylim(y_min_int - margin_int, y_max_int + margin_int)

# Add a single legend for the whole figure if possible, or rely on colors.
# Blue solid: Row Factor Low (-1)
# Red dashed: Row Factor High (+1)
handles, labels = axes2[1, 0].get_legend_handles_labels()
fig2.legend(handles, ['Nível -1 (Linha)', 'Nível +1 (Linha)'], loc='upper right', title="Nível do Fator da Linha")

fig2.suptitle('Gráfico de Interação de Efeitos no Reatividade', fontsize=16)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('interaction_matrix_plot.png')
plt.show()


print("Plots generated.")