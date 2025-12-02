import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy.stats as stats
import seaborn as sns
import numpy as np
from itertools import combinations


#import resultados_experimentos_1000_100_10 as r
#import resultados_experimentos_10000_100_10 as r
#import resultados_experimentos_1000_220_20 as r
#import resultados_experimentos_10000_220_20_2 as r
#import resultados_experimentos_1200_220_20 as r


#import resultados_experimentos_10000_220_20 as r
#import resultados_experimentos_100000_220_20 as r

#import resultados_experimentos_1200_220_20_2 as r
import results_sealer_20251202_011203.resultados_experimentos__12_793_18p5_793_10___0p1_500_0p5_500_1___1200_120_20_ as r
import results_sealer_20251202_014127.resultados_experimentos__12_1273_18p5_1273_10___0p2_1000_0p5_1000_1p5___1200_140_40_ as r
import results_sealer_20251202_020816.resultados_experimentos__12_1273_18p5_1273_10___0p2_1000_0p5_1000_1p5___50000_50_10_ as r
#import results_sealer_20251202_054041.resultados_experimentos__12_1273_18P5_1273_10___0P1_500_0P25_500_0P75___50000_50_10_ as r


reat = []
media = (r.vetor_keff[0]+r.vetor_keff[1]+r.vetor_keff[2]+r.vetor_keff[3])/4
for i in range(len(r.vetor_keff)):
        reat.append((r.vetor_keff[i]-media)*10**5)



vetor_keff = reat#r.vetor_keff 
matriz = r.matriz_planejamento_2k


# Nomes dos fatores conforme sua descrição
factor_names = ['Enriq', 'Temp1', 'Dens1', 'Temp2', 'Dens2']

# Criar DataFrame base
df = pd.DataFrame(matriz, columns=factor_names)

# 2. Adicionar Termos de Interação (2ª Ordem)
# Em DOE 2^k, geralmente modelamos Efeitos Principais + Interações de 2 fatores
for f1, f2 in combinations(factor_names, 2):
    interaction_name = f"{f1}*{f2}"
    df[interaction_name] = df[f1] * df[f2]

# Adicionar constante (Intercepto/Média Global) para o statsmodels
X = sm.add_constant(df)
y = vetor_keff

# 3. Ajustar o Modelo de Regressão (OLS)
model = sm.OLS(y, X).fit()

# 4. Imprimir a Equação do Modelo
print("="*80)
print("EQUAÇÃO DO MODELO DE REGRESSÃO (Keff)")
print("="*80)

equation_terms = []
# Pegar o intercepto primeiro
intercept = model.params['const']
equation_terms.append(f"{intercept:.5f}")

# Pegar os outros termos
for name, coef in model.params.items():
    if name == 'const': continue
    sign = "+" if coef >= 0 else "-"
    equation_terms.append(f"{sign} {abs(coef):.5f}*({name})")

# Quebrar a linha para visualização no terminal
equation_str = "Keff = " + equation_terms[0]
line_len = len(equation_str)
for term in equation_terms[1:]:
    if line_len + len(term) > 80:
        equation_str += "\n       " + term
        line_len = 7 + len(term)
    else:
        equation_str += " " + term
        line_len += 1 + len(term)

print(equation_str)
print("\nNota: O modelo considera efeitos principais e interações de 2ª ordem.")
print(f"R-squared: {model.rsquared:.4f}")
print("="*80)

# 5. Plotar Gráfico de Resíduos
residuals = model.resid
fitted_values = model.fittedvalues

plt.figure(figsize=(14, 6))

# Gráfico 1: Resíduos vs Valores Ajustados (Verificação de Homocedasticidade)
plt.subplot(1, 2, 1)
plt.scatter(fitted_values, residuals, color='blue', alpha=0.6, edgecolor='k')
plt.axhline(0, color='red', linestyle='--', linewidth=1)
plt.title('Resíduos vs. Valores Preditos')
plt.xlabel('Valores Preditos (Fitted Values)')
plt.ylabel('Resíduos')
plt.grid(True, linestyle=':', alpha=0.6)

# Gráfico 2: Histograma dos Resíduos (Verificação de Normalidade)
plt.subplot(1, 2, 2)
sns.histplot(residuals, kde=True, color='green', bins=10)
plt.title('Histograma dos Resíduos')
plt.xlabel('Resíduos')
plt.ylabel('Frequência')
plt.grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()
plt.savefig('grafico_residuos_modelo.png')
plt.show()

print("Gráfico de resíduos salvo como 'grafico_residuos_modelo.png'")



# 2. Criar modelo com Interações para calcular os efeitos
for f1, f2 in combinations(factor_names, 2):
    interaction_name = f"{f1}*{f2}"
    df[interaction_name] = df[f1] * df[f2]

X = sm.add_constant(df)
y = vetor_keff

# Ajustar modelo OLS
model = sm.OLS(y, X).fit()

# 3. Preparar dados para o "Normal Plot of Effects"
# Na codificação -1 a +1, o Efeito é 2x o Coeficiente
effects = model.params.drop('const') * 2
pvalues = model.pvalues.drop('const')

# Criar um DataFrame para facilitar a plotagem
df_effects = pd.DataFrame({
    'Effect': effects,
    'PValue': pvalues
})

# Ordenar pelos efeitos (do menor para o maior) para o plot de probabilidade
df_effects = df_effects.sort_values(by='Effect')

# Calcular os quantis teóricos (Eixo Y do gráfico normal)
n = len(df_effects)
# Fórmula para "Plotting Positions" (similar ao Minitab): (i - 0.5) / n
df_effects['Theoretical_Quantiles'] = stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n)

# Definir cores: Vermelho para significativo (p < 0.05), Azul para não significativo
# A apostila destaca os significativos fora da reta [cite: 1320]
df_effects['Color'] = ['red' if p < 0.05 else 'blue' for p in df_effects['PValue']]
df_effects['Label'] = [name if p < 0.05 else '' for name, p in zip(df_effects.index, df_effects['PValue'])]

# 4. Plotar o Gráfico
plt.figure(figsize=(10, 8))

# Desenhar a linha de referência (passando pelos pontos não significativos)
# Uma forma robusta é ligar o quartil 25% ao 75% dos dados, ou simplesmente ajustar uma reta
# nos pontos azuis (não significativos) que representam o erro puro.
if len(df_effects[df_effects['PValue'] >= 0.05]) > 1:
    # Ajuste linear apenas nos pontos não significativos (ruído) para traçar a linha de referência
    X_noise = df_effects[df_effects['PValue'] >= 0.05]['Effect']
    Y_noise = df_effects[df_effects['PValue'] >= 0.05]['Theoretical_Quantiles']
    z = np.polyfit(X_noise, Y_noise, 1)
    p = np.poly1d(z)
    
    # Extrapolar a linha para todo o range do gráfico
    x_range = np.linspace(df_effects['Effect'].min(), df_effects['Effect'].max(), 100)
    plt.plot(x_range, p(x_range), color='gray', linestyle='--', label='Distribuição Normal (Ruído)')

# Plotar os pontos
plt.scatter(df_effects['Effect'], df_effects['Theoretical_Quantiles'], 
            c=df_effects['Color'], s=60, edgecolor='k', zorder=3)

# Adicionar rótulos aos pontos significativos
for i, row in df_effects.iterrows():
    if row['PValue'] < 0.05:
        plt.text(row['Effect'], row['Theoretical_Quantiles'] + 0.1, i, 
                 fontsize=9, fontweight='bold', color='darkred')

# Estilização igual à da apostila (Figuras 26, 31)
plt.title('Normal Plot of the Effects\n(alpha = 0.05)', fontsize=14)
plt.xlabel('Efeito (Effect)', fontsize=12)
plt.ylabel('Quantil Normal / Probabilidade (Normal Score)', fontsize=12)
plt.axvline(0, color='black', linewidth=0.5)
plt.axhline(0, color='black', linewidth=0.5)
plt.grid(True, linestyle=':', alpha=0.6)

# Legenda manual para as cores
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='Significativo'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label='Não Significativo')
]
plt.legend(handles=legend_elements, loc='upper left')

plt.tight_layout()
plt.savefig('normal_plot_effects.png')
plt.show()

print("Gráfico 'Normal Plot of the Effects' gerado e salvo como 'normal_plot_effects.png'")