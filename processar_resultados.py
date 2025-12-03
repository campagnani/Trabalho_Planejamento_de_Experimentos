import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf  # Importação necessária para ANOVA funcionar
import seaborn as sns
from scipy import stats
from itertools import combinations
from statsmodels.stats.outliers_influence import variance_inflation_factor
import string
import libProcFatorial



# Variável auxiliar: Se qualquer opção de gráfico foi escolhida, entramos nos blocos de plotagem
args = libProcFatorial.args_command_line()
libProcFatorial.args = args
GERAR_GRAFICOS = args.show or args.savePDF or args.savePNG

# Importar dados do arquivo do usuário
import results_sealer_20251202_054041.resultados_experimentos__12_1273_18P5_1273_10___0P1_500_0P25_500_0P75___50000_50_10_ as r

# =============================================================================
# 1. PREPARAÇÃO DOS DADOS
# =============================================================================

# Cálculo da Reatividade (pcm)
# Considera os primeiros 4 pontos como pontos centrais (conforme padrão do arquivo)
pontos_centrais_keff = r.vetor_keff[:4]
media_keff_centro = np.mean(pontos_centrais_keff)

reatividade = []
for k in r.vetor_keff:
    # Conversão para pcm em relação à média do centro
    # rho = (k - k_ref) * 10^5
    rho = (k - media_keff_centro) * 10**5
    reatividade.append(rho)

# Nomes dos Fatores
factor_names = ['Enriq', 'TempComb', 'DensComb', 'TempRefri', 'DensRefri'] #deveria ser "r.factor_names"

# Gera automaticamente ['A', 'B', 'C', 'D', ...] baseado na quantidade de fatores
# Se len(factor_names) for 3, gera ['A', 'B', 'C']
cols_codificadas = list(string.ascii_uppercase[:len(factor_names)])

# DataFrame Codificado
df_cod = pd.DataFrame(r.matriz_planejamento_2k, columns=cols_codificadas)
df_cod['Reatividade'] = reatividade

# DataFrame Real (Não Codificado)
df_real = pd.DataFrame(r.matriz_real, columns=factor_names)
df_real['Reatividade'] = reatividade


### AVISO: Isso abaixo é útil apenas para o trabalho, para analisar a teoria (parte do trabalho)!! Não tem utilidade prática!! Remover posteriormente...
# Lembrete: a redução de experimentos tem que ser realizada ANTES do experimento ser realizado...
# A não ser que seja necessário provar que a teoria funciona de forma rápida, por exemplo para o Petri...
redu_p_central = input("\nReduzir pontos centrais: ").strip()
reduc_fatorial = input("\nReduzir fatoriais: ").strip()
[df_cod, df_real] = libProcFatorial.reduc_data(df_cod, df_real, int(redu_p_central), int(reduc_fatorial))


# Identificar pontos fatoriais e centrais
is_center = (df_cod[cols_codificadas] == 0).all(axis=1)
df_fatorial = df_cod[~is_center].copy()
df_centro = df_cod[is_center].copy()


# =============================================================================
# 1.5 CÁLCULO PRELIMINAR DE TODOS OS EFEITOS (PARA OS GRÁFICOS)
# =============================================================================
# Ajustamos um modelo completo (até 5ª ordem) apenas para extrair os efeitos
# para o Pareto e Gráfico Normal.
ordem_pre = input("Entre com a ordem dos efeitos prelimináres a serem analizados: ")
formula_maxima = f"Reatividade ~ (A + B + C + D + E)**{ordem_pre}"
print(formula_maxima)
model_preliminar = smf.ols(formula_maxima, data=df_cod).fit()

# Cálculos para os Gráficos Preliminares
effects_pre = model_preliminar.params.drop('Intercept') * 2
p_values_pre = model_preliminar.pvalues.drop('Intercept')
t_values_pre = model_preliminar.tvalues.drop('Intercept')
std_effects_pre = t_values_pre.abs()





# =============================================================================
# 3. PARTE 1 - GRÁFICOS PRELIMINARES
# =============================================================================

if GERAR_GRAFICOS:
    # --- 3.1 Gráfico de Efeitos Principais (com Ponto Central * Marrom) ---
    fig1, axes1 = plt.subplots(1, 5, figsize=(20, 5), sharey=True)
    global_mean_factorial = df_fatorial['Reatividade'].mean()
    mean_center_point = df_centro['Reatividade'].mean()

    # Calcular limites globais para o eixo Y
    y_vals = []
    for col in cols_codificadas:
        y_vals.extend(df_fatorial.groupby(col)['Reatividade'].mean().values)
    y_vals.append(mean_center_point)
    y_min, y_max = min(y_vals), max(y_vals)
    margin = (y_max - y_min) * 0.1

    for i, col in enumerate(cols_codificadas):
        ax = axes1[i]
        # Médias nos níveis -1 e +1
        means = df_fatorial.groupby(col)['Reatividade'].mean()
        
        # Plotar a linha de efeito (-1 a +1)
        ax.plot([-1, 1], [means[-1], means[1]], marker='o', linestyle='-', color='blue', label='Efeito Linear')
        
        # Plotar a linha de referência (Média Global Fatorial)
        ax.axhline(global_mean_factorial, color='gray', linestyle=':', linewidth=0.8)
        
        # Plotar Ponto Central (* Marrom)
        # Ele fica em x=0. Se estiver fora da linha azul, indica curvatura.
        ax.plot(0, mean_center_point, marker='*', color='brown', markersize=12, label='Ponto Central', linestyle='None')
        
        ax.set_title(factor_names[i])
        ax.set_xticks([-1, 0, 1])
        ax.set_xlim(-1.5, 1.5)
        ax.grid(True, linestyle=':', alpha=0.5)
        
        if i == 0:
            ax.set_ylabel('Média de Reatividade (pcm)')
            ax.legend(loc='best', fontsize='small')

    fig1.suptitle('Gráfico de Efeitos Principais (com Curvatura)', fontsize=16)
    plt.ylim(y_min - margin, y_max + margin)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # CHAMA A FUNÇÃO DE SALVAR/MOSTRAR
    libProcFatorial.finalizar_grafico(args, fig1, "Grafico_Efeitos_Principais")

    # --- 3.2 Gráfico de Interação ---
    fig2, axes2 = plt.subplots(5, 5, figsize=(18, 18), sharey=True)

    for i in range(5): # Linha (Fator Legenda)
        for j in range(5): # Coluna (Eixo X)
            ax = axes2[i, j]
            if i == j:
                ax.text(0.5, 0.5, factor_names[i], ha='center', va='center', fontsize=12, weight='bold')
                ax.set_axis_off()
            else:
                fi = cols_codificadas[i] # Fator que define as linhas (cor)
                fj = cols_codificadas[j] # Fator que define o eixo X
                
                # Calcular médias
                means = df_fatorial.groupby([fj, fi])['Reatividade'].mean().unstack()
                
                ax.plot([-1, 1], means[-1], 'b-o', label=f'{factor_names[i]} -1')
                ax.plot([-1, 1], means[1], 'r--s', label=f'{factor_names[i]} +1')
                
                ax.set_xticks([-1, 1])
                ax.grid(True, alpha=0.5)
                
                if j == 0: ax.set_ylabel('Reatividade')

    # Legenda global
    handles, labels = axes2[1, 0].get_legend_handles_labels()
    fig2.legend(handles, ['Nível -1', 'Nível +1'], loc='upper right', title="Nível do Fator")
    fig2.suptitle('Matriz de Interação', fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    libProcFatorial.finalizar_grafico(args, fig2, "Grafico_Interacoes")

    # --- 3.3 e 3.5 Gráficos Normais e de Pareto ---
    libProcFatorial.plot_pareto_normal(effects_pre, p_values=p_values_pre, model=model_preliminar, title_suffix="(Reatividade)")
    libProcFatorial.plot_pareto_normal(std_effects_pre, p_values=p_values_pre, model=model_preliminar, title_suffix="(Padronizados)", standardized=True)




# =============================================================================
# 3. PARTE 2 - MODELO INTERATIVO E ANOVA
# =============================================================================

print("\n" + "="*80)
print("DEFINIÇÃO DO MODELO")
print("="*80)
print("Escolha a ordem do modelo para análise:")
print("1 : Apenas efeitos principais (1ª ordem)")
print("2 : Efeitos principais e interações de 2ª ordem")
print("3 : Efeitos principais, 2ª e 3ª ordem")
print("4 : Efeitos principais até 4ª ordem")
print("5 : Modelo Completo (até 5ª ordem)")
print("0 : Personalizado (digite os termos ex: 'A B AB C DE')")
print("a : Automático (somente os termos de mais significância)")
print("x : Não fazer ANOVA")

choice = input("\nDigite sua opção: ").strip()

formula_str = ""
if choice == '0':
    print("Exemplos: 'A B AB' (Interações) ou 'A I(A**2)' (Quadráticos)")
    terms_input = input("Digite os termos separados por espaço: ").strip().split()
    
    terms_formatted = []
    for t in terms_input:
        # 1. Se for termo especial protegido por I(...) ou potência, mantém como está
        if "I(" in t or "**" in t:
            terms_formatted.append(t)
        # 2. Se for interação compacta ex: AB -> transforma em A:B
        elif len(t) > 1 and ':' not in t:
            t_split = ":".join(list(t))
            terms_formatted.append(t_split)
        # 3. Termos simples (A, B, C...)
        else:
            terms_formatted.append(t)
            
    formula_str = "Reatividade ~ " + " + ".join(terms_formatted)

elif choice.isdigit() and int(choice) in [1, 2, 3, 4, 5]:
    degree = int(choice)
    formula_str = "Reatividade ~ (" + " + ".join(cols_codificadas) + f")**{degree}"

elif choice=="a":
    exit(0)
elif choice=="x":
    exit(0)

else:
    print("Opção inválida. Usando modelo de 2ª ordem padrão.")
    formula_str = "Reatividade ~ (" + " + ".join(cols_codificadas) + ")**2"

print(f"\nFórmula utilizada: {formula_str}")


# --- AJUSTE DO MODELO (CODIFICADO) ---
model = smf.ols(formula_str, data=df_cod).fit()

# Cálculos para Efeitos
# Em DOE 2 níveis codificado, Efeito = 2 * Coeficiente
effects = model.params.drop('Intercept') * 2
t_values = model.tvalues.drop('Intercept')
p_values = model.pvalues.drop('Intercept')

# Efeitos Padronizados (t-values absolutos)
std_effects = t_values.abs()


# --- AJUSTE DO MODELO (NÃO CODIFICADO - REAIS) ---
# Precisamos adaptar a fórmula para usar as mesmas interações mas com o dataframe real
# é preciso converter os títulos das colunas do dataframe
# --- AJUSTE DO MODELO (NÃO CODIFICADO - REAIS) ---
import re  # Importando regex para substituição segura

# Mapeamento: {'A': 'Enriq', 'B': 'TempComb', ...}
mapping = dict(zip(cols_codificadas, factor_names))

# Função auxiliar que será chamada pelo regex para trocar pelo nome correto
def replace_var(match):
    return mapping[match.group(0)]

# Cria um padrão que busca EXATAMENTE as variáveis isoladas (A, B, C, D ou E)
# \b significa "borda da palavra", garantindo que não pegue letras dentro de outras palavras
pattern = r'\b(' + '|'.join(cols_codificadas) + r')\b'

# Faz a substituição em uma única passada (evita substituir o E de Enriq, por exemplo)
formula_real_str = re.sub(pattern, replace_var, formula_str)

print(f"Fórmula Reais: {formula_real_str}")

# Agora passamos a fórmula com os nomes CORRETOS para o df_real
model_real = smf.ols(formula_real_str, data=df_real).fit()




# =============================================================================
# 4. PARTE 2 - MODELO E ANOVA
# =============================================================================

print("\n" + "="*80)
print("ANÁLISE DO MODELO DE REGRESSÃO (REATIVIDADE)")
print("="*80)

# Plotar novo pareto
libProcFatorial.plot_pareto_normal(effects, p_values=p_values, model=model, title_suffix="(Reatividade)")
libProcFatorial.plot_pareto_normal(std_effects, p_values=p_values, model=model, title_suffix="(Padronizados)", standardized=True)

# --- 4.1 Equações ---

# Codificada
intercept_cod = model.params['Intercept']
eq_cod = f"Reatividade = {intercept_cod:.2f}"
for term, coef in model.params.items():
    if term == 'Intercept': continue
    sign = "+" if coef >= 0 else "-"
    # Substituir : por * para visualização estilo equação (ex: A:B -> A * B)
    term_clean = term.replace(':', ' * ')
    eq_cod += f"\n             {sign} {abs(coef):.2f} * {term_clean}"
print("\n--- Equação em Unidades Codificadas ---")
print(eq_cod)

# Não Codificada (Dinâmica - Baseada no model_real ajustado anteriormente)
# Agora usamos diretamente o model_real que foi criado com a fórmula traduzida (regex)
intercept_real = model_real.params['Intercept']
eq_real = f"Reatividade = {intercept_real:.2e}"

for term, coef in model_real.params.items():
    if term == 'Intercept': continue
    sign = "+" if coef >= 0 else "-"
    
    # O statsmodels usa ':' para interação nas fórmulas (ex: Enriq:TempComb)
    # Substituímos por ' * ' para ficar legível na equação final
    term_clean = term.replace(':', ' * ')
    
    eq_real += f"\n             {sign} {abs(coef):.2e} * {term_clean}"

print("\n--- Equação em Unidades Não Codificadas ---")
print(eq_real)

# --- 4.2 Tabela ANOVA ---
print("\n--- Tabela Análise de Variância (ANOVA) ---")
# Agora anova_lm funcionará porque usamos smf.ols (fórmula)
anova_table = sm.stats.anova_lm(model, typ=2) 

print(f"{'Fonte':<20} | {'GL':<5} | {'SQ (Aj.)':<12} | {'QM (Aj.)':<12} | {'Valor F':<10} | {'Valor-P':<10}")
print("-" * 85)
for idx, row in anova_table.iterrows():
    if idx == 'Residual': continue
    idx_clean = idx.replace(':', '*') # Limpar nome
    print(f"{idx_clean:<20} | {int(row['df']):<5} | {row['sum_sq']:<12.2f} | {row['sum_sq']/row['df']:<12.2f} | {row['F']:<10.2f} | {row['PR(>F)']:<10.4f}")

resid_row = anova_table.loc['Residual']
print(f"{'Erro':<20} | {int(resid_row['df']):<5} | {resid_row['sum_sq']:<12.2f} | {resid_row['sum_sq']/resid_row['df']:<12.2f} | {'':<10} | {'':<10}")

sst = model.centered_tss
dft = model.df_model + model.df_resid
print(f"{'Total':<20} | {int(dft):<5} | {sst:<12.2f} | {'':<12} | {'':<10} | {'':<10}")


# --- 4.3 Tabela Coeficientes Codificados ---
print("\n--- Tabela de Coeficientes Codificados ---")
# VIF Calculation - extraindo matrix exog do modelo de fórmula
exog_matrix = pd.DataFrame(model.model.exog, columns=model.model.exog_names)
vif_data = [variance_inflation_factor(exog_matrix.values, i) for i in range(exog_matrix.shape[1])]

print(f"{'Termo':<15} | {'Efeito':<10} | {'Coef':<10} | {'EP Coef':<10} | {'Valor-T':<10} | {'Valor-P':<10} | {'VIF':<10}")
print("-" * 90)
for i, term in enumerate(model.params.index):
    effect = model.params[term] * 2 if term != 'Intercept' else '-'
    coef = model.params[term]
    se = model.bse[term]
    t = model.tvalues[term]
    p = model.pvalues[term]
    vif = vif_data[i]
    
    eff_str = f"{effect:.2f}" if term != 'Intercept' else "-"
    term_clean = term.replace(':', '*') # Substituir : por *
    print(f"{term_clean:<15} | {eff_str:<10} | {coef:<10.2f} | {se:<10.2f} | {t:<10.2f} | {p:<10.4f} | {vif:<10.2f}")


# --- 4.4 Sumário do Modelo ---
def calculate_press_r2(model_obj):
    residuals = model_obj.resid
    hat_matrix_diag = model_obj.get_influence().hat_matrix_diag
    press = np.sum((residuals / (1 - hat_matrix_diag))**2)
    sst = np.sum((model_obj.model.endog - np.mean(model_obj.model.endog))**2)
    r2_pred = 1 - (press / sst)
    return r2_pred

r2 = model.rsquared
r2_adj = model.rsquared_adj
r2_pred = calculate_press_r2(model)
f_val = model.fvalue
p_val_model = model.f_pvalue

print("\n--- Sumário do Modelo ---")
print(f"R²          : {r2:.2%}")
print(f"R² Ajustado : {r2_adj:.2%}")
print(f"R² Predito  : {r2_pred:.2%}")
print(f"F-value     : {f_val:.2f}")
print(f"P-value     : {p_val_model:.4f}")


# =============================================================================
# 5. GRÁFICOS DE RESÍDUOS
# =============================================================================

if GERAR_GRAFICOS:
    residuos = model.resid
    ajustados = model.fittedvalues
    ordem = range(len(residuos))

    plt.figure(figsize=(14, 10))
    plt.suptitle("Gráficos de Resíduos para Reatividade", fontsize=16)

    # 1. Probabilidade Normal
    plt.subplot(2, 2, 1)
    stats.probplot(residuos, dist="norm", plot=plt)
    plt.title("Gráfico de Probabilidade Normal")
    plt.xlabel("Quantis Teóricos")
    plt.ylabel("Resíduos")

    # 2. Histograma
    plt.subplot(2, 2, 2)
    sns.histplot(residuos, kde=True, color='blue', edgecolor='k')
    plt.title("Histograma dos Resíduos")
    plt.xlabel("Resíduos")
    plt.ylabel("Frequência")

    # 3. Versus Ajustados
    plt.subplot(2, 2, 3)
    plt.scatter(ajustados, residuos, c='blue', edgecolor='k', alpha=0.7)
    plt.axhline(0, color='red', linestyle='--')
    plt.title("Resíduos vs. Ajustados")
    plt.xlabel("Valores Ajustados")
    plt.ylabel("Resíduos")

    # 4. Versus Ordem
    plt.subplot(2, 2, 4)
    plt.plot(ordem, residuos, marker='o', linestyle='-', color='blue', mec='k', alpha=0.7)
    plt.axhline(0, color='red', linestyle='--')
    plt.title("Resíduos vs. Ordem")
    plt.xlabel("Ordem de Observação")
    plt.ylabel("Resíduos")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()