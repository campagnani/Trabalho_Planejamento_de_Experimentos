import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import argparse

args = ""
def args_command_line():
    # --- CONFIGURAÇÃO DOS ARGUMENTOS DA LINHA DE COMANDO ---
    parser = argparse.ArgumentParser(description='Análise de DOE com controle de gráficos.')
    parser.add_argument('--show', action='store_true', help='Exibe os gráficos na tela')
    parser.add_argument('--savePDF', action='store_true', help='Salva os gráficos em formato PDF')
    parser.add_argument('--savePNG', action='store_true', help='Salva os gráficos em formato PNG')
    return parser.parse_args()



# Função auxiliar para salvar/mostrar
def finalizar_grafico(args, fig, nome_arquivo):
    if args.savePDF:
        fig.savefig(f"{nome_arquivo}.pdf", bbox_inches='tight')
        print(f"Salvo: {nome_arquivo}.pdf")
    if args.savePNG:
        fig.savefig(f"{nome_arquivo}.png", bbox_inches='tight', dpi=300)
        print(f"Salvo: {nome_arquivo}.png")
    
    if args.show:
        plt.show()
    else:
        plt.close(fig) # Fecha para liberar memória se não for exibir


def reduc_data(df_cod, df_real, redu_p_central=0, reduc_fatorial=0):
    """
    Reduz o conjunto de dados experimental removendo pontos centrais ou 
    transformando o Fatorial Completo em Fracionado.
    
    Args:
        df_cod: DataFrame com variáveis codificadas (A, B, C...)
        df_real: DataFrame com variáveis reais
        redu_p_central (int): Qtd de pontos centrais a remover.
        reduc_fatorial (int): 
            0 = Mantém completo
            1 = Fração 1/2 (E = ABCD)
            2 = Fração 1/4 (D = AB, E = BC)
            
    Returns:
        [df_cod_new, df_real_new]
    """
    # 1. Identificação das Colunas Codificadas (A, B, C...)
    # Assume que as colunas codificadas são aquelas presentes no df_cod menos a 'Reatividade'
    # ou fixamos nas primeiras colunas. Pelo seu código anterior, são A, B, C, D, E.
    cols_cod = [c for c in df_cod.columns if c != 'Reatividade']
    
    # Cópia para não alterar os originais fora da função
    df_c = df_cod.copy()
    df_r = df_real.copy()

    # =========================================================================
    # PARTE A: Redução de Pontos Centrais
    # =========================================================================
    if redu_p_central > 0:
        # Identifica índices onde todas as colunas de fatores são 0
        mask_center = (df_c[cols_cod] == 0).all(axis=1)
        indices_centro = df_c[mask_center].index
        
        # Seleciona os N primeiros índices para remover
        # Se pedir para remover 4 e tiver 4, remove todos.
        if len(indices_centro) > 0:
            to_drop = indices_centro[:redu_p_central]
            df_c = df_c.drop(to_drop)
            df_r = df_r.drop(to_drop)
            print(f"[INFO] {len(to_drop)} pontos centrais removidos.")
        else:
            print("[AVISO] Nenhum ponto central encontrado para remover.")

    # =========================================================================
    # PARTE B: Redução Fatorial (Fracionamento)
    # =========================================================================
    if reduc_fatorial > 0:
        # Nota: As relações de definição (ex: E = ABCD) funcionam matematicamente
        # para os níveis +/- 1. Para o ponto central (0,0,0...), a relação 
        # 0 = 0*0*0*0 também é verdadeira, então os pontos centrais restantes 
        # (se houver) serão preservados automaticamente pelo filtro.

        mask_keep = None
        
        if reduc_fatorial == 1:
            # Transformar em 2^(5-1): Gerador E = A * B * C * D
            # Mantemos apenas as linhas onde essa igualdade é verdadeira
            condicao = df_c['E'] == (df_c['A'] * df_c['B'] * df_c['C'] * df_c['D'])
            mask_keep = condicao
            print("[INFO] Reduzindo para Fatorial Fracionado 2^(5-1) (E=ABCD).")
            
        elif reduc_fatorial == 2:
            # Transformar em 2^(5-2): Geradores D = A*B  e  E = B*C
            condicao_D = df_c['D'] == (df_c['A'] * df_c['B'])
            condicao_E = df_c['E'] == (df_c['B'] * df_c['C'])
            mask_keep = condicao_D & condicao_E
            print("[INFO] Reduzindo para Fatorial Fracionado 2^(5-2) (D=AB, E=BC).")
            
        # Aplica o filtro
        if mask_keep is not None:
            df_c = df_c[mask_keep]
            # Sincroniza o df_real usando os índices que sobraram no df_c
            df_r = df_r.loc[df_c.index]

    return df_c, df_r




def plot_pareto_normal(effects_data,p_values, model, title_suffix, standardized=False, filename_prefix=""):
    # Preparar dados
    df_plot = pd.DataFrame({'Effect': effects_data})
    df_plot['Abs_Effect'] = df_plot['Effect'].abs()
    df_plot = df_plot.sort_values(by='Abs_Effect')
    
    # Cores: Significativo se p < 0.05
    p_vals_sorted = p_values.reindex(df_plot.index)
    colors = ['red' if p < 0.05 else 'blue' for p in p_vals_sorted]
    
    # 1. Pareto
    fig_par = plt.figure(figsize=(10, 8))
    plt.barh(df_plot.index, df_plot['Abs_Effect'], color=colors)
    if standardized:
        # Linha de corte do t-student
        plt.axvline(stats.t.ppf(0.975, model.df_resid), color='k', linestyle='--', label='t-crítico (0.05)')
        plt.legend()
    plt.title(f'Gráfico de Pareto dos Efeitos {title_suffix}')
    plt.xlabel('Magnitude do Efeito')
    plt.tight_layout()
    finalizar_grafico(args, fig_par, f"Pareto_{filename_prefix}")
    
    # 2. Normal Probability Plot
    df_norm = pd.DataFrame({'Effect': effects_data, 'PValue': p_values}).sort_values(by='Effect')
    n = len(df_norm)
    df_norm['Theoretical'] = stats.norm.ppf((np.arange(1, n + 1) - 0.5) / n)
    df_norm['Color'] = ['red' if p < 0.05 else 'blue' for p in df_norm['PValue']]
    
    fig_norm = plt.figure(figsize=(8, 8))
    plt.scatter(df_norm['Effect'], df_norm['Theoretical'], c=df_norm['Color'], s=50, edgecolors='k')
    
    # Linha de ruído (ajustada pelos não significativos)
    nonsig = df_norm[df_norm['PValue'] >= 0.05]
    if len(nonsig) > 2:
        z = np.polyfit(nonsig['Effect'], nonsig['Theoretical'], 1)
        p = np.poly1d(z)
        x_range = np.linspace(df_norm['Effect'].min(), df_norm['Effect'].max(), 100)
        plt.plot(x_range, p(x_range), 'k--', alpha=0.5, label='Ruído')
    
    # Labels
    for idx, row in df_norm.iterrows():
        if row['PValue'] < 0.05:
            # Substituir ':' por '*' para ficar mais legível no gráfico
            lbl = idx.replace(':', '*')
            plt.text(row['Effect'], row['Theoretical'], lbl, fontsize=8, ha='right', color='darkred')
            
    plt.title(f'Gráfico Normal dos Efeitos {title_suffix}')
    plt.xlabel('Efeito')
    plt.ylabel('Probabilidade Normal (%)')
    plt.axhline(0, color='k', lw=0.5)
    plt.axvline(0, color='k', lw=0.5)
    plt.grid(True)
    finalizar_grafico(args, fig_norm, f"NormalPlot_{filename_prefix}")