# Fatores do SEALER
fatores_0_sealer =    [ 12, 273+520,  18.5, 273+520, 10]
fatores_coef_sealer = [0.1,     500,   0.5,     500,  1]

# Fatores da SubCrítica
fatores_0_subcritica =    [1.10, 273+520,  20, 273+50, 0.90]
fatores_coef_subcritica = [0.04,     500,   2,     50, 0.09]


matriz_planejamento_2k = [
   [ 0,  0,  0,  0,  0],
   [ 0,  0,  0,  0,  0],
   [ 0,  0,  0,  0,  0],
   [ 0,  0,  0,  0,  0],
   [+1, +1, +1, +1, +1],
   [-1, +1, +1, +1, +1],
   [+1, -1, +1, +1, +1],
   [-1, -1, +1, +1, +1],
   [+1, +1, -1, +1, +1],
   [-1, +1, -1, +1, +1],
   [+1, -1, -1, +1, +1],
   [-1, -1, -1, +1, +1],
   [+1, +1, +1, -1, +1],
   [-1, +1, +1, -1, +1],
   [+1, -1, +1, -1, +1],
   [-1, -1, +1, -1, +1],
   [+1, +1, -1, -1, +1],
   [-1, +1, -1, -1, +1],
   [+1, -1, -1, -1, +1],
   [-1, -1, -1, -1, +1],
   [+1, +1, +1, +1, -1],
   [-1, +1, +1, +1, -1],
   [+1, -1, +1, +1, -1],
   [-1, -1, +1, +1, -1],
   [+1, +1, -1, +1, -1],
   [-1, +1, -1, +1, -1],
   [+1, -1, -1, +1, -1],
   [-1, -1, -1, +1, -1],
   [+1, +1, +1, -1, -1],
   [-1, +1, +1, -1, -1],
   [+1, -1, +1, -1, -1],
   [-1, -1, +1, -1, -1],
   [+1, +1, -1, -1, -1],
   [-1, +1, -1, -1, -1],
   [+1, -1, -1, -1, -1],
   [-1, -1, -1, -1, -1]
]

def conv_matriz_real(matriz_orig, fatores_0, fatores_coef):
    # Validação das dimensões (colunas)
    num_colunas = len(matriz_orig[0])
    if len(fatores_0) != num_colunas or len(fatores_coef) != num_colunas:
        raise ValueError(f"Erro: A matriz tem {num_colunas} colunas, mas as listas de fatores têm tamanhos diferentes.")

	# Calcular coeficiente
    matriz_real = []
    for i in range(len(matriz_orig)):
        nova_linha = []
        # 4. Loop pelas colunas (j)
        for j in range(num_colunas):
            # Aplica a fórmula: Valor * Coeficiente + Ponto Central
            valor_transformado = (matriz_orig[i][j] * fatores_coef[j]) + fatores_0[j]
            nova_linha.append(valor_transformado)
        matriz_real.append(nova_linha)

    return matriz_real

matriz_real_subcritica        = conv_matriz_real(matriz_planejamento_2k,fatores_0_subcritica,fatores_coef_subcritica)
matriz_real_sealer            = conv_matriz_real(matriz_planejamento_2k,fatores_0_sealer,fatores_coef_sealer)