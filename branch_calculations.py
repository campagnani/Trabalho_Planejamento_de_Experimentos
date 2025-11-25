
fatores_0 =    [1.11, 273+50,   18, 273+50, 0.95]

fatores_coef = [0.01,     50,  0.1,     50, 0.05]

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

matriz_real = conv_matriz_real(matriz_planejamento_2k,fatores_0,fatores_coef)