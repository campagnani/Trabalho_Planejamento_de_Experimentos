#######################################################################
####                                                               ####
####       CENTRO DE DESENVOLVIMENTO DA TECNOLOGIA NUCLEAR         ####
####           Thabalho de planejamento de experimentos            ####
####                Thalles Oliveira Campagnani                    ####
####                                                               ####
#######################################################################

import pprint
import branch_calculations
import os
import libOpenSealer

m = branch_calculations.matriz_real_sealer #apenas para simplificar



libOpenSealer.dir("results_sealer")
sealer = libOpenSealer.SealerArctic(particulas=50000,ciclos=440,inativo=40,atrasados=True)
#sealer.materiais(
#    enriquecimento=12,
#    tempComb=273+520, 
#    densidadeCombUO2=18.5, 
#    tempRefri=273+520, 
#    densidadeRefrigerante=10
#    )
#sealer.geometria()
#sealer.run()
#exit(0)

#Simular todos os casos
vetor_keff = []
vetor_keff_incerteza = []
for i in range(0, len(m)):
    print("################")
    print("####### ",i, " ######")
    print("################")
    sealer.materiais(
            enriquecimento=m[i][0],
            tempComb=m[i][1], 
            densidadeCombUO2=m[i][2], 
            tempRefri=m[i][3], 
            densidadeRefrigerante=m[i][4])
    sealer.geometria()
    sealer.settings.seed=i+1
    sealer.settings.export_to_xml()
    sealer.run()
    os.rename(f"statepoint.{sealer.settings.batches}.h5", f"statepoint.{i}.h5")

#Gerar vetor keff
for i in range(0, len(m)):
    sp = libOpenSealer.openmc.StatePoint(f"statepoint.{i}.h5")
    vetor_keff.append(sp.keff.n)
    vetor_keff_incerteza.append(sp.keff.s)
    sp.close()



#Escrever resultados no arquivo
config = [sealer.settings.particles, sealer.settings.batches, sealer.settings.inactive]
with open(f"resultados_experimentos_{branch_calculations.fatores_0_sealer}_{branch_calculations.fatores_coef_sealer}_{config}.py", "w") as f:
    f.write("# Resultados da Simulação OpenMC\n")
    f.write("# Arquivo gerado automaticamente\n\n")
    
    # Entradas
    f.write("config = ")
    pprint.pprint(config, stream=f)
    f.write("\n")

    f.write("fatores_0 = ")
    pprint.pprint(branch_calculations.fatores_0_sealer, stream=f)
    f.write("\n")

    f.write("fatores_coef = ")
    pprint.pprint(branch_calculations.fatores_coef_sealer, stream=f)
    f.write("\n")

    #Saídas
    f.write("vetor_keff = ")
    pprint.pprint(vetor_keff, stream=f)
    f.write("\n")
    
    f.write("vetor_keff_incerteza = ")
    pprint.pprint(vetor_keff_incerteza, stream=f)
    f.write("\n")

    #Matiz padrão
    f.write("matriz_planejamento_2k = ")
    pprint.pprint(branch_calculations.matriz_planejamento_2k, stream=f)
    f.write("\n")

    #Matriz gerada
    f.write("matriz_real = ")
    pprint.pprint(branch_calculations.matriz_real_sealer, stream=f)
    f.write("\n")
    
