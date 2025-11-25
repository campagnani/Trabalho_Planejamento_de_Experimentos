#######################################################################
####                                                               ####
####       CENTRO DE DESENVOLVIMENTO DA TECNOLOGIA NUCLEAR         ####
####           Thabalho de planejamento de experimentos            ####
####                Thalles Oliveira Campagnani                    ####
####                                                               ####
#######################################################################

import pprint
import branch_calculations
import libChicagoDenR1
libChicagoDenR1.simu = True


m = branch_calculations.matriz_real #apenas para simplificar

chicago = libChicagoDenR1.ChigagoDenR1(altura_fonte=None, particulas=100000, ciclos=220, inativo=20)
#teste manual
#chicago.u_nat(enriquecimento=1.11,tempCombustivel=273+50, densidadeCombustivel=18, tempModerador=273+50, densidadeModerador=0.95)
#chicago.run()


#Simular todos os casos
vetor_keff = []
vetor_keff_incerteza = []
for i in range(0, len(m)):
    print("################")
    print("####### ",i, " ######")
    print("################")
    chicago.u_nat(enriquecimento=m[i][0],tempCombustivel=m[i][1], densidadeCombustivel=m[i][2], tempModerador=m[i][3], densidadeModerador=m[i][4])
    chicago.run()
    sp = libChicagoDenR1.openmc.StatePoint("statepoint.220.h5")
    vetor_keff.append(sp.keff.n)
    vetor_keff_incerteza.append(sp.keff.s)
    sp.close()


#Escrever resultados no arquivo
with open("resultados_experimentos.py", "w") as f:
    f.write("# Resultados da Simulação OpenMC\n")
    f.write("# Arquivo gerado automaticamente\n\n")
    
    f.write("vetor_keff = ")
    pprint.pprint(vetor_keff, stream=f)
    f.write("\n")
    
    f.write("vetor_keff_incerteza = ")
    pprint.pprint(vetor_keff_incerteza, stream=f)
    f.write("\n")

    f.write("matriz_planejamento_2k = ")
    pprint.pprint(branch_calculations.matriz_planejamento_2k, stream=f)
    f.write("\n")

    f.write("fatores_0 = ")
    pprint.pprint(branch_calculations.fatores_0, stream=f)
    f.write("\n")

    f.write("fatores_coef = ")
    pprint.pprint(branch_calculations.fatores_coef, stream=f)
    f.write("\n")

    f.write("matriz_real = ")
    pprint.pprint(branch_calculations.matriz_real, stream=f)
    f.write("\n")
    
