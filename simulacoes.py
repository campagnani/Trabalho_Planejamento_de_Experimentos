#######################################################################
####                                                               ####
####       CENTRO DE DESENVOLVIMENTO DA TECNOLOGIA NUCLEAR         ####
####           Thabalho de planejamento de experimentos            ####
####                Thalles Oliveira Campagnani                    ####
####                                                               ####
#######################################################################

import libChicagoDenR1
libChicagoDenR1.simu = True


chicago = libChicagoDenR1.ChigagoDenR1(altura_fonte=None, particulas=50000, ciclos=100, inativo=10)

chicago.u_nat(enriquecimento=1.1,tempCombustivel=294,tempModerador=294, densidadeCombustivel=19.0, densidadeModerador=1.0)
chicago.run()