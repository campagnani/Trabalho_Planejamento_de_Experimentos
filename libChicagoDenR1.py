#!/usr/bin/python

from datetime import datetime
import openmc
import openmc.stats
import openmc.data
import numpy as np
import os
import math
from pprint import pprint
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker
import warnings

# Ignore warnings
warnings.filterwarnings("ignore")

os.system('clear')

simu = True

def mkdir(nome="teste_sem_nome",data=True,voltar=False):
    if (voltar==True):
        os.chdir("../")
    if (data==True):
        agora = datetime.now()
        nome = agora.strftime(nome+"_%Y%m%d_%H%M%S")
    if not os.path.exists(nome):
        os.makedirs(nome)
    os.chdir(nome)
    
def chdir(nome=None):
    if (nome != None):
        os.chdir(nome)
    else:
        diretorio_atual = os.getcwd()
        diretorios = [diretorio for diretorio in os.listdir(diretorio_atual) if os.path.isdir(os.path.join(diretorio_atual, diretorio))]

        data_mais_recente = 0
        pasta_mais_recente = None

        for diretorio in diretorios:
            data_criacao = os.path.getctime(os.path.join(diretorio_atual, diretorio))
            if data_criacao > data_mais_recente:
                data_mais_recente = data_criacao
                pasta_mais_recente = diretorio

        if pasta_mais_recente:
            os.chdir(os.path.join(diretorio_atual, pasta_mais_recente))
            print("Diretório mais recente encontrado:", pasta_mais_recente)
        else:
            print("Não foi possível encontrar um diretório mais recente.")


#####################################################
############ CARACTERÍSTICAS DA SUBCRÍTICA ##########
############          DO DEN-UFMG          ##########
#####################################################

class ChigagoDenR1:

    def __init__(self, material="u_nat", altura_fonte=0, particulas=1000, ciclos=100, inativo=10):
        #Definindo Material
        if(material=="u_nat"):
            self.material = self.u_nat()
        else:
            self.__del__(self)
        
        #Definindo geometria
        self.geometriaPadrao(altura_fonte)

        #Definindo simulação
        self.configuracoes(fonte=altura_fonte,particulas=particulas,ciclos=ciclos,inativo=inativo)
        
    def __del__(self):
        print(f"Objeto destruído.")

    def u_nat(self, enriquecimento=0.7, tempCombustivel=294,tempModerador=294, densidadeCombustivel=18.0, densidadeModerador=1.0):
        print("################################################")
        print("############ Definição dos Materiais ###########")
        print("############          U_nat          ###########")
        print("################################################")

        #self.combustivel = openmc.Material(name='Uránio Natural', material_id=1)
        #self.combustivel.add_nuclide('U234', 5.50000E-05, percent_type='ao')
        #self.combustivel.add_nuclide('U235', 7.20000E-03, percent_type='ao')
        #self.combustivel.add_nuclide('U238', 9.92745E-01, percent_type='ao')
        #self.combustivel.set_density('g/cm3', densidadeCombustivel)

        self.combustivel = openmc.Material(name='Uránio Natural', material_id=1)
        self.combustivel.add_element(element='U', percent=100, enrichment=enriquecimento)
        self.combustivel.set_density('g/cm3', densidadeCombustivel)
        self.combustivel.temperature=tempCombustivel

        self.moderador = openmc.Material(name='Água Leve', material_id=2)
        self.moderador.add_nuclide('H1' , 1.1187E-01 , percent_type='wo')
        self.moderador.add_nuclide('H2' , 3.3540E-05 , percent_type='wo')
        self.moderador.add_nuclide('O16', 8.8574E-01 , percent_type='wo')
        self.moderador.add_nuclide('O17', 3.5857E-04 , percent_type='wo')
        self.moderador.add_nuclide('O18', 1.9982E-03 , percent_type='wo')
        self.moderador.set_density('g/cm3', densidadeModerador)
        self.moderador.temperature=tempModerador
        
        self.ar = openmc.Material(name='Ar', material_id=3)
        self.ar.add_nuclide('N14' , 7.7826E-01 , percent_type='ao')
        self.ar.add_nuclide('N15' , 2.8589E-03 , percent_type='ao')
        self.ar.add_nuclide('O16' , 1.0794E-01 , percent_type='ao')
        self.ar.add_nuclide('O17' , 1.0156E-01 , percent_type='ao')
        self.ar.add_nuclide('O18' , 3.8829E-05 , percent_type='ao')
        self.ar.add_nuclide('Ar36', 2.6789E-03 , percent_type='ao')
        self.ar.add_nuclide('Ar38', 3.4177E-03 , percent_type='ao')
        self.ar.add_nuclide('Ar40', 3.2467E-03 , percent_type='ao')
        self.ar.set_density('g/cm3', 0.001225)

        self.aluminio = openmc.Material(name='Alúminio', material_id=4)
        self.aluminio.add_nuclide('Al27', 1 , percent_type ='wo')
        self.aluminio.set_density('g/cm3', 2.7)

        self.SS304 = openmc.Material(name='Aço INOX', material_id=5)
        self.SS304.add_element('C',  4.3000E-04 , percent_type = 'wo')
        self.SS304.add_element('Cr', 1.9560E-01 , percent_type = 'wo')
        self.SS304.add_element('Ni', 9.6600E-02 , percent_type = 'wo')
        self.SS304.add_element('Mo', 8.9000E-03 , percent_type = 'wo')
        self.SS304.add_element('Mn', 5.4000E-04 , percent_type = 'wo')
        self.SS304.add_element('Si', 5.0000E-04 , percent_type = 'wo')
        self.SS304.add_element('Cu', 2.0000E-05 , percent_type = 'wo')
        self.SS304.add_element('Co', 3.0000E-05 , percent_type = 'wo')
        self.SS304.add_element('P',  2.7000E-04 , percent_type = 'wo')
        self.SS304.add_element('S',  1.0000E-04 , percent_type = 'wo')
        self.SS304.add_element('N',  1.4000E-04 , percent_type = 'wo')
        self.SS304.add_element('Fe', 6.9687E-01 , percent_type = 'wo')
        self.SS304.set_density('g/cm3', 7.92)

        self.fonte = openmc.Material(name='Plutonio Berilio', material_id=6)
        self.fonte.add_nuclide('Pu239', percent=0.6204, percent_type='wo')
        self.fonte.add_nuclide('Pu240', percent=0.0453, percent_type='wo')
        self.fonte.add_nuclide('Pu241', percent=0.0015, percent_type='wo')
        self.fonte.add_nuclide('Be9',   percent=0.3328, percent_type='wo')
        self.fonte.set_density('g/cm3', 4.6781)

        self.materiais = openmc.Materials([self.combustivel,self.moderador,self.ar,self.aluminio,self.SS304,self.fonte,])
        self.materiais.cross_sections = '/opt/nuclear-data/endfb-viii.0-hdf5/cross_sections.xml' 
        if simu:
            self.materiais.export_to_xml()
        
        #Já definir as cores dos materiais para futuros plots
        self.colors = {
            self.moderador: 'blue',
            self.combustivel: 'yellow',
            self.ar: 'pink',
            self.aluminio: 'black',
            self.SS304: 'gray',
            self.fonte: 'red'
        }


    def geometriaPadrao(self, altura_fonte = 0):
        print("################################################")
        print("############ Definição da Geometria ############")
        print("################################################")

        #Fabricante: Nuclear Chicago Corporation
        #Código: NC-9000
        #Finalidade: pesquisa e fins didáticos.
        #Inauguração: 17 de dezembro de 1960 no Instituto Tecnológico de Aeronáutica.
        #Descomissionamento: 1990?
        ############ Características construtivas

        ######### Barra combustível ##########

            #Composição: urânio natural compactado
            #Revestimento: alumínio (0.1 cm de espessura)
        clad_combustivel_diametro_externo = 3.07  #cm
        clad_combustivel_diametro_interno = 1.27  #cm
            #Quantidade: 			1410
            #Densidade 
        Barra_combustivel_Densidade = 			18.0	#± 0,1 g/cm³.
            #Geometria: Barra anelar
            #Dimensões:
        Barra_combustivel_Comprimento       =   21.45   #± 0,5 cm.
        Barra_combustivel_Diametro_externo  =   2.87    #± 0,1 cm.
        Barra_combustivel_Diametro_interno  =   1.47    #± 0,1 cm.
            #Peso: 				1,867	Kg
                #Peso total: 			2,63	t

        #Vareta combustível:
            #Composição: Alumínio
            #Quantidade: 283 (282 com combustível e 1 com a fonte de nêutrons).
            #Dimensões:
        Vareta_combustivel_Espessura        = 	0.1     #cm
        Vareta_combustivel_Diametro_interno =   3.32    #± 0,2 cm.
        Vareta_combustivel_Comprimento      =   59*2.54  #59"

        #Interior:
        #Parte Superior: 5 barras de combustível em série e ar.
        #Parte Inferior: Água
        #Suporte/separador:
            #Material: alumínio
            #Altura:
            
        #Reticulado:
            #Formato: hexagonal
        Reticulado_Distancia=   5 #cm entre os centros de das barras mais próximas.
            
        #Tanque:
            #Material: aço inoxidável
            #Dimensões:
        Tanque_Altura       =	152.5
        self.tanque_altura  =   152.5
        Tanque_Diametro     =	122
        Tanque_Espessura    =	1.2

        ######### Moderador e refletor ##########

            #Tipo: água leve.
            #Obs:  sofre um contínuo tratamento por resinas de deionização. São usadas, em série, duas colunas de resina, uma para os aniontes e outra para os cationtes.
            #Nível: 1350 mm de altura do fundo do tanque
            #Espessura:
                #Fundo: 165mm (distância do fundo a parte mais baixa do comb.)
                #Lateral: 200mm (distância da última vareta a lateral do tanque)
                #Superior: 0mm (não é refletido).

        #Grade:
            #Função: definir o arranjo
            #Composição: alumínio
            #Dimensões:
        Grade_Espessura  =   1   # em altura
        Grade_ressalto   =   0.4 # Altura do ressalto / apoio da barra
        Grade_Diametro   =   101 #cm 
            #Quantidade: duas.
            #Posicionamento:
                #Primeira: fundo do tanque
        Grade_Posicionamento = 26.0        # Em relação a grade de baixo (cm)

        #Fonte de nêutrons:
            #Composição: Plutônio + Berílio
            #Localização: Vareta central
            #80g de plutonio total
            #Atividade de 5 curie total
            #3 cilindros de aço inox com plutônio-berilio
            #Colocados dentro de recipiente aluminio de 3 cm de diametro e 20cm de altura
            #
            #Atividades cada mini-cilindro:
            # 1.87 * 10^6
            # 3.53 * 10^6
            # 3.86 * 10^6

        #Diâmetro médio do núcleo: 800 mm

        ###
        ###
        ###

        ###
        # ZERO DA GEOMETRIA É O MEIO DO TANQUE
        ###

        #Variáveis de Dimenções dos planos
        fundo_tanque_inferior   = -Tanque_Altura/2
        fundo_tanque_superior   = fundo_tanque_inferior + Tanque_Espessura
        lateral_tanque_interna  = Tanque_Diametro - (2*Tanque_Espessura)
        altura_tanque           = fundo_tanque_inferior + Tanque_Altura
        nivel_dagua             = altura_tanque - 6*2.54
        grade_inferior_down     = fundo_tanque_superior + 2.54
        grade_interior_ressalto = grade_inferior_down + Grade_ressalto
        grade_inferior_up       = grade_inferior_down + Grade_Espessura
        grade_superior_down     = grade_inferior_up + Grade_Posicionamento
        grade_superior_up       = grade_superior_down + Grade_Espessura
        vareta_altura           = grade_interior_ressalto + Vareta_combustivel_Comprimento

        ##Planos internos a vareta
        refletor_interno_superior = grade_inferior_down + 6.5*2.54
        suporte_interno_superior = refletor_interno_superior + 0.2
        elemento_combustivel = suporte_interno_superior + Barra_combustivel_Comprimento*5

        clad_comb_1               = suporte_interno_superior + 0.5    #
        clad_comb_2               = clad_comb_1  + 20.45     # Fuel
        clad_comb_3               = clad_comb_2  + 0.5       
        clad_comb_4               = clad_comb_3  + 0.5       #
        clad_comb_5               = clad_comb_4  + 20.45     # Fuel
        clad_comb_6               = clad_comb_5  + 0.5  
        clad_comb_7               = clad_comb_6  + 0.5       #
        clad_comb_8               = clad_comb_7  + 20.45     # Fuel
        clad_comb_9               = clad_comb_8  + 0.5  
        clad_comb_10              = clad_comb_9  + 0.5       #
        clad_comb_11              = clad_comb_10 + 20.45     # Fuel
        clad_comb_12              = clad_comb_11 + 0.5  
        clad_comb_13              = clad_comb_12 + 0.5       #
        clad_comb_14              = clad_comb_13 + 20.45     # Fuel
        clad_comb_15              = clad_comb_14 + 0.5  
        self.centro_fonte = clad_comb_7 + 10.225
        #Criação das formas geométricas

        # Planos Horizontais
        plano_grade_superior_1          = openmc.ZPlane(z0=grade_superior_down,)
        plano_grade_superior_2          = openmc.ZPlane(z0=grade_superior_up,)
        plano_ressalto_grade            = openmc.ZPlane(z0=grade_interior_ressalto)
        plano_grade_inferior_1          = openmc.ZPlane(z0=grade_inferior_down,)
        plano_grade_inferior_2          = openmc.ZPlane(z0=grade_inferior_up,)
        plano_altura_tanque             = openmc.ZPlane(z0=altura_tanque,)
        plano_beirada_altura_tanque     = openmc.ZPlane(z0=altura_tanque - 1.2)
        plano_divisao_altura_tanque_sup = openmc.ZPlane(z0=altura_tanque - 76 + 1.2)
        plano_divisao_altura_tanque_inf = openmc.ZPlane(z0=altura_tanque - 76 - 1.2)
        plano_fundo_tanque_superior     = openmc.ZPlane(z0=fundo_tanque_superior,)
        plano_fundo_tanque_inferior     = openmc.ZPlane(z0=fundo_tanque_inferior,)
        plano_refletor_interno          = openmc.ZPlane(z0=refletor_interno_superior)
        plano_suporte_interno           = openmc.ZPlane(z0=suporte_interno_superior)
        plano_elemento_combustivel      = openmc.ZPlane(z0=elemento_combustivel)
        plano_vareta_altura             = openmc.ZPlane(z0=vareta_altura)
        plano_refletor_lateral_superior = openmc.ZPlane(z0=nivel_dagua)                       

        # Especial divisões no combustível
        plano_clad_comb_1               = openmc.ZPlane(z0=clad_comb_1)      #
        plano_clad_comb_2               = openmc.ZPlane(z0=clad_comb_2 )     # Fuel
        plano_clad_comb_3               = openmc.ZPlane(z0=clad_comb_3 )     
        plano_clad_comb_4               = openmc.ZPlane(z0=clad_comb_4 )     #
        plano_clad_comb_5               = openmc.ZPlane(z0=clad_comb_5 )     # Fuel
        plano_clad_comb_6               = openmc.ZPlane(z0=clad_comb_6 )
        plano_clad_comb_7               = openmc.ZPlane(z0=clad_comb_7 )     #
        plano_clad_comb_8               = openmc.ZPlane(z0=clad_comb_8 )     # Fuel
        plano_clad_comb_9               = openmc.ZPlane(z0=clad_comb_9 )
        plano_clad_comb_10              = openmc.ZPlane(z0=clad_comb_10)     #
        plano_clad_comb_11              = openmc.ZPlane(z0=clad_comb_11)     # Fuel
        plano_clad_comb_12              = openmc.ZPlane(z0=clad_comb_12)
        plano_clad_comb_13              = openmc.ZPlane(z0=clad_comb_13)     #
        plano_clad_comb_14              = openmc.ZPlane(z0=clad_comb_14)     # Fuel
        plano_clad_comb_15              = openmc.ZPlane(z0=clad_comb_15)


        #plano_suporte_interno_central   = openmc.ZPlane(z0=suporte_interno_central)
        #plano_refletor_interno_central  = openmc.ZPlane(z0=refletor_interno_central)

        # Superfícies de cilindros

        # Combustivel
        cilindro_raio_interno_combustivel = openmc.ZCylinder(r=Barra_combustivel_Diametro_interno/2)
        cilindro_raio_externo_combustivel = openmc.ZCylinder(r=Barra_combustivel_Diametro_externo/2)
        clad_raio_interno_combustivel     = openmc.ZCylinder(r=clad_combustivel_diametro_interno/2)
        clad_raio_externo_combustivel     = openmc.ZCylinder(r=clad_combustivel_diametro_externo/2)

        # Radial fora do combustivel
        cilindro_raio_interno_vareta    = openmc.ZCylinder(r=Vareta_combustivel_Diametro_interno/2)
        cilindro_raio_externo_vareta    = openmc.ZCylinder(r=Vareta_combustivel_Diametro_interno/2+Vareta_combustivel_Espessura)
        cilindro_raio_interno_tanque    = openmc.ZCylinder(r=lateral_tanque_interna/2)
        cilindro_raio_externo_tanque    = openmc.ZCylinder(r=Tanque_Diametro/2)
        cilindro_raio_externo_grade     = openmc.ZCylinder(r=Grade_Diametro/2)

        # Ressalto da grade
        cilindro_ressalto = openmc.ZCylinder(r=(Vareta_combustivel_Diametro_interno/2+Vareta_combustivel_Espessura)/2-0.2)

        # Beirada externa do tanque
        self.cilindro_beirada_tanque = openmc.ZCylinder(r=Tanque_Diametro/2+3.0)

        # Células Vareta
        self.celula_moderador1               = openmc.Cell(fill=self.moderador,   region=+plano_fundo_tanque_superior&-plano_grade_inferior_1&+cilindro_raio_externo_vareta)
        self.celula_moderador2               = openmc.Cell(fill=self.moderador,   region=+plano_grade_inferior_2&-plano_grade_superior_1&+cilindro_raio_externo_vareta)
        self.celula_moderador3               = openmc.Cell(fill=self.moderador,   region=+plano_grade_superior_2&-plano_refletor_lateral_superior&+cilindro_raio_externo_vareta)
        self.celula_refletor_interno         = openmc.Cell(fill=self.moderador,   region=+plano_ressalto_grade&-plano_refletor_interno&-cilindro_raio_interno_vareta
                                                          | +plano_grade_inferior_1&-plano_ressalto_grade&-cilindro_ressalto
                                                          | +plano_fundo_tanque_superior&-plano_grade_inferior_1&-cilindro_raio_externo_vareta)

        self.celula_clad_vareta              = openmc.Cell(fill=self.aluminio,    region=+plano_ressalto_grade&-plano_vareta_altura&+cilindro_raio_interno_vareta&-cilindro_raio_externo_vareta)
        self.celula_suporte_interno          = openmc.Cell(fill=self.aluminio,    region=+plano_refletor_interno&-plano_suporte_interno&-cilindro_raio_interno_vareta)
        self.celula_clad_combustivel_interno = openmc.Cell(fill=self.aluminio,    region=+plano_suporte_interno&-plano_elemento_combustivel&+clad_raio_interno_combustivel&-cilindro_raio_interno_combustivel)
        self.celula_clad_combustivel_externo = openmc.Cell(fill=self.aluminio,    region=+plano_suporte_interno&-plano_elemento_combustivel&+cilindro_raio_externo_combustivel&-clad_raio_externo_combustivel)

        self.celula_ar_interno_elemento      = openmc.Cell(fill=self.ar,          region=+plano_suporte_interno&-plano_elemento_combustivel&-clad_raio_interno_combustivel)
        self.celula_ar_externo_elemento      = openmc.Cell(fill=self.ar,          region=+plano_suporte_interno&-plano_elemento_combustivel&+clad_raio_externo_combustivel&-cilindro_raio_interno_vareta)
        self.celula_ar_superior_interno      = openmc.Cell(fill=self.ar,          region=+plano_elemento_combustivel&-plano_vareta_altura&-cilindro_raio_interno_vareta)
        self.celula_ar_externo_vareta        = openmc.Cell(fill=self.ar,          region=+plano_refletor_lateral_superior&-plano_vareta_altura&+cilindro_raio_externo_vareta)

        self.celula_combustivel = openmc.Cell(fill=self.combustivel, region=+plano_clad_comb_1&-plano_clad_comb_2&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel
                                                   | +plano_clad_comb_4&-plano_clad_comb_5&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel
                                                   | +plano_clad_comb_7&-plano_clad_comb_8&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel
                                                   | +plano_clad_comb_10&-plano_clad_comb_11&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel
                                                   | +plano_clad_comb_13&-plano_clad_comb_14&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)

        self.celula_clad_bottom = openmc.Cell(fill=self.aluminio, region=+plano_suporte_interno&-plano_clad_comb_1&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_1      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_1&-plano_clad_comb_2&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_2      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_2&-plano_clad_comb_3&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_3      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_3&-plano_clad_comb_4&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_4      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_5&-plano_clad_comb_6&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_5      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_6&-plano_clad_comb_7&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_6      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_8&-plano_clad_comb_9&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_7      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_9&-plano_clad_comb_10&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_8      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_11&-plano_clad_comb_12&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_9      = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_12&-plano_clad_comb_13&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)
        self.celula_clad_10     = openmc.Cell(fill=self.aluminio, region=+plano_clad_comb_14&-plano_clad_comb_15&+cilindro_raio_interno_combustivel&-cilindro_raio_externo_combustivel)

        # Celulas grade de espaçamento (construída de cima para baixo)
        self.celula_ressalto_inferior = openmc.Cell(fill=self.aluminio, region=+plano_grade_inferior_1&-plano_ressalto_grade&+cilindro_ressalto)
        self.celula_grade_inferior    = openmc.Cell(fill=self.aluminio, region=+plano_ressalto_grade&-plano_grade_inferior_2&+cilindro_raio_externo_vareta)        
        self.celula_grade_superior    = openmc.Cell(fill=self.aluminio, region=+plano_grade_superior_1&-plano_grade_superior_2&+cilindro_raio_externo_vareta)        

        ####################################################################################################################################
        ############################### Celulas e superfícies específicas para o Universo Vareta Central (fonte)############################
        ####################################################################################################################################
        if altura_fonte is not None:
        
            ##Planos internos a vareta central
            suporte_interno_central            = fundo_tanque_superior   + 5*2.54
            refletor_inferior_interno_central  = suporte_interno_central - 0.2
            limite_clad_fonte_inferior         = self.centro_fonte - 10 + altura_fonte
            self.limite_fonte_inferior         = limite_clad_fonte_inferior + 0.1
            self.limite_fonte_superior         = self.limite_fonte_inferior + 19.8
            limite_clad_fonte_superior         = self.limite_fonte_superior + 0.1
            self.diametro_cilindro_fonte       = 1.28433406196 # cm
            self.diametro_cilindro_ar_fonte    = 2.8 # cm
            diametro_clad_fonte                = 3.0 # cm

            plano_superior_suporte_fonte = openmc.ZPlane(z0=suporte_interno_central)
            plano_inferior_suporte_fonte = openmc.ZPlane(z0=refletor_inferior_interno_central)
            plano_superior_fonte         = openmc.ZPlane(z0=self.limite_fonte_superior)
            plano_inferior_fonte         = openmc.ZPlane(z0=self.limite_fonte_inferior)
            plano_clad_fonte_superior    = openmc.ZPlane(z0=limite_clad_fonte_superior)
            plano_clad_fonte_inferior    = openmc.ZPlane(z0=limite_clad_fonte_inferior)

            superficie_radial_fonte      = openmc.ZCylinder(r=self.diametro_cilindro_fonte/2)
            superficie_radial_ar_fonte   = openmc.ZCylinder(r=self.diametro_cilindro_ar_fonte/2)
            superficie_radial_clad_fonte = openmc.ZCylinder(r=diametro_clad_fonte/2)

            self.celula_refletor_inferior_suporte_fonte  = openmc.Cell(fill=self.moderador,   region=+plano_fundo_tanque_superior&-plano_inferior_suporte_fonte&-cilindro_raio_interno_vareta)
            self.celula_suporte_interno_fonte    = openmc.Cell(fill=self.aluminio,    region=+plano_inferior_suporte_fonte&-plano_superior_suporte_fonte&-cilindro_raio_interno_vareta)
            self.celula_fonte                    = openmc.Cell(fill=self.fonte,   region=+plano_inferior_fonte&-plano_superior_fonte&-superficie_radial_fonte)  
            self.celula_clad_fonte               = openmc.Cell(fill=self.aluminio,    region=+plano_inferior_fonte&-plano_superior_fonte&-superficie_radial_clad_fonte&+superficie_radial_ar_fonte
                                                                            | +plano_clad_fonte_inferior&-plano_inferior_fonte&-superficie_radial_clad_fonte
                                                                            | +plano_superior_fonte&-plano_clad_fonte_superior&-superficie_radial_clad_fonte)  
            self.celula_ar_fonte                 = openmc.Cell(fill=self.ar,          region=+plano_inferior_fonte&-plano_superior_fonte&-superficie_radial_ar_fonte&+superficie_radial_fonte)
            self.celula_refletor_lateral_fonte   = openmc.Cell(fill=self.moderador,   region=+plano_clad_fonte_inferior&-plano_clad_fonte_superior&-cilindro_raio_interno_vareta&+superficie_radial_clad_fonte)
            self.celula_refletor_superior_fonte  = openmc.Cell(fill=self.moderador,   region=+plano_clad_fonte_superior&-plano_refletor_lateral_superior&-cilindro_raio_interno_vareta)
            self.celula_refletor_inferior_fonte  = openmc.Cell(fill=self.moderador,   region=-plano_clad_fonte_inferior&+plano_superior_suporte_fonte&-cilindro_raio_interno_vareta)

            # CÉLULAS QUE FORAM DUPLICADAS PELO FATO DE NÃO PODER REUTILIZAR CÉLULAS PARA OUTROS UNIVERSOS
            self.celula_moderador1_fonte          = openmc.Cell(fill=self.moderador,   region=+plano_fundo_tanque_superior&-plano_grade_inferior_1&+cilindro_raio_externo_vareta)
            self.celula_moderador2_fonte          = openmc.Cell(fill=self.moderador,   region=+plano_grade_inferior_2&-plano_grade_superior_1&+cilindro_raio_externo_vareta)
            self.celula_moderador3_fonte          = openmc.Cell(fill=self.moderador,   region=+plano_grade_superior_2&-plano_refletor_lateral_superior&+cilindro_raio_externo_vareta)
            self.celula_grade_inferior_fonte      = openmc.Cell(fill=self.aluminio,    region=+plano_grade_inferior_1&-plano_grade_inferior_2&+cilindro_raio_externo_vareta)
            self.celula_grade_superior_fonte      = openmc.Cell(fill=self.aluminio,    region=+plano_grade_superior_1&-plano_grade_superior_2&+cilindro_raio_externo_vareta)
            self.celula_clad_vareta_fonte         = openmc.Cell(fill=self.aluminio,    region=+plano_fundo_tanque_superior&-plano_vareta_altura&+cilindro_raio_interno_vareta&-cilindro_raio_externo_vareta)
            self.celula_ar_externo_vareta_fonte   = openmc.Cell(fill=self.ar,          region=+plano_refletor_lateral_superior&-plano_vareta_altura&+cilindro_raio_externo_vareta)
            self.celula_ar_superior_interno_fonte = openmc.Cell(fill=self.ar,          region=+plano_elemento_combustivel&-plano_vareta_altura&-cilindro_raio_interno_vareta)


            self.universo_fonte = openmc.Universe(cells=(self.celula_clad_vareta_fonte,self.celula_refletor_inferior_fonte,self.celula_suporte_interno_fonte,
                                                    self.celula_fonte,self.celula_clad_fonte,self.celula_refletor_lateral_fonte, self.celula_ar_fonte,
                                                    self.celula_refletor_superior_fonte,self.celula_ar_externo_vareta_fonte,self.celula_refletor_inferior_suporte_fonte,
                                                    self.celula_grade_superior_fonte, self.celula_grade_inferior_fonte,self.celula_ar_superior_interno_fonte,
                                                    self.celula_moderador1_fonte, self.celula_moderador2_fonte, self.celula_moderador3_fonte,))  # Identidade da fonte
            
        ####################################################################################################################################
        ####################################################################################################################################

        ####################################################################################################################################
        ###############################     Celulas e superfícies específicas para o Universo agua (central)    ############################
        ####################################################################################################################################

        self.celula_agua               = openmc.Cell(fill=self.moderador,   region=+plano_fundo_tanque_superior&-plano_refletor_lateral_superior&-cilindro_raio_externo_vareta)
        self.celula_agua_2             = openmc.Cell(fill=self.moderador,   region=+plano_fundo_tanque_superior&-plano_grade_inferior_1&+cilindro_raio_externo_vareta)
        self.celula_agua_3             = openmc.Cell(fill=self.moderador,   region=+plano_grade_inferior_2&-plano_grade_superior_1&+cilindro_raio_externo_vareta)
        self.celula_agua_4             = openmc.Cell(fill=self.moderador,   region=+plano_grade_superior_2&-plano_refletor_lateral_superior&+cilindro_raio_externo_vareta)
        self.celula_ar_central         = openmc.Cell(fill=self.ar,          region=+plano_refletor_lateral_superior&-plano_vareta_altura&-cilindro_raio_externo_vareta)
        self.celula_ar_central_externa = openmc.Cell(fill=self.ar,          region=+plano_refletor_lateral_superior&-plano_vareta_altura&+cilindro_raio_externo_vareta)

        self.celula_central_grade_inferior  = openmc.Cell(fill=self.aluminio,    region=+plano_grade_inferior_1&-plano_grade_inferior_2&+cilindro_raio_externo_vareta)
        self.celula_central_grade_superior  = openmc.Cell(fill=self.aluminio,    region=+plano_grade_superior_1&-plano_grade_superior_2&+cilindro_raio_externo_vareta)
        
        self.universo_agua = openmc.Universe(cells=(self.celula_agua,self.celula_agua_2, self.celula_agua_3, self.celula_agua_4, self.celula_ar_central, 
                                                    self.celula_ar_central_externa, self.celula_central_grade_inferior, self.celula_central_grade_superior))  
                                                    # Identidade da agua

        ####################################################################################################################################
        ####################################################################################################################################

        #Celula para universo apenas com refletor
        self.celula_refletor_matriz          = openmc.Cell(fill=self.moderador, region=+plano_fundo_tanque_superior&-plano_grade_inferior_1
                                                                                     | +plano_grade_inferior_2&-plano_grade_superior_1
                                                                                     | +plano_grade_superior_2&-plano_refletor_lateral_superior                                                                                     )
        self.celula_grade_externa_inferior   = openmc.Cell(fill=self.aluminio,  region=+plano_grade_inferior_1&-plano_grade_inferior_2)
        self.celula_grade_externa_superior   = openmc.Cell(fill=self.aluminio,  region=+plano_grade_superior_1&-plano_grade_superior_2)
        self.celula_ar_externo               = openmc.Cell(fill=self.ar,        region=+plano_refletor_lateral_superior&-plano_vareta_altura)

        #Universos
        self.universo_vareta_combustível     = openmc.Universe(cells=(self.celula_clad_vareta,self.celula_refletor_interno, self.celula_suporte_interno, self.celula_ar_interno_elemento,
                                                                self.celula_combustivel,self.celula_ar_externo_elemento, self.celula_moderador1, self.celula_moderador2, self.celula_moderador3, 
                                                                self.celula_ar_superior_interno, self.celula_clad_combustivel_interno,self.celula_clad_combustivel_externo,
                                                                self.celula_clad_1, self.celula_clad_2,self.celula_clad_3,self.celula_clad_4,self.celula_clad_5,self.celula_clad_6,self.celula_clad_7,
                                                                self.celula_clad_8,self.celula_clad_9,self.celula_clad_10,self.celula_clad_bottom,self.celula_ar_externo_vareta,
                                                                self.celula_grade_superior, self.celula_grade_inferior, self.celula_ressalto_inferior))

        #self.universo_vareta_central         = openmc.Universe(cells=(self.celula_vareta,self.celula_refletor_fonte,self.celula_suporte_interno_fonte, \
        #                                                         self.celula_moderador_interno_fonte, self.celula_moderador))

        self.universo_agua_ar               = openmc.Universe(cells=(self.celula_refletor_matriz,self.celula_ar_externo,self.celula_grade_externa_inferior,
                                                                self.celula_grade_externa_superior))

        #Criação da Matriz Hexagonal

        matriz_hexagonal = openmc.HexLattice()
        matriz_hexagonal.center = (0., 0.)
        matriz_hexagonal.pitch = (2*2.54,)
        matriz_hexagonal.outer = self.universo_agua_ar
        matriz_hexagonal.orientation = 'x'

        #print(matriz_hexagonal.show_indices(num_rings=13))

        anel_misturado_10 = [self.universo_agua_ar]*4 + [self.universo_vareta_combustível] + [self.universo_agua_ar] + [self.universo_vareta_combustível] + \
                            [self.universo_agua_ar]*7 + [self.universo_vareta_combustível] + [self.universo_agua_ar] + [self.universo_vareta_combustível] + \
                            [self.universo_agua_ar]*7 + [self.universo_vareta_combustível] + [self.universo_agua_ar] + [self.universo_vareta_combustível] + \
                            [self.universo_agua_ar]*7 + [self.universo_vareta_combustível] + [self.universo_agua_ar] + [self.universo_vareta_combustível] + \
                            [self.universo_agua_ar]*7 + [self.universo_vareta_combustível] + [self.universo_agua_ar] + [self.universo_vareta_combustível] + \
                            [self.universo_agua_ar]*7 + [self.universo_vareta_combustível] + [self.universo_agua_ar] + [self.universo_vareta_combustível] + \
                            [self.universo_agua_ar]*3
        anel_comb_mod_9  = [self.universo_vareta_combustível]*54
        anel_comb_mod_8  = [self.universo_vareta_combustível]*48
        anel_comb_mod_7  = [self.universo_vareta_combustível]*42
        anel_comb_mod_6  = [self.universo_vareta_combustível]*36
        anel_comb_mod_5  = [self.universo_vareta_combustível]*30
        anel_comb_mod_4  = [self.universo_vareta_combustível]*24
        anel_comb_mod_3  = [self.universo_vareta_combustível]*18
        anel_comb_mod_2  = [self.universo_vareta_combustível]*12
        anel_comb_mod_1  = [self.universo_vareta_combustível]*6
        if altura_fonte != None:#Se altura for = None, então não existe fonte
            anel_font_mod_0  = [self.universo_fonte]
        else:
            anel_font_mod_0  = [self.universo_agua]

        matriz_hexagonal.universes = [anel_misturado_10, anel_comb_mod_9, anel_comb_mod_8, anel_comb_mod_7, anel_comb_mod_6, anel_comb_mod_5, anel_comb_mod_4, anel_comb_mod_3, anel_comb_mod_2, anel_comb_mod_1, anel_font_mod_0]
        print(matriz_hexagonal)
        #A celula do reator é preenchida com a matriz_hexagonal e depois água, até chegar na superficie lateral do refletor
        self.celula_reator_matriz_hexagonal  = openmc.Cell(fill=matriz_hexagonal, region=-cilindro_raio_externo_grade&-plano_vareta_altura&+plano_fundo_tanque_superior)

        ## Celulas externas a matriz hexagonal

        #Tanque
        self.celula_refletor = openmc.Cell(fill=self.moderador, region=+plano_fundo_tanque_superior&-plano_refletor_lateral_superior&+cilindro_raio_externo_grade&-cilindro_raio_interno_tanque)
        self.celula_ar_externo_matriz = openmc.Cell(fill=self.ar, region=+plano_refletor_lateral_superior&-plano_altura_tanque&+cilindro_raio_externo_grade&-cilindro_raio_interno_tanque)
        self.celula_tanque   = openmc.Cell(fill=self.SS304, region=+cilindro_raio_interno_tanque&-cilindro_raio_externo_tanque&+plano_fundo_tanque_superior&-plano_divisao_altura_tanque_inf
                                                                | +cilindro_raio_interno_tanque&-self.cilindro_beirada_tanque&+plano_divisao_altura_tanque_inf&-plano_divisao_altura_tanque_sup
                                                                | +cilindro_raio_interno_tanque&-cilindro_raio_externo_tanque&+plano_divisao_altura_tanque_sup&-plano_beirada_altura_tanque
                                                                | +plano_fundo_tanque_inferior&-plano_fundo_tanque_superior&-cilindro_raio_externo_tanque
                                                                | +plano_beirada_altura_tanque&-plano_altura_tanque&+cilindro_raio_interno_tanque&-self.cilindro_beirada_tanque)

        self.celula_ar_interna_tanque = openmc.Cell(fill=self.ar, region=-cilindro_raio_interno_tanque&+plano_vareta_altura&-plano_altura_tanque)

        # Boundary conditions
        fronteira = 10
        self.fronteira_ar_lateral = Tanque_Diametro/2 + fronteira
        self.fronteira_ar_superior = altura_tanque + fronteira
        self.fronteira_ar_inferior = fundo_tanque_inferior - fronteira
        cilindro_boundary        = openmc.ZCylinder (r=self.fronteira_ar_lateral, boundary_type='vacuum')
        plano_superior_boundary  = openmc.ZPlane    (z0=self.fronteira_ar_superior, boundary_type='vacuum')
        plano_inferior_boundary  = openmc.ZPlane    (z0=self.fronteira_ar_inferior, boundary_type='vacuum')

        self.celula_ar_externa_tanque = openmc.Cell(fill=self.ar, region=-cilindro_boundary&-plano_superior_boundary&+plano_altura_tanque
                                                            | -cilindro_boundary&+self.cilindro_beirada_tanque&+plano_beirada_altura_tanque&-plano_altura_tanque
                                                            | -cilindro_boundary&+cilindro_raio_externo_tanque&-plano_beirada_altura_tanque&+plano_divisao_altura_tanque_sup
                                                            | -cilindro_boundary&+self.cilindro_beirada_tanque&-plano_divisao_altura_tanque_sup&+plano_divisao_altura_tanque_inf
                                                            | -cilindro_boundary&+cilindro_raio_externo_tanque&-plano_divisao_altura_tanque_inf&+plano_fundo_tanque_inferior
                                                            | -cilindro_boundary&-plano_fundo_tanque_inferior&+plano_inferior_boundary)

        # Universo outer
        self.universo_fora = openmc.Universe(cells=(self.celula_reator_matriz_hexagonal, self.celula_refletor, self.celula_ar_externo_matriz,self.celula_tanque,self.celula_ar_externa_tanque,self.celula_ar_interna_tanque))
        # Célula outer
        self.celula_outer = openmc.Cell(fill=self.universo_fora, region=-cilindro_boundary&+plano_inferior_boundary&-plano_superior_boundary)

        ############ Exportar Geometrias
        self.geometria = openmc.Geometry([self.celula_outer])
        if simu:
            self.geometria.export_to_xml()

    
    def plot2D_secao_transversal(self,basis="xz",width=[142,142],pixels=[15000,15000],origin=(0,0,0)):
        print("################################################")
        print("############        Plot 2D         ############")
        print("################################################")
        ############ Plotar Secão Transversal
        secao_transversal = openmc.Plot.from_geometry(self.geometria)
        secao_transversal.type = 'slice'
        secao_transversal.basis = basis
        secao_transversal.width = width
        secao_transversal.origin = origin
        secao_transversal.filename = 'plot_secao_transversal_' + basis
        secao_transversal.pixels = pixels
        secao_transversal.color_by = 'material'
        secao_transversal.colors = self.colors
        ############ Exportar Plots e Plotar
        plotagem = openmc.Plots((secao_transversal,))
        if simu:
            plotagem.export_to_xml()  
            openmc.plot_geometry()

    def plot3D(self):
        print("################################################")
        print("############        Plot 3D         ############")
        print("################################################")
        ############ Plotar em 3D
        plot_3d = openmc.Plot.from_geometry(self.geometria)
        plot_3d.type = 'voxel'
        plot_3d.filename = 'plot_voxel'
        plot_3d.pixels = (500, 500, 500)
        plot_3d.color_by = 'material'
        plot_3d.colors = self.colors
        plot_3d.width = (150., 150., 150.)
        ############ Exportar Plots e Plotar
        plotagem = openmc.Plots((plot_3d,))
        if simu:
            plotagem.export_to_xml()  
            openmc.plot_geometry()
            openmc.voxel_to_vtk(plot_3d.filename+'.h5', plot_3d.filename)


    def configuracoes(self,fonte=0,particulas=1000,ciclos=100,inativo=10,atrasados=True):
        print("################################################")
        print("########### Definição da Simulação  ############")
        print("################################################")
        # Volume calculation

        self.settings = openmc.Settings()
        self.settings.particles = particulas
        self.settings.batches = ciclos
        self.settings.create_delayed_neutrons = atrasados
        self.settings.photon_transport = True
        if fonte is not None:
            self.atividade = 1.01 * 10**7
            space = openmc.stats.CylindricalIndependent(
                r=openmc.stats.Uniform(0, self.diametro_cilindro_fonte/2),
                phi=openmc.stats.Uniform(0, 2 * np.pi),
                z=openmc.stats.Uniform(self.limite_fonte_inferior, self.limite_fonte_superior),
                #origin=(0.0, 0.0, self.limite_fonte_inferior+9.9)
            )
            #print("r: ", self.diametro_cilindro_fonte/2)
            #print("phi: ", 2 * np.pi)
            #print("z2: ", self.limite_fonte_superior)
            #print("z1: ", self.limite_fonte_inferior)
            #print("origin: ", space.origin)
            angle = openmc.stats.Isotropic()
            energy = openmc.stats.Discrete(
                [
                    3.70E+04,
                    7.80E+04,
                    1.70E+05,
                    3.10E+05,
                    7.00E+05,
                    1.00E+06,
                    2.90E+06,
                    5.60E+06,
                    1.00E+07
                ],
                [
                    1.15E-03,
                    2.89E-03,
                    6.41E-03,
                    1.48E-02,
                    7.05E-02,
                    1.09E-01,
                    1.48E-01,
                    2.69E-01,
                    3.78E-01
                ]
            )
            #time = openmc.stats.Uniform(0, self.atividade)
            self.settings.source = openmc.IndependentSource(space=space,angle=angle,energy=energy)#,time=time)
            self.settings.inactive = 0
            self.settings.run_mode = 'fixed source'
        else:
            self.settings.inactive = inativo
            self.settings.source = openmc.IndependentSource(space=openmc.stats.Point())
        self.settings.output = {'tallies': False}
        if simu:
            self.settings.export_to_xml()
        #Faça ciclos ser acessada de fora
        self.ciclos = ciclos
        print(self.settings)
            
    def run(self):
        print("################################################")
        print("###########        Executando       ############")
        print("################################################")
        if simu:
            openmc.run()

    def tallies(self,):

        energy_filter_thermal  = openmc.EnergyFilter([1.0E-05, 1.0  ]) 
        energy_filter_fast = openmc.EnergyFilter([1.0    , 11E06]) 
        #energy_filter_fonte    = openmc.EnergyFilter([10E06  , 20E06]) 

        neutron_particle_filter = openmc.ParticleFilter(bins='neutron')
        photon_particle_filter = openmc.ParticleFilter(bins='photon')

        print("################################################")
        print("###########          Básicos        ############")
        print("################################################")
        tally_nu = openmc.Tally(name='nu',tally_id=1)
        tally_nu.filters.append(neutron_particle_filter)
        tally_nu.scores.append('nu-fission')
        tally_fission = openmc.Tally(name='reaction rate',tally_id=2)
        tally_fission.filters.append(neutron_particle_filter)
        tally_fission.scores.append('fission')

        print("################################################")
        print("###########       Taxa de dose      ############")
        print("################################################")

        ### TALLY PARA CALCULAR TAXA DE DOSE ###

        # TOPO CENTRAL NEUTRONS
        energy_bins_n, dose_coeffs_n = openmc.data.dose_coefficients(particle='neutron', geometry='ISO')
        energy_function_filter_n = openmc.EnergyFunctionFilter(energy_bins_n, dose_coeffs_n)
        energy_filter_dose_neutron = openmc.EnergyFilter(energy_bins_n)
        
        self.r_divisions_central_leak = [0,3.32/2]
        self.z_divisions_central_leak = [(self.tanque_altura/2),(self.tanque_altura/2)+1]
        mesh_leak_central_filter = openmc.MeshFilter(
            openmc.CylindricalMesh(
                r_grid = (self.r_divisions_central_leak),
                z_grid = (self.z_divisions_central_leak)
                )
            )
        
        dose_tally_neutron_central = openmc.Tally(name='neutron_dose_mesh_leak_central',tally_id=4)
        dose_tally_neutron_central.scores = ['flux']
        dose_tally_neutron_central.filters = [mesh_leak_central_filter , neutron_particle_filter, energy_function_filter_n, energy_filter_dose_neutron]

        # TOPO CENTRAL FÓTONS
        if self.settings.photon_transport:
            energy_bins_p, dose_coeffs_p = openmc.data.dose_coefficients(particle='photon', geometry='ISO')
            energy_function_filter_p = openmc.EnergyFunctionFilter(energy_bins_p, dose_coeffs_p)
            energy_filter_dose_photon  = openmc.EnergyFilter(energy_bins_p)
            dose_tally_photon_central = openmc.Tally(name='photon_dose_mesh_leak_central',tally_id=3)
            dose_tally_photon_central.scores = ['flux']
            dose_tally_photon_central.filters = [mesh_leak_central_filter , photon_particle_filter, energy_function_filter_p, energy_filter_dose_photon]

        # LATERAL ALTURA FONTE NEUTRONS
        r_divisions_lateral_leak = [122/2, (122/2)+1]
        z_divisions_lateral_leak = [self.centro_fonte-10, self.centro_fonte+10.0]
        mesh_leak_lateral        = openmc.CylindricalMesh(r_grid = (r_divisions_lateral_leak), z_grid = (z_divisions_lateral_leak))
        mesh_leak_lateral_filter = openmc.MeshFilter(mesh_leak_lateral)
        dose_tally_neutron_lateral = openmc.Tally(name='neutron_dose_mesh_leak_lateral',tally_id=15) 
        dose_tally_neutron_lateral.filters.append(energy_filter_dose_neutron) 
        dose_tally_neutron_lateral.filters.append(neutron_particle_filter) 
        dose_tally_neutron_lateral.filters.append(energy_function_filter_n)
        dose_tally_neutron_lateral.filters.append(energy_filter_dose_neutron)
        dose_tally_neutron_lateral.filters.append(mesh_leak_lateral_filter)
        dose_tally_neutron_lateral.scores.append('flux')

        if self.settings.photon_transport:
            # LATERAL ALTURA FONTE FÓTONS
            dose_tally_photon_lateral = openmc.Tally(name='photon_dose_mesh_leak_lateral',tally_id=16) 
            dose_tally_photon_lateral.filters.append(energy_filter_dose_photon) 
            dose_tally_photon_lateral.filters.append(photon_particle_filter) 
            dose_tally_photon_lateral.filters.append(energy_function_filter_p)
            dose_tally_photon_lateral.filters.append(energy_filter_dose_photon)
            dose_tally_photon_lateral.filters.append(mesh_leak_lateral_filter)
            dose_tally_photon_lateral.scores.append('flux')

        # TOPO ACIMA DO COMBUSTIVEL NEUTRONS
        r_divisions_top_comb_leak = [0,3.32/2]
        mesh_axial_comb = openmc.CylindricalMesh(origin=(2*2.54,0,0),r_grid = (r_divisions_top_comb_leak), z_grid = (self.z_divisions_central_leak))
        mesh_leak_top_comb_filter = openmc.MeshFilter(mesh_axial_comb)
        dose_tally_neutron_top_comb = openmc.Tally(name='neutron_dose_mesh_leak_top_comb',tally_id=17) 
        dose_tally_neutron_top_comb.filters.append(energy_filter_dose_neutron) 
        dose_tally_neutron_top_comb.filters.append(neutron_particle_filter) 
        dose_tally_neutron_top_comb.filters.append(energy_function_filter_n)
        dose_tally_neutron_top_comb.filters.append(energy_filter_dose_neutron)
        dose_tally_neutron_top_comb.filters.append(mesh_leak_top_comb_filter)
        dose_tally_neutron_top_comb.scores.append('flux')

        if self.settings.photon_transport:
            # TOPO ACIMA DO COMBUSTIVEL FÓTONS
            dose_tally_photon_top_comb = openmc.Tally(name='photon_dose_mesh_leak_top_comb',tally_id=18) 
            dose_tally_photon_top_comb.filters.append(energy_filter_dose_photon) 
            dose_tally_photon_top_comb.filters.append(photon_particle_filter) 
            dose_tally_photon_top_comb.filters.append(energy_function_filter_p)
            dose_tally_photon_top_comb.filters.append(energy_filter_dose_photon)
            dose_tally_photon_top_comb.filters.append(mesh_leak_top_comb_filter)
            dose_tally_photon_top_comb.scores.append('flux')

        # TOPO ACIMA DO REATOR MÉDIA NEUTRONS
        self.r_divisions_avarage_leak = [0,6.5]
        self.z_divisions_avarage_leak = [(self.tanque_altura/2),(self.tanque_altura/2)+1]
        mesh_leak_avarage_filter = openmc.MeshFilter(
            openmc.CylindricalMesh(
                r_grid = (self.r_divisions_avarage_leak),
                z_grid = (self.z_divisions_avarage_leak)
                )
            )
        
        dose_tally_neutron_avarage = openmc.Tally(name='neutron_dose_mesh_leak_avarage',tally_id=80)
        dose_tally_neutron_avarage.scores = ['flux']
        dose_tally_neutron_avarage.filters = [mesh_leak_avarage_filter , neutron_particle_filter, energy_function_filter_n, energy_filter_dose_neutron]

        if self.settings.photon_transport:
            # TOPO ACIMA DO REATOR MÉDIA FÓTONS
            dose_tally_photon_avarage = openmc.Tally(name='photon_dose_mesh_leak_avarage',tally_id=70) 
            dose_tally_photon_avarage.filters.append(energy_filter_dose_photon) 
            dose_tally_photon_avarage.filters.append(photon_particle_filter) 
            dose_tally_photon_avarage.filters.append(energy_function_filter_p)
            dose_tally_photon_avarage.filters.append(energy_filter_dose_photon)
            dose_tally_photon_avarage.filters.append(mesh_leak_avarage_filter)
            dose_tally_photon_avarage.scores.append('flux')

        print("################################################")
        print("###########         Espectro        ############")
        print("################################################")
        
        energy = np.logspace(-5, np.log10(11E06), num=151)
        #energy = [1.0000E-05, 1.0, 5.0E+03, 20.0E+06]
        energy_filter = openmc.EnergyFilter(energy)

        # ESPECTRO COMUM MÉDIO INTERNO AO COMBUSTÍVEL
        tally_spectrum_fuel = openmc.Tally(name='Fluxo espectro interno comb',tally_id=5) # F34
        tally_spectrum_fuel.scores.append('flux')
        tally_spectrum_fuel.filters.append(energy_filter)
        tally_spectrum_fuel.filters.append(neutron_particle_filter)
        tally_spectrum_fuel.filters.append(
            openmc.MeshFilter(
                openmc.CylindricalMesh(
                    r_grid=(0, 1.47/2), #De 0 ao raio interno do combustivel
                    z_grid=(self.centro_fonte-10, self.centro_fonte+10), #No nivel da fonte
                    origin=(2*2.54,0,0)
                    )
                )
            )

        # ESPECTRO COMUM ACIMA FONTE
        
        tally_spectrum_central = openmc.Tally(name='Fluxo espectro acima fonte',tally_id=50) # F34
        tally_spectrum_central.scores.append('flux')
        tally_spectrum_central.filters.append(energy_filter)
        tally_spectrum_central.filters.append(neutron_particle_filter)
        tally_spectrum_central.filters.append(
            openmc.MeshFilter(
                openmc.CylindricalMesh(
                    r_grid=(0, 1.5), #De 0 ao raio interno do tubo central
                    z_grid=(self.centro_fonte+10, self.centro_fonte+10+1),#do topo da fonte a 1cm acima da fonte
                    )
                )
            )


        print("################################################")
        print("###########        Mesh Radial      ############")
        print("################################################")
        divisions_radial = 1000
        r_divisions_radial = np.linspace(0.0,self.fronteira_ar_lateral,divisions_radial+1).tolist()    
        z_divisions_radial = [self.limite_fonte_inferior,self.limite_fonte_superior]
        mesh_radial = openmc.CylindricalMesh(r_grid = (r_divisions_radial), z_grid = (z_divisions_radial))
        mesh_filter_radial = openmc.MeshFilter(mesh_radial)

        tally_radial_thermal = openmc.Tally(name='MESH_Radial_Termico',tally_id=6)
        tally_radial_thermal.filters.append(mesh_filter_radial)
        tally_radial_thermal.filters.append(energy_filter_thermal)
        tally_radial_thermal.filters.append(neutron_particle_filter)
        tally_radial_thermal.scores.append('flux')

        tally_radial_fast = openmc.Tally(name='MESH_Radial_Rapido',tally_id=7)
        tally_radial_fast.filters.append(mesh_filter_radial)
        tally_radial_fast.filters.append(energy_filter_fast)
        tally_radial_fast.filters.append(neutron_particle_filter)
        tally_radial_fast.scores.append('flux')

        #tally_radial_fonte = openmc.Tally(name='MESH_Radial_Fonte',tally_id=8)
        #tally_radial_fonte.filters.append(mesh_filter_radial)
        #tally_radial_fonte.filters.append(energy_filter_fonte)
        #tally_radial_fonte.scores.append('flux')

        print("################################################")
        print("###########        Mesh Axial       ############")
        print("################################################")
        divisions_axial = 1000
        z_divisions_axial = np.linspace(self.fronteira_ar_inferior,self.fronteira_ar_superior,divisions_axial+1).tolist()    
        r_divisions_axial_central = [0,3.32/2]
        mesh_axial_central = openmc.CylindricalMesh(r_grid = (r_divisions_axial_central), z_grid = (z_divisions_axial))
        mesh_filter_axial_central = openmc.MeshFilter(mesh_axial_central)

        tally_axial_central_thermal = openmc.Tally(name='MESH_Axial_Central_Thermal',tally_id=9)
        tally_axial_central_thermal.filters.append(energy_filter_thermal)
        tally_axial_central_thermal.filters.append(mesh_filter_axial_central)
        tally_axial_central_thermal.filters.append(neutron_particle_filter)
        tally_axial_central_thermal.scores.append('flux')

        r_divisions_axial_comb = [0,1.47/2]
        mesh_axial_comb = openmc.CylindricalMesh(origin=(2*2.54,0,0),r_grid = (r_divisions_axial_comb), z_grid = (z_divisions_axial))
        mesh_filter_axial_comb = openmc.MeshFilter(mesh_axial_comb)

        tally_axial_comb_thermal = openmc.Tally(name='MESH_Axial_Comb_Thermal',tally_id=10)
        tally_axial_comb_thermal.filters.append(energy_filter_thermal)
        tally_axial_comb_thermal.filters.append(mesh_filter_axial_comb)
        tally_axial_comb_thermal.filters.append(neutron_particle_filter)
        tally_axial_comb_thermal.scores.append('flux')

        tally_axial_central_fast = openmc.Tally(name='MESH_Axial_Central_Fast',tally_id=11)
        tally_axial_central_fast.filters.append(energy_filter_fast)
        tally_axial_central_fast.filters.append(mesh_filter_axial_central)
        tally_axial_central_fast.filters.append(neutron_particle_filter)
        tally_axial_central_fast.scores.append('flux')

        tally_axial_comb_fast = openmc.Tally(name='MESH_Axial_Comb_Fast',tally_id=12)
        tally_axial_comb_fast.filters.append(energy_filter_fast)
        tally_axial_comb_fast.filters.append(mesh_filter_axial_comb)
        tally_axial_comb_fast.filters.append(neutron_particle_filter)
        tally_axial_comb_fast.scores.append('flux')

        print("################################################")
        print("###########        Mesh Cubico      ############")
        print("################################################")
        divisions_cubico = 284
        x_divisions = np.linspace(-self.fronteira_ar_lateral,self.fronteira_ar_lateral,divisions_cubico+1).tolist()
        y_divisions = np.linspace(-self.fronteira_ar_lateral,self.fronteira_ar_lateral,divisions_cubico+1).tolist()  
        z_divisions = [self.limite_fonte_inferior,self.limite_fonte_superior] #np.linspace( self.fronteira_ar_inferior,self.fronteira_ar_superior,101).tolist() 
        mesh_cubico = openmc.RectilinearMesh()
        mesh_cubico.x_grid = x_divisions
        mesh_cubico.y_grid = y_divisions
        mesh_cubico.z_grid = z_divisions
        mesh_filter_cubico = openmc.MeshFilter(mesh_cubico)
        # TERMICO
        tally_cubico_termico = openmc.Tally(name='MESH_Cubico_Termico',tally_id=13)
        tally_cubico_termico.filters.append(mesh_filter_cubico)
        tally_cubico_termico.filters.append(energy_filter_thermal)
        tally_cubico_termico.filters.append(neutron_particle_filter)
        tally_cubico_termico.scores.append('flux')
        # RAPIDO
        tally_cubico_rapido = openmc.Tally(name='MESH_Cubico_Rapido',tally_id=19)
        tally_cubico_rapido.filters.append(mesh_filter_cubico)
        tally_cubico_rapido.filters.append(energy_filter_fast)
        tally_cubico_rapido.filters.append(neutron_particle_filter)
        tally_cubico_rapido.scores.append('flux')

        ###############################################################################################################
        divisions_rc = 1000
        x_divisions_rc = np.linspace(0,self.fronteira_ar_lateral,divisions_rc+1).tolist()
        y_divisions_rc = [-1.47/2, 1.47/2] #diametro interno do anel do combustivel
        z_divisions_rc = [self.limite_fonte_inferior,self.limite_fonte_superior] 
        mesh_cubico_rc = openmc.RectilinearMesh()
        mesh_cubico_rc.x_grid = x_divisions_rc
        mesh_cubico_rc.y_grid = y_divisions_rc
        mesh_cubico_rc.z_grid = z_divisions_rc
        mesh_filter_cubico_rc = openmc.MeshFilter(mesh_cubico_rc)
        # TERMICO
        tally_cubico_termico_rc = openmc.Tally(name='MESH_Cubico_Termico_rc',tally_id=20)
        tally_cubico_termico_rc.filters.append(mesh_filter_cubico_rc)
        tally_cubico_termico_rc.filters.append(energy_filter_thermal)
        tally_cubico_termico_rc.filters.append(neutron_particle_filter)
        tally_cubico_termico_rc.scores.append('flux')
        # RAPIDO
        tally_cubico_rapido_rc = openmc.Tally(name='MESH_Cubico_Rapido_rc',tally_id=21)
        tally_cubico_rapido_rc.filters.append(mesh_filter_cubico_rc)
        tally_cubico_rapido_rc.filters.append(energy_filter_fast)
        tally_cubico_rapido_rc.filters.append(neutron_particle_filter)
        tally_cubico_rapido_rc.scores.append('flux')

        ############## Coleção de tallies ##############
        vetor_tallies = []
        vetor_tallies.append(tally_axial_central_thermal)
        vetor_tallies.append(tally_axial_central_fast)
        vetor_tallies.append(tally_axial_comb_thermal)
        vetor_tallies.append(tally_axial_comb_fast)
        vetor_tallies.append(tally_cubico_termico)
        vetor_tallies.append(tally_cubico_rapido)
        vetor_tallies.append(tally_cubico_termico_rc)
        vetor_tallies.append(tally_cubico_rapido_rc)
        vetor_tallies.append(tally_spectrum_central)
        vetor_tallies.append(tally_spectrum_fuel)
        vetor_tallies.append(tally_radial_thermal)
        vetor_tallies.append(tally_radial_fast)
        #vetor_tallies.append(tally_radial_fonte)
        vetor_tallies.append(tally_fission)
        vetor_tallies.append(tally_nu)
        vetor_tallies.append(dose_tally_neutron_central)
        vetor_tallies.append(dose_tally_neutron_lateral)
        vetor_tallies.append(dose_tally_neutron_top_comb)
        vetor_tallies.append(dose_tally_neutron_avarage)
        if self.settings.photon_transport:
            vetor_tallies.append(dose_tally_photon_central)
            vetor_tallies.append(dose_tally_photon_lateral)
            vetor_tallies.append(dose_tally_photon_top_comb)
            vetor_tallies.append(dose_tally_photon_avarage)
        tallies = openmc.Tallies(vetor_tallies)
        if simu:
            tallies.export_to_xml()
            self.run()
        if self.fonte:
            print("################################################")
            print("########  Trabalhando dados com fonte  #########")
            print("################################################")
            sp = openmc.StatePoint('statepoint.'+str(self.ciclos)+'.h5')

            uncertainty = 0.05

            # Acesse os resultados do tally radial
            nu                   = sp.get_tally(scores=['nu-fission'])
            fission              = sp.get_tally(scores=['fission'])
            flux_radial_thermal  = sp.get_tally(scores=['flux'], name='MESH_Radial_Termico')
            flux_radial_fast     = sp.get_tally(scores=['flux'], name='MESH_Radial_Rapido')
            #flux_radial_fonte     = sp.get_tally(scores=['flux'], name='MESH_Radial_Fonte')

            flux_espectro_fuel   = sp.get_tally(scores=['flux'], name='Fluxo espectro interno comb')
            flux_espectro_central   = sp.get_tally(scores=['flux'], name='Fluxo espectro acima fonte')

            flux_cubico_thermal        = sp.get_tally(scores=['flux'], name='MESH_Cubico_Termico')
            flux_cubico_fast           = sp.get_tally(scores=['flux'], name='MESH_Cubico_Rapido')
            flux_cubico_thermal_rc     = sp.get_tally(scores=['flux'], name='MESH_Cubico_Termico_rc')
            flux_cubico_fast_rc        = sp.get_tally(scores=['flux'], name='MESH_Cubico_Rapido_rc')

            flux_axial_comb_thermal    = sp.get_tally(scores=['flux'], name='MESH_Axial_Comb_Thermal')
            flux_axial_central_thermal = sp.get_tally(scores=['flux'], name='MESH_Axial_Central_Thermal')
            flux_axial_comb_fast       = sp.get_tally(scores=['flux'], name='MESH_Axial_Comb_Fast')
            flux_axial_central_fast    = sp.get_tally(scores=['flux'], name='MESH_Axial_Central_Fast')

            dose_n_central             = sp.get_tally(scores=['flux'], name='neutron_dose_mesh_leak_central')
            dose_n_lateral             = sp.get_tally(scores=['flux'], name='neutron_dose_mesh_leak_lateral')
            dose_n_top_comb            = sp.get_tally(scores=['flux'], name='neutron_dose_mesh_leak_top_comb')
            dose_n_avarage             = sp.get_tally(scores=['flux'], name='neutron_dose_mesh_leak_avarage')
            if self.settings.photon_transport:
                dose_p_central         = sp.get_tally(scores=['flux'], name='photon_dose_mesh_leak_central')
                dose_p_lateral         = sp.get_tally(scores=['flux'], name='photon_dose_mesh_leak_lateral')
                dose_p_top_comb        = sp.get_tally(scores=['flux'], name='photon_dose_mesh_leak_top_comb')
                dose_p_avarage         = sp.get_tally(scores=['flux'], name='photon_dose_mesh_leak_avarage')

            # Tallies básicos
            nu_mean                 = nu.mean
            nu_std_dev              = nu.std_dev
            fission_mean            = fission.mean
            fission_std_dev         = fission.std_dev
            # Espectro
            flux_espectro_fuel_mean    = flux_espectro_fuel.mean
            flux_espectro_fuel_dev     = flux_espectro_fuel.std_dev
            flux_espectro_central_mean = flux_espectro_central.mean
            flux_espectro_central_dev  = flux_espectro_central.std_dev
            # Mesh Radial
            #flux_rad_fonte_mean     = flux_radial_fonte.mean
            #flux_rad_fonte_dev      = flux_radial_fonte.std_dev
            flux_rad_fast_mean      = flux_radial_fast.mean
            flux_rad_fast_dev       = flux_radial_fast.std_dev
            flux_rad_thermal_mean   = flux_radial_thermal.mean
            flux_rad_thermal_dev    = flux_radial_thermal.std_dev
            # Mesh Cúbico
            flux_cub_thermal        = flux_cubico_thermal.mean
            flux_cub_thermal_dev    = flux_cubico_thermal.std_dev
            flux_cub_fast           = flux_cubico_fast.mean
            flux_cub_fast_dev       = flux_cubico_fast.std_dev
            flux_cub_thermal_rc     = flux_cubico_thermal_rc.mean
            flux_cub_thermal_dev_rc = flux_cubico_thermal_rc.std_dev
            flux_cub_fast_rc        = flux_cubico_fast_rc.mean
            flux_cub_fast_dev_rc    = flux_cubico_fast_rc.std_dev
            # Mesh Axial
            flux_axial_comb_thermal_mean    = flux_axial_comb_thermal.mean
            flux_axial_comb_thermal_dev     = flux_axial_comb_thermal.std_dev
            flux_axial_central_thermal_mean = flux_axial_central_thermal.mean
            flux_axial_central_thermal_dev  = flux_axial_central_thermal.std_dev
            flux_axial_comb_fast_mean       = flux_axial_comb_fast.mean
            flux_axial_comb_fast_dev        = flux_axial_comb_fast.std_dev
            flux_axial_central_fast_mean    = flux_axial_central_fast.mean
            flux_axial_central_fast_dev     = flux_axial_central_fast.std_dev
            # Taxas de dose
            dose_n_central_mean             = dose_n_central.mean
            dose_n_central_dev              = dose_n_central.std_dev
            dose_n_lateral_mean             = dose_n_lateral.mean
            dose_n_lateral_dev              = dose_n_lateral.std_dev
            dose_n_top_comb_mean            = dose_n_top_comb.mean
            dose_n_top_comb_dev             = dose_n_top_comb.std_dev
            dose_n_avarage_mean             = dose_n_avarage.mean
            dose_n_avarage_dev              = dose_n_avarage.std_dev
            if self.settings.photon_transport:
                dose_p_central_mean             = dose_p_central.mean
                dose_p_central_dev              = dose_p_central.std_dev
                dose_p_lateral_mean             = dose_p_lateral.mean
                dose_p_lateral_dev              = dose_p_lateral.std_dev
                dose_p_top_comb_mean            = dose_p_top_comb.mean
                dose_p_top_comb_dev             = dose_p_top_comb.std_dev
                dose_p_avarage_mean             = dose_p_avarage.mean
                dose_p_avarage_dev              = dose_p_avarage.std_dev

            ### TAXAS DE DOSE ###
            # Volume de fuga central no topo do reator
            volume_dose_central = np.pi * (self.z_divisions_central_leak[1]-self.z_divisions_central_leak[0]) * (self.r_divisions_central_leak[1]**2)
            # Volume de fuga lateral 
            r2 = (122/2)+1
            r1 = 122/2
            h = 20
            volume_dose_lateral = h*(np.pi*r2**2-np.pi*r1**2)
            # Volume de fuga acima do combustível mais perto do centro
            volume_dose_top_comb = np.pi * (self.z_divisions_central_leak[1]-self.z_divisions_central_leak[0]) * (self.r_divisions_central_leak[1]**2)
            # Volume de fuga mpedio do topo
            volume_dose_top_avarage = np.pi * (self.z_divisions_central_leak[1]-self.z_divisions_central_leak[0]) * (self.r_divisions_avarage_leak[1]**2)

            dose_central_n = []
            dose_central_n_dev = []
            dose_central_n_energy = []

            # NÊUTRONS
            for i in range(0,len(dose_n_central_mean)):
                dose=self.atividade*dose_n_central_mean[i][0][0]/volume_dose_central
                incerteza=self.atividade*dose_n_central_dev[i][0][0]/volume_dose_central
                if incerteza/dose < uncertainty:
                    dose_central_n.append(dose)
                    dose_central_n_dev.append(incerteza)
                    dose_central_n_energy.append(energy_bins_n[i])
            print("dose_central_n", np.sum(dose_central_n))

            dose_lateral_n = []
            dose_lateral_n_dev = []
            dose_lateral_n_energy = []

            for i in range(0,len(dose_n_lateral_mean)):
                dose=self.atividade*dose_n_lateral_mean[i][0][0]/volume_dose_lateral
                incerteza=self.atividade*dose_n_lateral_dev[i][0][0]/volume_dose_lateral
                if incerteza/dose < uncertainty:
                    dose_lateral_n.append(dose)
                    dose_lateral_n_dev.append(incerteza)
                    dose_lateral_n_energy.append(energy_bins_n[i])
            print("dose_lateral_n", np.sum(dose_lateral_n))

            dose_top_comb_n = []
            dose_top_comb_n_dev = []
            dose_top_comb_n_energy = []

            for i in range(0,len(dose_n_top_comb_mean)):
                dose=self.atividade*dose_n_top_comb_mean[i][0][0]/volume_dose_top_comb
                incerteza=self.atividade*dose_n_top_comb_dev[i][0][0]/volume_dose_top_comb
                if incerteza/dose < uncertainty:
                    dose_top_comb_n.append(dose)
                    dose_top_comb_n_dev.append(incerteza)
                    dose_top_comb_n_energy.append(energy_bins_n[i])
            print("dose_top_comb_n", np.sum(dose_top_comb_n))

            dose_avarage_n = []
            dose_avarage_n_dev = []
            dose_avarage_n_energy = []

            for i in range(0,len(dose_n_avarage_mean)):
                dose=self.atividade*dose_n_avarage_mean[i][0][0]/volume_dose_top_avarage
                incerteza=self.atividade*dose_n_avarage_dev[i][0][0]/volume_dose_top_avarage
                if incerteza/dose < uncertainty:
                    dose_avarage_n.append(dose)
                    dose_avarage_n_dev.append(incerteza)
                    dose_avarage_n_energy.append(energy_bins_n[i])
            print("dose_avarage_n", np.sum(dose_avarage_n))

            # FÓTONS
            if self.settings.photon_transport:

                dose_central_p = []
                dose_central_p_dev = []
                dose_central_p_energy = []

                for i in range(0,len(dose_p_central_mean)):
                    dose=self.atividade*dose_p_central_mean[i][0][0]/volume_dose_central
                    incerteza=self.atividade*dose_p_central_dev[i][0][0]/volume_dose_central
                    if incerteza/dose < uncertainty:
                        dose_central_p.append(dose)
                        dose_central_p_dev.append(incerteza)
                        dose_central_p_energy.append(energy_bins_p[i])
                print("dose_central_p", np.sum(dose_central_p))

                dose_lateral_p = []
                dose_lateral_p_dev = []
                dose_lateral_p_energy = []

                for i in range(0,len(dose_p_lateral_mean)):
                    dose=self.atividade*dose_p_lateral_mean[i][0][0]/volume_dose_lateral
                    incerteza=self.atividade*dose_p_lateral_dev[i][0][0]/volume_dose_lateral
                    if incerteza/dose < uncertainty:
                        dose_lateral_p.append(dose)
                        dose_lateral_p_dev.append(incerteza)
                        dose_lateral_p_energy.append(energy_bins_p[i])
                print("dose_lateral_p", np.sum(dose_lateral_p))
                
                dose_top_comb_p = []
                dose_top_comb_p_dev = []
                dose_top_comb_p_energy = []

                for i in range(0,len(dose_p_top_comb_mean)):
                    dose=self.atividade*dose_p_top_comb_mean[i][0][0]/volume_dose_top_comb
                    incerteza=self.atividade*dose_p_top_comb_dev[i][0][0]/volume_dose_top_comb
                    if incerteza/dose < uncertainty:
                        dose_top_comb_p.append(dose)
                        dose_top_comb_p_dev.append(incerteza)
                        dose_top_comb_p_energy.append(energy_bins_p[i])
                print("dose_top_comb_p", np.sum(dose_top_comb_p))

                dose_avarage_p = []
                dose_avarage_p_dev = []
                dose_avarage_p_energy = []

                for i in range(0,len(dose_p_avarage_mean)):
                    dose=self.atividade*dose_p_avarage_mean[i][0][0]/volume_dose_top_avarage
                    incerteza=self.atividade*dose_p_avarage_dev[i][0][0]/volume_dose_top_avarage
                    if incerteza/dose < uncertainty:
                        dose_avarage_p.append(dose)
                        dose_avarage_p_dev.append(incerteza)
                        dose_avarage_p_energy.append(energy_bins_p[i])
                print("dose_avarage_p", np.sum(dose_avarage_p))

            ### MESH AXIAL ###
            # Volumes das areas do mesh axial central
            volume_axial_central = []
            for i in range(0, divisions_axial):  # Use o número apropriado de intervalos
                z1 = z_divisions_axial[i]
                z2 = z_divisions_axial[i + 1]
                r = r_divisions_axial_central[1] - r_divisions_axial_central[0]
                volume_axial_central.append(np.pi * (z2 - z1) * r**2)

            # Volumes das areas do mesh axial comb
            volume_axial_comb = []
            for i in range(0, divisions_axial):  # Use o número apropriado de intervalos
                z1 = z_divisions_axial[i]
                z2 = z_divisions_axial[i + 1]
                r = r_divisions_axial_comb[1] - r_divisions_axial_comb[0]
                volume_axial_comb.append(np.pi * (z2 - z1) * r**2)

            fluxo_axial_central_thermal = []
            fluxo_axial_central_thermal_dev = []
            fluxo_z_central_thermal = []

            print()
            print(' MESH axial central thermal:')
            for i in range(0,divisions_axial):
                fluxo=self.atividade*flux_axial_central_thermal_mean[i][0][0]/volume_axial_central[i]
                incerteza=self.atividade*flux_axial_central_thermal_dev[i][0][0]/volume_axial_central[i]
                if incerteza/fluxo < uncertainty:
                    if z_divisions_axial[i]<self.centro_fonte+10 and z_divisions_axial[i]>self.centro_fonte-10:
                        fluxo_axial_central_thermal.append(fluxo)
                        fluxo_axial_central_thermal_dev.append(incerteza)
                        fluxo_z_central_thermal.append(z_divisions_axial[i])

            print(self.centro_fonte+10 , "\t", self.centro_fonte-10)
            print(self.centro_fonte)
            fluxo_axial_comb_thermal = []
            fluxo_axial_comb_thermal_dev = []
            fluxo_z_comb_thermal = []

            print()
            print(' MESH axial comb thermal:')
            for i in range(0,divisions_axial):
                fluxo=self.atividade*flux_axial_comb_thermal_mean[i][0][0]/volume_axial_comb[i]
                incerteza=self.atividade*flux_axial_comb_thermal_dev[i][0][0]/volume_axial_comb[i]
                if incerteza/fluxo < uncertainty:
                    fluxo_axial_comb_thermal.append(fluxo)
                    fluxo_axial_comb_thermal_dev.append(incerteza)
                    fluxo_z_comb_thermal.append(z_divisions_axial[i])
            
            fluxo_axial_central_fast = []
            fluxo_axial_central_fast_dev = []
            fluxo_z_central_fast = []

            print()
            print(' MESH axial central fast:')
            for i in range(0,divisions_axial):
                fluxo=self.atividade*flux_axial_central_fast_mean[i][0][0]/volume_axial_central[i]
                incerteza=self.atividade*flux_axial_central_fast_dev[i][0][0]/volume_axial_central[i]
                if incerteza/fluxo < uncertainty:
                    if z_divisions_axial[i]<self.centro_fonte+10 and z_divisions_axial[i]>self.centro_fonte-10:
                        fluxo_axial_central_fast.append(fluxo)
                        fluxo_axial_central_fast_dev.append(incerteza)
                        fluxo_z_central_fast.append(z_divisions_axial[i])

            fluxo_axial_comb_fast = []
            fluxo_axial_comb_fast_dev = []
            fluxo_z_comb_fast = []

            print()
            print(' MESH axial comb fast:')
            for i in range(0,divisions_axial):
                fluxo=self.atividade*flux_axial_comb_fast_mean[i][0][0]/volume_axial_comb[i]
                incerteza=self.atividade*flux_axial_comb_fast_dev[i][0][0]/volume_axial_comb[i]
                if incerteza/fluxo < uncertainty:
                    fluxo_axial_comb_fast.append(fluxo)
                    fluxo_axial_comb_fast_dev.append(incerteza)
                    fluxo_z_comb_fast.append(z_divisions_axial[i])

            ### MESH CÚBICO ###
            print()
            print(' MESH cúbico termico:')
            V = (x_divisions[1]-x_divisions[0])*(y_divisions[1]-y_divisions[0])*(z_divisions[1]-z_divisions[0])
            amplitude_thermal = self.atividade*flux_cub_thermal.reshape((divisions_cubico, divisions_cubico))/V
            amplitude_thermal_dev = self.atividade*flux_cub_thermal_dev.reshape((divisions_cubico, divisions_cubico))/V
            amplitude_thermal_dev_per = np.zeros((len(amplitude_thermal_dev),len(amplitude_thermal_dev)))
            for i in range(0,divisions_cubico-1):
                for j in range(0,divisions_cubico-1):
                    if amplitude_thermal[i][j] != 0:
                        amplitude_thermal_dev_per[i][j] = amplitude_thermal_dev[i][j] / amplitude_thermal[i][j]
                        if amplitude_thermal_dev_per[i][j] > uncertainty:
                            amplitude_thermal[i][j] = None
                    else:
                        amplitude_thermal_dev_per[i][j] = 1
                        amplitude_thermal[i][j] = None

            print()
            print(' MESH cúbico rapido:')
            V = (x_divisions[1]-x_divisions[0])*(y_divisions[1]-y_divisions[0])*(z_divisions[1]-z_divisions[0])
            amplitude_fast = self.atividade*flux_cub_fast.reshape((divisions_cubico, divisions_cubico))/V
            amplitude_fast_dev = self.atividade*flux_cub_fast_dev.reshape((divisions_cubico, divisions_cubico))/V
            amplitude_fast_dev_per = np.zeros((len(amplitude_fast_dev),len(amplitude_fast_dev)))
            for i in range(0,divisions_cubico-1):
                for j in range(0,divisions_cubico-1):
                    if amplitude_fast[i][j] != 0:
                        amplitude_fast_dev_per[i][j] = amplitude_fast_dev[i][j] / amplitude_fast[i][j]
                        if amplitude_fast_dev_per[i][j] > uncertainty:
                            amplitude_fast[i][j] = None
                    else:
                        amplitude_fast_dev_per[i][j] = 1
                        amplitude_fast[i][j] = None

            ### ESPECTRO DO COMBUSTÍVEL ###
            print()
            print(' Espectro de fluxo:')
            print()

            flux_spec_fuel_mean      = []
            flux_spec_fuel_dev       = []
            flux_spec_fuel_energy    = []



            V=20*np.pi*(1.47/2)**2

            for i in range(0,len(energy)-1):
                fluxo = self.atividade*flux_espectro_fuel_mean[i][0][0]/V
                incerteza = self.atividade*flux_espectro_fuel_dev[i][0][0]/V
                if incerteza/fluxo < uncertainty:
                    flux_spec_fuel_mean.append(fluxo)
                    flux_spec_fuel_dev.append(incerteza)
                    flux_spec_fuel_energy.append(energy[i+1])
                    #print(" Intervalo ", i,": ","\t", format(flux_espectro[-1], '.4e'), "+/-", format(flux_dev_espectro[-1], '.4e'), "[neutron/cm².s]")

            V=1*np.pi*1.5**2
            flux_spec_central_mean   = []
            flux_spec_central_dev    = []
            flux_spec_central_energy = []
            for i in range(0,len(energy)-1):
                fluxo = self.atividade*flux_espectro_central_mean[i][0][0]/V
                incerteza = self.atividade*flux_espectro_central_dev[i][0][0]/V
                if incerteza/fluxo < uncertainty:
                    flux_spec_central_mean.append(fluxo)
                    flux_spec_central_dev.append(incerteza)
                    flux_spec_central_energy.append(energy[i+1])
                    #print(" Intervalo ", i,": ","\t", format(flux_espectro[-1], '.4e'), "+/-", format(flux_dev_espectro[-1], '.4e'), "[neutron/cm².s]")


            ### MESH RADIAL ###
            # Retirando o mesh radial
            print('')
            print("MESH RADIAL:")
            print('')

            # Volumes das areas do mesh radial
            volume_radial = []
            for i in range(0, divisions_radial):  # Use o número apropriado de intervalos
                r1 = r_divisions_radial[i]
                r2 = r_divisions_radial[i + 1]
                h = z_divisions_radial[1] - z_divisions_radial[0]
                volume_radial.append(3.14159265359 * (r2**2 - r1**2) * h)

            flux_rad_termico = []  # Vetores para armazenar resultados
            flux_dev_termico = []
            flux_r_termico   = []

            #flux_rad_fonte   = []
            #flux_dev_fonte   = []
            #flux_r_fonte     = []

            flux_rad_rapido  = []  
            flux_dev_rapido  = []
            flux_r_rapido    = []

            print()
            print(' Fluxo em intervalo térmico:')
            print()
            # Fluxo em intervalo térmico 
            for i in range(0,divisions_radial):
                fluxo=self.atividade*flux_rad_thermal_mean[i][0][0]/volume_radial[i]
                incerteza=self.atividade*flux_rad_thermal_dev[i][0][0]/volume_radial[i]
                if incerteza/fluxo < uncertainty:
                    if r_divisions_radial[i]>3.52/2:
                        flux_rad_termico.append(fluxo)
                        flux_dev_termico.append(incerteza)
                        flux_r_termico.append(r_divisions_radial[i])

                #print(" Intervalo ", i,": ","\t", format(flux_rad_termico[i], '.4e'), "+/-", format(flux_dev_termico[i], '.4e'), "[neutron/cm².s]")

            #print()
            #print(' Fluxo em intervalo de ressonancia:')
            #print()
            ## Fluxo em intervalo de ressonancia 
            #for i in range(0,divisions_radial-1):
            #    fluxo=self.atividade*flux_rad_fonte_mean[i][0][0]/volume_radial[i]
            #    incerteza=self.atividade*flux_rad_fonte_dev[i][0][0]/volume_radial[i]
            #    if incerteza/fluxo < 0.05:
            #        flux_rad_fonte.append(fluxo)
            #        flux_dev_fonte.append(incerteza)
            #        flux_r_fonte.append(r_divisions_radial[i])
            #
            #    #print(" Intervalo ", i,": ","\t", format(flux_rad_fonte[i], '.4e'), "+/-", format(flux_dev_fonte[i], '.4e'), "[neutron/cm².s]")

            print()
            print(' Fluxo em intervalo rapido:')
            print()
            # Fluxo em intervalo térmico 
            for i in range(0,divisions_radial):
                fluxo=self.atividade*flux_rad_fast_mean[i][0][0]/volume_radial[i]
                incerteza=self.atividade*flux_rad_fast_dev[i][0][0]/volume_radial[i]
                if incerteza/fluxo < uncertainty:
                    if r_divisions_radial[i]>3.52/2:
                        flux_rad_rapido.append(fluxo)
                        flux_dev_rapido.append(incerteza)
                        flux_r_rapido.append(r_divisions_radial[i])
                #print(" Intervalo ", i,": ","\t", format(flux_rad_rapido[i], '.4e'), "+/-", format(flux_dev_rapido[i], '.4e'), "[neutron/cm².s]")



            ###################################
            #Radial cubico
            ###################################

            flux_cub_thermal_rc
            flux_cub_thermal_dev_rc
            flux_cub_fast_rc
            flux_cub_fast_dev_rc


            flux_rc_rad_termico = []  # Vetores para armazenar resultados
            flux_rc_dev_termico = []
            flux_rc_r_termico   = []

            flux_rc_rad_rapido  = []  
            flux_rc_dev_rapido  = []
            flux_rc_r_rapido    = []

            volume_cubico_radial=(
                (x_divisions_rc[1]-x_divisions_rc[0])*
                (y_divisions_rc[1]-y_divisions_rc[0])*
                (z_divisions_rc[1]-z_divisions_rc[0])
                )

            print()
            print('RC Fluxo em intervalo térmico:')
            print()
            # Fluxo em intervalo térmico 
            for i in range(0,divisions_rc):
                fluxo=self.atividade*flux_cub_thermal_rc[i][0][0]/volume_cubico_radial
                incerteza=self.atividade*flux_cub_thermal_dev_rc[i][0][0]/volume_cubico_radial
                if incerteza/fluxo < uncertainty:
                    if x_divisions_rc[i]>self.diametro_cilindro_ar_fonte/2:
                        flux_rc_rad_termico.append(fluxo)
                        flux_rc_dev_termico.append(incerteza)
                        flux_rc_r_termico.append(x_divisions_rc[i])

            print()
            print('RC Fluxo em intervalo rápido:')
            print()
            # Fluxo em intervalo térmico 
            for i in range(0,divisions_rc):
                fluxo=self.atividade*flux_cub_fast_rc[i][0][0]/volume_cubico_radial
                incerteza=self.atividade*flux_cub_fast_dev_rc[i][0][0]/volume_cubico_radial
                if incerteza/fluxo < uncertainty:
                    if x_divisions_rc[i]>self.diametro_cilindro_ar_fonte/2:
                        flux_rc_rad_rapido.append(fluxo)
                        flux_rc_dev_rapido.append(incerteza)
                        flux_rc_r_rapido.append(x_divisions_rc[i])



            # PLT styles

            # ['Solarize_Light2', '_classic_test_patch', '_mpl-gallery', '_mpl-gallery-nogrid', 'bmh', 'classic', 'dark_background', 'fast', 'fivethirtyeight', 
            # 'ggplot', 'grayscale', 'seaborn-v0_8', 'seaborn-v0_8-bright', 'seaborn-v0_8-colorblind', 'seaborn-v0_8-dark', 'seaborn-v0_8-dark-palette', 
            # 'seaborn-v0_8-darkgrid', 'seaborn-v0_8-deep', 'seaborn-v0_8-muted', 'seaborn-v0_8-notebook', 'seaborn-v0_8-paper', 'seaborn-v0_8-pastel', 
            # 'seaborn-v0_8-poster', 'seaborn-v0_8-talk', 'seaborn-v0_8-ticks', 'seaborn-v0_8-white', 'seaborn-v0_8-whitegrid', 'tableau-colorblind10']

            # Colors

            # xkcd color survey, prefixed with 'xkcd:' (e.g., 'xkcd:sky blue'; case insensitive) https://xkcd.com/color/rgb/
            # Tableau Colors from the 'T10' categorical palette:
            # {'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan'}

            ################################################################################################################################
            #Seção de choque
            #U235 = openmc.data.IncidentNeutron.from_hdf5('/opt/nuclear-data/endfb-viii.0-hdf5/neutron/U235.h5')
            #pprint(list(U235.reactions.values()))
            #pprint(U235.energy)
            #U235_f = U235[18]
            #plt.plot(U235.energy['294K']), U235_f.xs['294K'](U235.energy['294K'])

            ################################################################################################################################
            # MESH CÚBICO THERMAL - GRAFICO 3D
            # Função para escurecer e aumentar a vivacidade do colormap
            def enhance_and_darken_cmap(cmap, gamma, scale):
                colors = cmap(np.linspace(0, 1, 256)) ** gamma
                darkened_colors = colors * scale
                darkened_colors[:, -1] = 1  # Manter a opacidade
                return LinearSegmentedColormap.from_list('enhanced_darkened_coolwarm', darkened_colors)

            # Criar um colormap com cores mais vivas e escurecidas
            coolwarm = plt.get_cmap('coolwarm')
            enhanced_darkened_coolwarm = enhance_and_darken_cmap(coolwarm, gamma=1.0, scale=1.0)

            plt.style.use('seaborn-v0_8-paper')
            from matplotlib.ticker import LogFormatter

            # Criando a figura e os eixos 3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            #Formatação cientifica
            def scientific_format(x, pos):
                return f'{x:.1e}'.replace('e+0', 'e').replace('e+','e')

            # Plotando a superfície 3D
            X, Y = np.meshgrid(x_divisions[1:], y_divisions[1:])
            surf = ax.plot_surface(X, Y, amplitude_thermal, cmap=enhanced_darkened_coolwarm, edgecolor='black')
            cbar = fig.colorbar(surf, aspect=10.0, fraction=0.02, pad=0.018)
            cbar.set_label('Flux (neutrons.cm⁻².s⁻¹)', size=20, labelpad=8)
            cbar.ax.tick_params(axis='y', labelsize=12)
            cbar.ax.yaxis.set_major_formatter(FuncFormatter(scientific_format))
            # Adicionando rótulos aos eixos
            ax.set_xlabel('X position (cm)', size=20, labelpad=6)
            ax.set_ylabel('Y position (cm)', size=20, labelpad=6)
            #ax.set_zlabel('Flux (neutrons.cm⁻².s⁻¹)',size=20,labelpad=15)
            ax.tick_params(axis='x', labelsize=12, pad=1)
            ax.tick_params(axis='y', labelsize=12, pad=1)
            ax.tick_params(axis='z', labelsize=11, pad=9)
            #ax.set_zlim(0, 27000)
            ax.zaxis.set_major_formatter(FuncFormatter(scientific_format))
            

            # Exibindo o plot
            plt.show()

            # MESH CÚBICO FAST - GRAFICO 3D
            # Função para escurecer e aumentar a vivacidade do colormap
            def enhance_and_darken_cmap(cmap, gamma, scale):
                colors = cmap(np.linspace(0, 1, 256)) ** gamma
                darkened_colors = colors * scale
                darkened_colors[:, -1] = 1  # Manter a opacidade
                return LinearSegmentedColormap.from_list('enhanced_darkened_coolwarm', darkened_colors)

            # Criar um colormap com cores mais vivas e escurecidas
            coolwarm = plt.get_cmap('coolwarm')
            enhanced_darkened_coolwarm = enhance_and_darken_cmap(coolwarm, gamma=1.0, scale=1.0)

            plt.style.use('seaborn-v0_8-paper')

            # Criando a figura e os eixos 3D
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

            # Plotando a superfície 3D
            X, Y = np.meshgrid(x_divisions[1:], y_divisions[1:])
            surf = ax.plot_surface(X, Y, amplitude_fast, cmap=enhanced_darkened_coolwarm, edgecolor='black')
            cbar = fig.colorbar(surf, aspect=10.0, fraction=0.02, pad=0.018)
            cbar.set_label('Flux (neutrons.cm⁻².s⁻¹)', size=20, labelpad=8)
            cbar.ax.tick_params(axis='y', labelsize=12)
            cbar.ax.yaxis.set_major_formatter(FuncFormatter(scientific_format))
            # Adicionando rótulos aos eixos
            ax.set_xlabel('X position (cm)', size=20, labelpad=6)
            ax.set_ylabel('Y position (cm)', size=20, labelpad=6)
            #ax.set_zlabel('Flux (neutrons.cm⁻².s⁻¹)', size=20, labelpad=30)
            ax.tick_params(axis='x', labelsize=12, pad=1)
            ax.tick_params(axis='y', labelsize=12, pad=1)
            ax.tick_params(axis='z', labelsize=11, pad=9)
            ax.set_zlim(0, 175000)
            ax.zaxis.set_major_formatter(FuncFormatter(scientific_format))
            
            # Exibindo o plot
            plt.show()
            ################################################################################################################################
            # ESPECTRO DE ENERGIA

            plt.style.use('seaborn-v0_8-paper')
            # Escala logarítimica
            plt.xscale('log')
            plt.yscale('log')

            plt.plot(flux_spec_fuel_energy, flux_spec_fuel_mean, color='xkcd:red', linestyle='-', linewidth=1,marker='^', markersize=5, label='Inside the annular fuel ')
            plt.plot(flux_spec_central_energy, flux_spec_central_mean, color='xkcd:royal blue', linestyle='-', linewidth=1,marker='.', markersize=7, label='Above source')

            # Títulos e legenda
            plt.title('Flux Spectrum at Irradiation Positions', fontsize=24)
            plt.ylabel('Flux (neutrons.cm⁻².s⁻¹)', fontsize=20)
            plt.xlabel('Energy (eV)', fontsize=20)
            plt.legend(fontsize=20, loc='upper left')

            # Gridlines 
            plt.grid(True, which='both', axis='y', linestyle='--', linewidth=0.2, color='gray')
            plt.grid(True, which='both', axis='x', linestyle='--', linewidth=0.2, color='gray') # which='major'
            plt.tick_params(axis='both', which='major', labelsize=16)

            plt.tight_layout()
            plt.show() 

            ########################################################################################################################################
            ####           RADIAL           ######
            plt.style.use('seaborn-v0_8-paper')

            fig = plt.figure()

            ########################################################################################################################################
            # Convert the lists to NumPy arrays for element-wise negation
            radius_np = np.array(r_divisions_radial[1:])

            # Invert the values in radius_np
            inverted_radius = -radius_np

            # Criando o espelhamento
            #grid = gridspec.GridSpec(nrows=1, ncols=2, wspace=0, figure=fig) # hspace=0 (horizontal space)

            grafico = fig.add_subplot()

            # Gráficos
            #grafico = fig.add_subplot(grid[0, 1], zorder=3)
            #grafico.margins(0)
            #invert_grafico = fig.add_subplot(grid[0, 0], zorder=2, sharey=grafico) 
            #invert_grafico.margins(0)
            linewidth=1.5
            ########################################################################################################################################
            # Termico
            grafico.plot( np.array(flux_rc_r_termico), flux_rc_rad_termico, color='xkcd:red', linestyle='-', linewidth=linewidth, label='Thermal group')
            grafico.plot(-np.array(flux_rc_r_termico), flux_rc_rad_termico, color='xkcd:red', linestyle='-', linewidth=linewidth)

            # Fonte
            #grafico.plot( np.array(flux_r_fonte), flux_rad_fonte, color='xkcd:darkish green', linestyle='-', linewidth=0.5, label='Fonte')
            #grafico.plot(-np.array(flux_r_fonte), flux_rad_fonte, color='xkcd:darkish green', linestyle='-', linewidth=0.5)

            # Rapido
            grafico.plot( np.array(flux_rc_r_rapido), flux_rc_rad_rapido, color='xkcd:royal blue', linestyle='-',  linewidth=linewidth, label='Fast group')
            grafico.plot(-np.array(flux_rc_r_rapido), flux_rc_rad_rapido, color='xkcd:royal blue', linestyle='-', linewidth=linewidth)

            ########################################################################################################################################

            # y-axis label and legend
            plt.ylabel('Flux (neutrons.cm⁻².s⁻¹)', fontsize=20)
            plt.yscale('log')
            grafico.legend(fontsize=20)

            plt.tick_params(axis='both', which='major', labelsize=12)
            #plt.gca().yaxis.get_offset_text().set_size(12)
            #plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter())
            #plt.gca().yaxis.get_major_formatter().set_powerlimits((0, 0))

            # Centered x-axis label
            plt.xlabel( 'Radial Position (cm)', fontsize=20)

            # Centered title between the two subplots
            #plt.suptitle('Radial Flux Distribution', x=0.56, y=0.88, ha='center', fontsize=24)
            plt.title('Radial Flux Distribution',loc='center', fontsize=24)

            # expanding the y-axis limit
            #grafico.set_ylim(0, 1050)
            #invert_grafico.set_ylim(0, 5.0e+13)

            # clean up overlapping ticks and set axis to not visible
            #grafico.tick_params(axis='y', labelleft=False, left=False)
            #grafico.spines['left'].set_visible(False)
            #invert_grafico.spines['right'].set_visible(False)

            # Set x-axis limit in invert_mcnp to include zero point
            #invert_grafico.set_xlim(-radius_np[-1], 0)  # Adjust the range to include the zero point

            # Add gridlines to make zero point more visible
            #invert_grafico.grid(True, which='both', linestyle='--', linewidth=0.2, color='gray')
            grafico.grid(True, which='both', linestyle='--', linewidth=0.2, color='gray')
            

            plt.tight_layout()
            plt.show()


            ########################################################################################################################################
            #####         AXIAL         ######
            
            plt.style.use('seaborn-v0_8-paper')

            # Converter listas para arrays NumPy
            fluxo_axial_central_thermal = np.array(fluxo_axial_central_thermal)
            fluxo_axial_comb_thermal = np.array(fluxo_axial_comb_thermal)
            fluxo_z_central_thermal = np.array(fluxo_z_central_thermal)
            fluxo_z_comb_thermal = np.array(fluxo_z_comb_thermal)

            fluxo_axial_central_fast = np.array(fluxo_axial_central_fast)
            fluxo_axial_comb_fast = np.array(fluxo_axial_comb_fast)
            fluxo_z_central_fast = np.array(fluxo_z_central_fast)
            fluxo_z_comb_fast = np.array(fluxo_z_comb_fast)

            condicao_fuel_thermal_avarage = (fluxo_z_comb_thermal > -13.58) & (fluxo_z_comb_thermal < 9.16)
            condicao_central_thermal_avarage = (fluxo_z_central_thermal > 7.943) & (fluxo_z_central_thermal < 12.547)
            condicao_fuel_fast_avarage = (fluxo_z_comb_fast > -13.58) & (fluxo_z_comb_fast < 9.16)
            condicao_central_fast_avarage = (fluxo_z_central_fast > 7.943) & (fluxo_z_central_fast < 12.547)
            fluxo_thermal_fuel_avarage = fluxo_axial_comb_thermal[condicao_fuel_thermal_avarage]
            fluxo_thermal_central_avarage = fluxo_axial_central_thermal[condicao_central_thermal_avarage]
            fluxo_fast_fuel_avarage = fluxo_axial_comb_fast[condicao_fuel_fast_avarage]
            fluxo_fast_central_avarage = fluxo_axial_central_fast[condicao_central_fast_avarage]
            print()
            print('Altura considerada acima da fonte:', (12.547-7.943), 'Posições: 7.943 para 12.547')
            print()
            print('Altura considerada centro do combustível:', (13.58+9.16), 'Posições: -13.58 para 9.16')
            print()
            print('Média de fluxo térmico dentro do combustível:', "{:e}".format(np.mean(fluxo_thermal_fuel_avarage)))
            print('')
            print('Média de fluxo rápido dentro do combustível:', "{:e}".format(np.mean(fluxo_fast_fuel_avarage)))
            print('')
            print('Média de fluxo térmico acima da fonte:', "{:e}".format(np.mean(fluxo_thermal_central_avarage)))
            print('')
            print('Média de fluxo rápido acima da fonte:', "{:e}".format(np.mean(fluxo_fast_central_avarage)))
            print()

            condicao_thermal = fluxo_z_central_thermal > 0
            fluxo_axial_central_thermal_positivos = fluxo_axial_central_thermal[condicao_thermal]
            fluxo_z_central_thermal_positivos = fluxo_z_central_thermal[condicao_thermal]
            fluxo_axial_central_thermal_negativos = fluxo_axial_central_thermal[~condicao_thermal]
            fluxo_z_central_thermal_negativos = fluxo_z_central_thermal[~condicao_thermal]

            condicao_fast = fluxo_z_central_fast > 0
            fluxo_axial_central_fast_positivos = fluxo_axial_central_fast[condicao_fast]
            fluxo_z_central_fast_positivos = fluxo_z_central_fast[condicao_fast]
            fluxo_axial_central_fast_negativos = fluxo_axial_central_fast[~condicao_fast]
            fluxo_z_central_fast_negativos = fluxo_z_central_fast[~condicao_fast]

            plt.plot(fluxo_axial_central_thermal_positivos, fluxo_z_central_thermal_positivos, color='xkcd:burnt orange', linestyle='-',  linewidth=linewidth, label='Central Thermal')
            plt.plot(fluxo_axial_central_thermal_negativos, fluxo_z_central_thermal_negativos, color='xkcd:burnt orange', linestyle='-',  linewidth=linewidth)
            plt.plot(fluxo_axial_central_fast_positivos, fluxo_z_central_fast_positivos, color='xkcd:cerulean', linestyle='-',  linewidth=linewidth, label='Central Fast')
            plt.plot(fluxo_axial_central_fast_negativos, fluxo_z_central_fast_negativos, color='xkcd:cerulean', linestyle='-',  linewidth=linewidth)
            
            
            plt.plot(fluxo_axial_comb_thermal, fluxo_z_comb_thermal, color='xkcd:red', linestyle='--',  linewidth=linewidth, label='Comb Thermal')
            plt.plot(fluxo_axial_comb_fast, fluxo_z_comb_fast, color='xkcd:royal blue', linestyle='--',  linewidth=linewidth, label='Comb Fast')

            # Legendas
            plt.legend(fontsize=20)

            plt.tick_params(axis='both', which='major', labelsize=12)
            #plt.gca().xaxis.set_major_formatter(ticker.ScalarFormatter())
            #plt.gca().xaxis.get_major_formatter().set_powerlimits((0, 0))

            # Títulos e escala
            plt.xscale('log')
            plt.title('Axial flux distribution at Irradiation Positions', fontsize=24)
            plt.ylabel('Axial position (cm)', fontsize=20)
            plt.xlabel('Flux (neutrons.cm⁻².s⁻¹)', fontsize=20)
            plt.gca().xaxis.get_offset_text().set_size(12)

            # Add gridlines to make zero point more visible
            plt.grid(True, which='both', linestyle='--', linewidth=0.2, color='gray')

            plt.tight_layout()
            plt.show()