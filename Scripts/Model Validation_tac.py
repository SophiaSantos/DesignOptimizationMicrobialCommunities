import os
import itertools

from framed import load_cbmodel
from optimModels.utils.utils import fix_exchange_reactions_model
from optimModels.simulation.simul_problems import StoicSimulationProblem


basePath = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/PycharmProjects/DD_DeCaF/Tests/examples/models/sophia_models"
basePath_r = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/PycharmProjects/DD_DeCaF/Tests/examples/"

met_list = []
rxn_list = []
met_dict = {}
flux_dict = {}
met_list_name=[]
exchange_dict = {}
constraints = {}
exchange_list = []
header = ['Rxn', 'Flux', 'Condition']
Rxn = []
Rxn_lists = []
Flux_lists = []
Flux = []
dict_data = []
len_list=[]
conditions = []

KO = []

SBML_FILE = (os.path.join(basePath, 'model_tac_v7.xml'))
model = load_cbmodel(SBML_FILE, exchange_detection_mode='R_EX_')
newModel = fix_exchange_reactions_model(model)

print(newModel.id)

for r_id, rxn in newModel.reactions.items():
    if rxn.is_exchange:
        rxn.lb = -1000 if rxn.lb is None else rxn.lb
        rxn.ub = 1000 if rxn.ub is None else rxn.ub

for r_id, rxn in newModel.reactions.items():
    if r_id.startswith('R_ATPM'):
        rxn.lb = 8.39
        rxn.ub = 8.39
    if rxn.is_exchange:
        rxn.lb = 0
#        rxn.lb = -10
        exchange_list.append(r_id)
    if r_id.startswith('R_sink_'):
        rxn.lb = 0
        exchange_list.append(r_id)

bRXN = "R_e_Biomass_aerobic__cytop"

carbon_source = ['R_EX_C00031__dra']

essential = ['R_EX_C00047__dra', 'R_EX_C00123__dra', 'R_EX_C00135__dra']

isoleucine = ['R_EX_C00407__dra']

valine = ['R_EX_C00183__dra']

minimal_medium = ['R_EX_C00014__dra', 'R_EX_C00059__dra', 'R_EX_C00001__dra', 'R_EX_C00009__dra', 'R_EX_C00153__dra', 'R_EX_C00175__dra', 'R_EX_C00255__dra']

aerobic = ['R_EX_C00007__dra', 'R_EX_C14818__dra']

newModel.biomass_reaction = bRXN

#for i in optim_medium:
#    constraints.update({i: (-0.5, 1000)})

for i in minimal_medium:
    constraints.update({i: (-1000, 1000)})

for i in essential:
    constraints.update({i: (-0.05, 1000)})

for i in isoleucine:
    constraints.update({i: (-0.07, 0)})

for i in valine:
    constraints.update({i: (-0.05, 0)})

for i in aerobic:
    constraints.update({i: (-3, 0)})

for i in KO:
    constraints.update({i: (0, 0)})


print(constraints)

met = newModel.metabolites.items()

print(len(met))

for j in carbon_source:

    constraints.update({j: (-0.3, 0)})
    print(constraints)

    simulProb = StoicSimulationProblem(newModel, objective={bRXN: 1}, method="pFBA", constraints=constraints)
    pfba = simulProb.simulate()

    if pfba.solverStatus == 1:

        print("===============================================================")
        print("Biomass: " + str(pfba.get_fluxes_distribution()[bRXN]))
        print("===============================================================")

        #fva = FVA(newModel, obj_percentage=0.1,reactions=['R_EX_C00760__dra', 'R_EX_C00014__dra', 'R_EX_C00087__dra', 'R_EX_C14818__dra','R_EX_C00080__dra', 'R_EX_C00007__dra', 'R_EX_C00009__dra', 'R_EX_C00288__dra','R_EX_C00993__dra', 'R_EX_C00033__dra'], constraints=constraints)

        essential = simulProb.find_essential_drains()

        constraints.update({j: (0, 0)})

        for r_id, rxn in newModel.reactions.items():
            # print(r_id + ": " + str(pfba.get_fluxes_distribution()[r_id]))
            if rxn.is_exchange or r_id.startswith('R_sink'):
                met_list.append(rxn.get_substrates())
                rxn_list.append(pfba.get_fluxes_distribution()[r_id])
                # print(r_id + ": " + str(pfba.get_fluxes_distribution()[r_id]))

        met_list_merged = list(itertools.chain(*met_list))
        met_list_name = [1] * len(met_list_merged)

        for m_id, met in newModel.metabolites.items():
            for k in met_list_merged:
                for l in range(0, len(met_list_merged)):
                    if met_list_merged[l] == m_id:
                        met_list_name[l] = str(met)

        for m in range(0, len(met_list_name)):
            exchange_dict.update({met_list_name[m]: rxn_list[m]})

        for key in exchange_dict:
            if exchange_dict[key] != 0:
                Rxn.append(key)
                Flux.append(exchange_dict[key])
                print('Rxn: ' + str(key) + ', Flux: ' + str(exchange_dict[key]))

    print('=========== Flux distribuition ================')
    print(pfba.get_fluxes_distribution())
    print('============ Essential Reactions ===============')
    print(essential)
    print('============ FVA results ===============')
    #print(fva)


#for r_id, rxn in newModel.reactions.items():
#    print(rxn)


