import os
import itertools

from framed import load_cbmodel, FVA, pFBA
from optimModels.utils.utils import fix_exchange_reactions_model
from optimModels.simulation.simul_problems import StoicSimulationProblem


basePath = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/PycharmProjects/Minimal_Medium_Community_Analysis/examples/models/sophia_models/"

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


SBML_FILE = (os.path.join(basePath, 'model_tye_v1.xml'))
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
        #        rxn.lb = 0
        exchange_list.append(r_id)

bRXN = "R_e_Biomass__in"

pyruvate = 'R_EX_C00022__in'

carbon_source_list = ['R_EX_C00022__in', 'R_EX_C00033__in', 'R_EX_C00256__in']
carbon_source = ['R_EX_C00022__in']

essential = ['R_EX_C14818__in', 'R_EX_C14819__in', 'R_EX_C06232__in', 'R_EX_C00117__in']

minimal_medium = ['R_EX_C00014__in', 'R_EX_C00009__in', 'R_EX_C00087__in', 'R_EX_C00059__in', 'R_EX_C00094__in', 'R_EX_C00022__in']


newModel.biomass_reaction = bRXN


for i in minimal_medium:
    constraints.update({i: (-1000, 1000)})

#for i in carbon_source_list:
#    constraints.update({i: (-10, 1000)})

for i in essential:
    constraints.update({i: (-1000, 1000)})

print(constraints)

for j in carbon_source:

    constraints.update({j: (-12, 0)})
    print(constraints)

    simulProb = StoicSimulationProblem(newModel, objective={bRXN: 1}, method="pFBA", constraints=constraints)
    pfba = simulProb.simulate()

    if pfba.solverStatus == 1:

        print("===============================================================")
        print("Biomass: " + str(pfba.get_fluxes_distribution()[bRXN]))
        print("===============================================================")

        constraints.update({i: (0, 0)})

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


