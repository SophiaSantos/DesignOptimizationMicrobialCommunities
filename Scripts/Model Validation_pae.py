import os
import itertools

from framed import load_cbmodel, FVA, blocked_reactions
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

SBML_FILE = (os.path.join(basePath, 'model_pae_v7.xml'))
model = load_cbmodel(SBML_FILE, exchange_detection_mode='R_EX_')
newModel = fix_exchange_reactions_model(model)

print(newModel.id)

for r_id, rxn in newModel.reactions.items():
    if rxn.is_exchange:
        rxn.lb = -1000 if rxn.lb is None else rxn.lb
        rxn.ub = 1000 if rxn.ub is None else rxn.ub

for r_id, rxn in newModel.reactions.items():
    if r_id.startswith('R_ATPM'):
        rxn.lb = 1.5
        rxn.ub = 1.5
    if rxn.is_exchange:
        rxn.lb = 0
 #       rxn.lb = -10
        exchange_list.append(r_id)
    if r_id.startswith('R_Sink_'):
        rxn.lb = 0
        exchange_list.append(r_id)

bRXN = "R_e_Biomass_anaerobic__cytop"

carbon_source = ['R_EX_C00011__dra']

carbon_source_list = ['R_EX_C00025__dra', 'R_EX_C00037__dra', 'R_EX_C00041__dra', 'R_EX_C00047__dra', 'R_EX_C00049__dra', 'R_EX_C00062__dra', 'R_EX_C00064__dra', 'R_EX_C00065__dra','R_EX_C00073__dra', 'R_EX_C00079__dra', 'R_EX_C00082__dra', 'R_EX_C00097__dra', 'R_EX_C00123__dra', 'R_EX_C00135__dra', 'R_EX_C00148__dra', 'R_EX_C00152__dra', 'R_EX_C00183__dra', 'R_EX_C00188__dra', 'R_EX_C00407__dra']

minimal_medium = ['R_EX_C00014__dra', 'R_EX_C00320__dra', 'R_EX_C00013__dra', 'R_EX_C00177__dra', 'R_EX_C00282__dra']

essential = ['R_EX_C00504__dra', 'R_EX_C00253__dra']

aerobic = ['R_EX_C00007__dra', 'R_EX_C14818__dra']

KO = ['R_EX_C00013__dra', 'R_EX_C00080__dra', 'R_EX_C01755__dra', 'R_EX_C00469__dra', 'R_EX_C00221__dra', 'R_EX_C00267__dra', 'R_EX_C00246__dra']


newModel.biomass_reaction = bRXN

#for i in optim_medium:
#    constraints.update({i: (-0.5, 1000)})

for i in minimal_medium:
    constraints.update({i: (-1000, 1000)})

for i in aerobic:
    constraints.update({i: (0, 0)})

for i in essential:
    constraints.update({i: (-5, 0)})

for i in carbon_source:
    constraints.update({i: (-2, 0)})
    
#for i in KO:
#    constraints.update({i: (0, 1000)})

#for i in carbon_source_list:
#    constraints.update({i: (0, 1000)})

print(constraints)

met = newModel.metabolites.items()

print(len(met))

for j in carbon_source:

    constraints.update({j: (-2, 0)})
    print(constraints)

    simulProb = StoicSimulationProblem(newModel, objective={bRXN: 1}, method="pFBA", constraints=constraints)
    pfba = simulProb.simulate()

    if pfba.solverStatus == 1:

        print("===============================================================")
        print("Biomass: " + str(pfba.get_fluxes_distribution()[bRXN]))
        print("===============================================================")

        fva = FVA(newModel, obj_percentage=0.1, constraints=constraints)

        essential = simulProb.find_essential_drains()

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

    constraints.update({j: (0, 0)})

    print('=========== Flux distribuition ================')
    print(pfba.get_fluxes_distribution())
    print('============ Essential Reactions ===============')
    print(essential)
    print('============ FVA results ===============')
    #for r, f in fva.items():
    #    if f == [0, 0]:
    #        print(r)


#for r_id, rxn in newModel.reactions.items():
#    print(rxn)

blocked = blocked_reactions(newModel)
print(blocked)

