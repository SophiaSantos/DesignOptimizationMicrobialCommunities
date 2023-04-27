import csv
import os

from reframed import load_cbmodel, save_cbmodel, Environment, pFBA
from reframed.community.model import Community
from optimModels.utils.utils import fix_exchange_reactions_model
from optimModels.simulation.simul_problems import StoicSimulationProblem
from reframed import SteadyCom

community_name = []
organism_list = []
o = []

basePath = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/PycharmProjects/DD_DeCaF/Tests/examples/models/nitro"
list_models = []
constraints = {}

for file in os.listdir(basePath):
    if file.endswith("smetana_fbc2.xml"):
        SBML_FILE = (os.path.join(basePath, file))
        Id = file
        model = load_cbmodel(SBML_FILE, exchange_detection='R_EX_')
        model.biomass_reaction = "R_e_Biomass__cytop"
        list_models.append(model)

#print(list_models)

community = Community('Nitro', models=list_models)
merged = community.merge_models()
Environment.empty(merged, inplace=True)

print(merged.id)


save_cbmodel(merged, basePath + "/commModel_nitro.xml", flavor='fbc2')

minimal_medium = ['R_EX_M_so4_e', 'R_EX_M_o2_e', 'R_EX_M_pi_e', 'R_EX_M_fe2_e', 'R_EX_M_nh4_e']

carbon_source = ['R_EX_M_co2_e']

so4 = ['R_EX_M_so4_e']

KO = ['R_EX_M_no2_e']


for r_id, rxn in merged.reactions.items():
    if r_id.startswith('R_EX_'):
        rxn.lb = -1000 if rxn.lb is None else rxn.lb
        rxn.ub = 1000 if rxn.ub is None else rxn.ub

print(merged.biomass_reaction)

#simulProb = StoicSimulationProblem(merged, objective={'R_Community_Growth': 1}, method="pFBA")

#pfba = simulProb.simulate()

#essential = simulProb.find_essential_drains()
#print("Essential; ", essential)

for i in minimal_medium:
    constraints.update ({i:(-1000, 0)})

for i in KO:
    constraints.update({i: (0, 0)})

for i in carbon_source:
    constraints.update({i: (-0.435, 0)})

#for i in so4:
#    constraints.update({i: (-5, 0)})


print(constraints)

sol = pFBA(merged, objective={'community_growth': 1}, constraints=constraints)

print("Biomass: " + str(sol.show_values(pattern="Biomass", sort=True)))


#for r_id, rxn in merged.reactions.items():
#    if r_id.startswith('R_EX_'):
#        if sol.show_values(pattern=r_id) != 0:
#            print('Rxn: ' + str(r_id) + " Flux: " + str(sol.show_values(pattern=r_id)))

print('================================================================')

#comm_sol = pFBA(merged)
#org_fluxes = community.split_fluxes(comm_sol.values)

for x in list_models:
    name = x.id
    if name in community.organisms:
        print(sol.show_values(pattern='_e_'+ name))
