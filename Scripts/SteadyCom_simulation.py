import os

from reframed import load_cbmodel
from reframed import Community, Environment
from reframed import SteadyCom, SteadyComVA

basePath = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/Doutoramento/Dissertation/models/community_models"

constraints = {}
list_models = []

for file in os.listdir(basePath):
    if file.endswith(".xml"):
        SBML_FILE = (os.path.join(basePath, file))
        Id = file
        print(Id)
        model = load_cbmodel(SBML_FILE, flavor="fbc2")
        model.biomass_reaction = "R_e_Biomass__cytop"

        for r_id, rxn in model.reactions.items():
            if r_id.startswith('R_EX_'):
                rxn.lb = 0

        list_models.append(model)


community = Community('CV', models = list_models)
merged = community.merge_models()
Environment.empty(merged, inplace = True)

print(merged.metabolites.items())


minimal_medium = ['R_EX_M_s_e', 'R_EX_M_o2_e', 'R_EX_M_ppi_e',  'R_EX_M_fe2_e', 'R_EX_M_nh4_e', 'R_EX_M__e']

carbon_source = ['R_EX_M_co2_e']

biomass_exchange = ['R_EX_M_e_Biomass_e']

no2_pool = ['R_TR0000075_EXTMEM__outme_model_nvulgaris']


for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

#for i in no2_pool:
#    constraints.update({i:(0,0)})

#for i in community.reaction_map:
#    print(i[1])

print(constraints)

print('===============================================')
print('=== SteadyComm ===')

steady_sol = SteadyCom(community)
print(steady_sol)
print(steady_sol.exchange_map)

for i, j in steady_sol.exchange_map.items():
    constraints.update({'R_EX_'+str(i[1]):(0, 1000)})

for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

steady_sol2 = SteadyCom(community, constraints=constraints)

print(steady_sol2)
print(steady_sol2.exchange_map)

for i, j in steady_sol2.exchange_map.items():
    constraints.update({'R_EX_'+str(i[1]):(0, 1000)})

for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

steady_sol3 = SteadyCom(community, constraints=constraints)

print(steady_sol3)
print(steady_sol3.exchange_map)


for i, j in steady_sol3.exchange_map.items():
    constraints.update({'R_EX_'+str(i[1]):(0, 1000)})

for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

steady_sol4 = SteadyCom(community, constraints=constraints)

print(steady_sol4)
print(steady_sol4.exchange_map)

for i, j in steady_sol4.exchange_map.items():
    constraints.update({'R_EX_'+str(i[1]):(0, 1000)})

for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

steady_sol5 = SteadyCom(community, constraints=constraints)

print(steady_sol5)
print(steady_sol5.exchange_map)

for i, j in steady_sol5.exchange_map.items():
    constraints.update({'R_EX_'+str(i[1]):(0, 1000)})

for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

steady_sol6 = SteadyCom(community, constraints=constraints)

print(steady_sol6)
print(steady_sol6.exchange_map)

for i, j in steady_sol6.exchange_map.items():
    constraints.update({'R_EX_'+str(i[1]):(0, 1000)})

for i in minimal_medium:
    constraints.update ({i:(-1000,0)})

for i in carbon_source:
    constraints.update({i:(-0.421,0)})

print(constraints)

steady_sol7 = SteadyCom(community, constraints=constraints)

print(steady_sol7)
print(steady_sol7.print_internal_fluxes('model_neuropaeaV4'))
print(steady_sol7.exchange_map)

steady_solFVA = SteadyComVA(community, obj_frac=0.8, constraints=constraints)
print('===============================================')
print('=== SteadyCommFVA ===')
print(steady_solFVA)