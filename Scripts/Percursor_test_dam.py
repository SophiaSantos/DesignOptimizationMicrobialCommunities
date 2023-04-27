import os

import cobra
from cobra import Model, Reaction, Metabolite
from cobra.io import read_sbml_model

basePath = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/PycharmProjects/DD_DeCaF/Tests/examples/models/sophia_models"

SBML_FILE = (os.path.join(basePath, 'model_dam_v3.xml'))
model = read_sbml_model(SBML_FILE)

constraints = {}

KO=[]

for i in model.reactions:
    if '__dra' in i.id:
        i.bounds = (0, 1000)

bRXN = "R_e_Biomass__cytop"

carbon_source = ['EX_C00095__dra']

minimal_medium = ['EX_C00014__dra', 'EX_C00087__dra', 'EX_C00009__dra', 'EX_C00001__dra']

essential = ['EX_C00504__dra', 'EX_C00253__dra', 'EX_C00255__dra']

model.biomass_reaction = bRXN

ATPM = ['ATPM__cytop']

for i in minimal_medium:
    model.reactions.get_by_id(i).lower_bound = -1000
    model.reactions.get_by_id(i).upper_bound = 1000

for i in essential:
    model.reactions.get_by_id(i).lower_bound = -1
    model.reactions.get_by_id(i).upper_bound = 0

for i in carbon_source:
    model.reactions.get_by_id(i).lower_bound = -0.2
    model.reactions.get_by_id(i).upper_bound = 0

for i in ATPM:
    model.reactions.get_by_id(i).upper_bound = 1.39
    model.reactions.get_by_id(i).lower_bound = 1.39

for i in KO:
    model.reactions.get_by_id(i).upper_bound = 0
    model.reactions.get_by_id(i).lower_bound = 0

model.biomass_reaction = bRXN

#model.objective = "e_Biomass_anaerobic__cytop"

solution = model.optimize()

print(model.summary())

print(model.metabolites.C00080__pmf_.summary())

#cobra.flux_analysis.find_blocked_reactions(model)