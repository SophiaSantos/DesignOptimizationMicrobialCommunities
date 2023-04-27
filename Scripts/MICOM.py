import os
import pandas as pd
import numpy as np

from micom.workflows import fix_medium
from micom.media import minimal_medium

from micom import Community

basePath = "C:/Users/Sophia Santos/OneDrive - Universidade do Minho/CEB/PycharmProjects/DD_DeCaf/Tests/examples/models/nitro/"

for file in os.listdir(basePath):
    if file.endswith("taxonomy.csv"):
        CSVFILE = (os.path.join(basePath, file))
        taxonomy = pd.read_csv(CSVFILE, delimiter=";")

com = Community(taxonomy)
#print("Build a community with a total of {}
#print(com.reactions)

reactions = ["EX_co2_m", "EX_o2_m", "EX_so4_m", "EX_nh4_m", "EX_pi_m", "EX_fe2_m"]
fluxes = [0.04, 1000, 1000, 1000, 1000, 1000]

candidate_medium= pd.Series(fluxes, reactions)

print(candidate_medium)

#print(candidate_medium)

print(com.objective.expression)
fba = com.optimize()
pfba = com.optimize(fluxes = True, pfba = True)

sol = com.cooperative_tradeoff(fluxes = True, pfba = True)
sol3 = com.cooperative_tradeoff()

sol1 = com.cooperative_tradeoff(min_growth=0.0075, fluxes=True)  # single value
sol2 = com.cooperative_tradeoff(min_growth=[0.0075, 0.01], fluxes=True)  # one value for each taxa

#sols = com.cooperative_tradeoff(fraction=np.arange(0.001, 1.01, 0.001))

#rates = sols.solution.apply(lambda x: x.members.growth_rate)

med = minimal_medium(com, 0.8)
minMed = minimal_medium(com, 0.8, minimize_components=True)

rates_1 = sol1.members.growth_rate.drop("medium")  # extracellular medium has no growth rate

med_1 = minimal_medium(com, 0.1*sol1.growth_rate, min_growth=rates_1, exports=True)

com.medium = med_1
pfba_medium = com.optimize(fluxes = True, pfba = True)

print(fba.members)

sol1.fluxes.to_csv(basePath+'fluxes_MICOM.csv', sep=';')

print('Solução 1 e Solução 2')
print(sol1.members, sol2.members)
print(sol1.growth_rate, sol2.growth_rate)
#print(sols)

#print(rates)
print(minMed)
#print(minMed)

print(med_1)
