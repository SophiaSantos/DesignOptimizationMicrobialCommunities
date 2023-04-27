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
#print("Build a community with a total of {} reactions.".format(len(com.reactions)))

reactions = ["EX_co2_m", "EX_o2_m", "EX_so4_m", "EX_nh4_m", "EX_pi_m", "EX_fe2_m"]
fluxes = [0.425, 1000, 1000, 1000, 1000, 1000]
#fluxes = [12, 50, 1000, 1000, 1000, 1000, 1000, 1000]

candidate_medium = pd.DataFrame({"reaction": reactions, "flux": fluxes})

medium1 = pd.Series(fluxes, reactions)

print(medium1)

#print(candidate_medium)

print(com.objective.expression)

single = com.optimize_all()
print(single)


print('========= MICOM =========================')
micom_sol = com.optimize(fluxes=True, pfba=True)
print(micom_sol.growth_rate)
print(micom_sol.members.growth_rate)
print('====================================')
rates = micom_sol.members.growth_rate.drop("medium")  # extracellular medium has no growth rate

med = minimal_medium(com, 0.1*micom_sol.growth_rate, min_growth=rates, exports=True)

com.medium=med

print(com.medium)

micom_sol_1 = com.optimize(fluxes=True, pfba=True)
print(micom_sol.growth_rate)
print(micom_sol.members.growth_rate)
print('====================================')

strategies = ['original', 'lmoma']


for strategy in strategies:
    print('========= OPTCOM =========================')
    optcom_sol = com.optcom(strategy=strategy,  fluxes=True)
    print(strategy)
    print(optcom_sol.growth_rate)
    print(optcom_sol.members.growth_rate)
    print('====================================')
    rates_1 = optcom_sol.members.growth_rate.drop("medium")  # extracellular medium has no growth rate

    med_1 = minimal_medium(com, 0.1*optcom_sol.growth_rate, min_growth=rates_1, exports=True)

    com.medium=med_1

    optcom_sol1 = com.optcom(strategy = strategy, min_growth=0.0075, fluxes=True)

    print(med_1)
    print(optcom_sol1.growth_rate)
    print(optcom_sol1.members.growth_rate)

    sol = com.cooperative_tradeoff(fraction=0.5, fluxes=True, pfba=True)

    print('====================================')
