import calliope

model = calliope.read_netcdf("results/model_results_B1.nc")

print("Fossil output:", model.results.fossil_output_debug.values)
print("Total generation:", model.results.total_generation_debug.values)
# print("Non-fossil output:", model.results.non_fossil_output_debug.values)
# print("Domestic demand:", model.results.domestic_demand_debug.values)

