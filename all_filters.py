filter_list = ['SA', '','AF', "AU","MA", "PO","OR","NA","EU", "IO"]

all_data = []

for f in filter_list:
    df = perform_simulation(10000, f)

    all_data.append(df)
df = pd.concat(all_data)

df.to_csv('all_sim.csv', index = False)
