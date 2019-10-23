import pandas as pd
import numpy as np
filter_list = ['SA', '','AF', "AU","MA", "PO","OR","NA","EU", "IO"]

all_data = []

def build_df():
    df = pd.read_excel("IOC_Names_File_Plus-9.2.xlsx")
    df.columns = ["index", "rank","junk_1", "junk_2", "junk_3", "name", "author", "location", "junk_4","junk_5", "junk_6"]
    df = df.drop([c for c in df.columns if c[:4] == "junk"], axis = 'columns')
    df = df.drop("index", axis = 'columns')
    df = df.drop([0,1], axis = 'rows')
    df = df.drop(df[df['rank'].isin(['Blank', 'TAXON', 'ssp'])].index, axis = 'rows')

    info_dict = []
    for idx in df.index:
        row = df.loc[idx, :]
        rank = row['rank']

        if rank == "ORDER":
            order = row['name'].split(" ")[1]

        if rank == "Family":
            family = row['name'].split(" ")[1]

        if rank == "Genus":
            genus = row['name']

        if rank == "Species":
            name = row['name'].split(" ")[1]
            location = row['location']

            record = {"order":order, "family":family, "genus":genus, "species": name, 'location': location}
            info_dict.append(record)

    df = pd.DataFrame(info_dict)

    return df

def filter_df(all_df, filt):
    df = all_df.copy()
    def find_filter(string, filt):
        return string.find(filt) >-1

    df['in_area'] = df['location'].apply(find_filter, args = (filt,))

    return df

def aggregate_info(filtered, threshold):
    col_name = 'in_area'
    gen = filtered.groupby("genus")[col_name].mean().sort_values()
    # threshes = gen[gen > 0].value_counts().sort_index().index
    # plt.hist(threshes, bins = 30)

    fam = filtered.groupby("family")[col_name].mean()
    # fam.hist(bins = 30)
    sa_fam = fam[fam > threshold]
    s_f = filtered[filtered['family'].isin(sa_fam.index)].groupby('family')['species'].count()
    s_f_g = []
    for f in s_f.index:
        gen = filtered[filtered['family'] == f]['genus'].value_counts().count()
        s_f_g.append({'name':f, 'sub_taxa': gen})
    df_sfg = pd.DataFrame(s_f_g)
    df_sfg['species'] = s_f.values
    df_sfg['steps'] = 2
    df_sfg['kind'] = 'Family'
    # df_sfg

    order = filtered.groupby("order")[col_name].mean()
    # order.hist(bins = 30)
    sa_order = order[order > threshold]
    s_o = filtered[filtered['order'].isin(sa_order.index)].groupby('order')['species'].count()

    s_o_g = []
    for o in s_o.index:
        fam = filtered[filtered['order'] == o]['family'].value_counts().count()
        s_o_g.append({'name':o, 'sub_taxa': fam})
    df_sof = pd.DataFrame(s_o_g)
    df_sof['species'] = s_o.values
    df_sof['steps'] = 3
    df_sof['kind'] = 'Order'

    new_df = pd.concat([df_sfg, df_sof])
    new_df = new_df.reset_index()
    return new_df

def find_distribution(filtered, threshold):
    col_name = "in_area"#+ filter.lower()

    gen = filtered.groupby("genus")[col_name].mean().sort_values()
    keep = list(gen[gen >= threshold].index)
    return filtered[filtered['genus'].isin(keep)].groupby('genus')['species'].count().value_counts()
def sample_distribution( obs, dist):
    return np.random.choice(dist.index, p = dist.values, size = obs, replace = True)

def num_species(steps, dist):
    if steps == 1:
        return np.random.choice(dist)
    else:
        sub_t = np.random.choice(dist)
        return sum([num_species(steps-1, dist) for i in range(sub_t)])


def predict_sub_taxa(steps, species, obs, default_dist):
    sub_taxa = 0
    total = 0
    while total < species:
        sub_taxa += 1
        dist = sample_distribution(obs, default_dist)
        total += num_species((steps - 1), dist)
    return sub_taxa


def sub_taxa_distribution(n, steps, species, obs, default_dist):
    preds = []
    for i in range(n):
        preds.append(predict_sub_taxa(steps, species, obs, default_dist))
    return preds

def perform_simulation(all_df, filter, n = 10000):
    print("Filt", filter)
    to_ret = []
    filtered = filter_df(all_df, filter)
    for t in np.arange(.5, 1, .05):
        print("Thresh", t)
        aggregated = aggregate_info(filtered, t)
        default_dist = find_distribution(filtered, t)
        obs = sum(default_dist)
        default_dist = default_dist / obs
        preds = []
        for i, idx in enumerate(aggregated.index):
            if i% 5 == 0:
                print(i)
            row = aggregated.loc[idx, :]
            species = row['species']
            steps = row['steps']
            pred = sub_taxa_distribution(n, steps, species, obs, default_dist)
            preds.append(pred)
        aggregated['preds'] = preds
        aggregated['thresh'] = t
        aggregated['filter'] = filter
        to_ret.append(aggregated)
    return pd.concat(to_ret)


all_df = build_df()
all_data = []
for f in filter_list:
    sim = perform_simulation(all_df, f, 10000)
    all_data.append(sim)

df = pd.concat(all_data)
df.to_csv('all_sim.csv', index = False)

# col_name = 'in_area'
# for f in filter_list:
#     filtered = filter_df(all_df, f)
#     for t in np.arange(.5,1, .05):
#         gen = df.groupby("genus")[col_name].mean().sort_values()

#
# for f in filter_list:
#     df = perform_simulation(10000, f)
#
#     all_data.append(df)
# df = pd.concat(all_data)
#
# df.to_csv('all_sim.csv', index = False)
