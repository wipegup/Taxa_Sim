import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def build_info_df(filter = 'SA'):
    df = pd.read_excel("IOC_Names_File_Plus-9.2.xlsx")
    df.columns = ["index", "rank","junk_1", "junk_2", "junk_3", "name", "author", "location", "junk_4","junk_5", "junk_6"]
    df = df.drop([c for c in df.columns if c[:4] == "junk"], axis = 'columns')
    df = df.drop("index", axis = 'columns')
    df = df.drop([0,1], axis = 'rows')
    df = df.drop(df[df['rank'].isin(['Blank', 'TAXON', 'ssp'])].index, axis = 'rows')

    info_dict = []
    col_name = "in_area"#+ filter.lower()
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
            in_area = row['location'].find(filter) >-1

            record = {"order":order, "family":family, "genus":genus, "species": name, col_name: in_area}
            info_dict.append(record)

    df = pd.DataFrame(info_dict)
    return df

def find_distribution(threshold, filter):
    col_name = "in_area"#+ filter.lower()

    df = build_info_df()
    gen = df.groupby("genus")[col_name].mean().sort_values()
    keep = list(gen[gen >= threshold].index)
    return df[df['genus'].isin(keep)].groupby('genus')['species'].count().value_counts()

def aggregate_info(filter = 'SA', threshold = .75):
    col_name = "in_area"#+ filter.lower()

    df = build_info_df(filter)

    gen = df.groupby("genus")[col_name].mean().sort_values()
    threshes = gen[gen > 0].value_counts().sort_index().index
    plt.hist(threshes, bins = 30)

    fam = df.groupby("family")[col_name].mean()
    # fam.hist(bins = 30)
    sa_fam = fam[fam > threshold]
    s_f = df[df['family'].isin(sa_fam.index)].groupby('family')['species'].count()
    s_f_g = []
    for f in s_f.index:
        gen = df[df['family'] == f]['genus'].value_counts().count()
        s_f_g.append({'name':f, 'sub_taxa': gen})
    df_sfg = pd.DataFrame(s_f_g)
    df_sfg['species'] = s_f.values
    df_sfg['steps'] = 2
    df_sfg['kind'] = 'Family'
    # df_sfg

    order = df.groupby("order")[col_name].mean()
    # order.hist(bins = 30)
    sa_order = order[order > .75]
    s_o = df[df['order'].isin(sa_order.index)].groupby('order')['species'].count()

    s_o_g = []
    for o in s_o.index:
        fam = df[df['order'] == o]['family'].value_counts().count()
        s_o_g.append({'name':o, 'sub_taxa': fam})
    df_sof = pd.DataFrame(s_o_g)
    df_sof['species'] = s_o.values
    df_sof['steps'] = 3
    df_sof['kind'] = 'Order'

    df = pd.concat([df_sfg, df_sof])
    df = df.reset_index()
    return df

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

def perform_simulation(n = 10, filter = 'SA'):
    to_ret = []
    for t in np.arange(.5,1,.05):

        df = aggregate_info(filter, t)

        default_dist = find_distribution(t, filter)
        obs = sum(default_dist)
        default_dist = default_dist / obs
        print(df.shape)
        preds = []
        for i, idx in enumerate(df.index):
            if i% 5 == 0:
                print(i)
            row = df.loc[idx, :]
            species = row['species']
            steps = row['steps']
            pred = sub_taxa_distribution(n, steps, species, obs, default_dist)
            preds.append(pred)
        df['preds'] = preds
        df['thresh'] = t
        df['filter'] = filter
        to_ret.append(df)
    return pd.concat(to_ret)

def print_output(lg_df):
    fig, axs = plt.subplots(4,5, figsize = (10,10))

    for idx, i in enumerate(lg_df.index):
        row = idx // 5
        col = idx % 5
        ax = axs[row, col]
        ax.hist(lg_df.loc[i,'preds'], bins = 30)
        ax.vlines(lg_df.loc[i, 'sub_taxa'], *ax.get_ylim(), color = 'red')
        ax.set_title(lg_df.loc[i,'name'])
    plt.tight_layout()

dist_dict = {1: 330,
 2: 130,
 3: 56,
 4: 44,
 5: 36,
 6: 22,
 7: 19,
 8: 15,
 9: 11,
 11: 11,
 10: 8,
 15: 6,
 12: 5,
 14: 5,
 13: 3,
 16: 2,
 17: 2,
 18: 2,
 22: 2,
 27: 2,
 30: 2,
 44: 1,
 19: 1,
 20: 1,
 21: 1,
 24: 1,
 25: 1,
 31: 1,
 32: 1,
 34: 1,
 35: 1,
 49: 1}


# df.to_csv("SA_info.csv", index = False)
