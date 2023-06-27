import pandas as pd
HP_PEC_NToma = pd.read_csv("HP_PEC_NToma_matrix.csv", index_col = 0)
HP_PEC_PToma = pd.read_csv("HP_PEC_PToma_matrix.csv", index_col = 0)
IL4_PEC_NToma = pd.read_csv("IL4_PEC_NToma_matrix.csv", index_col = 0)
IL4_PEC_PToma = pd.read_csv("IL4_PEC_PToma_matrix.csv", index_col = 0)
SS_PEC_NToma = pd.read_csv("SS_PEC_NToma_matrix.csv", index_col = 0)
SS_PEC_PToma = pd.read_csv("SS_PEC_PToma_matrix.csv", index_col = 0)
ThioIL4PCD206_PEC_PToma = pd.read_csv("ThioIL4PCD206_PEC_PToma_matrix.csv", index_col = 0)
ThioPCD206_PEC_PToma = pd.read_csv("ThioPCD206_PEC_PToma_matrix.csv", index_col = 0)

intersect_genes = list(set(HP_PEC_NToma.index) & set(HP_PEC_PToma.index) & set(IL4_PEC_NToma.index) & set(IL4_PEC_PToma.index) & 
                       set(SS_PEC_NToma.index) & set(SS_PEC_PToma.index) & set(ThioIL4PCD206_PEC_PToma.index) & set(ThioPCD206_PEC_PToma.index))

HP_PEC_NToma_inter = HP_PEC_NToma.loc[intersect_genes].T
HP_PEC_NToma_inter["condition"] = "HP_PEC_NToma"
HP_PEC_PToma_inter = HP_PEC_PToma.loc[intersect_genes].T
HP_PEC_PToma_inter["condition"] = "HP_PEC_PToma"

IL4_PEC_NToma_inter = IL4_PEC_NToma.loc[intersect_genes].T
IL4_PEC_NToma_inter["condition"] = "IL4_PEC_NToma"
IL4_PEC_PToma_inter = IL4_PEC_PToma.loc[intersect_genes].T
IL4_PEC_PToma_inter["condition"] = "IL4_PEC_PToma"

SS_PEC_NToma_inter = SS_PEC_NToma.loc[intersect_genes].T
SS_PEC_NToma_inter["condition"] = "SS_PEC_NToma"
SS_PEC_PToma_inter = SS_PEC_PToma.loc[intersect_genes].T
SS_PEC_PToma_inter["condition"] = "SS_PEC_PToma"

ThioIL4PCD206_PEC_PToma_inter = ThioIL4PCD206_PEC_PToma.loc[intersect_genes].T
ThioIL4PCD206_PEC_PToma_inter["condition"] = "ThioIL4PCD206_PEC_PToma"
ThioPCD206_PEC_PToma_inter = ThioPCD206_PEC_PToma.loc[intersect_genes].T
ThioPCD206_PEC_PToma_inter["condition"] = "ThioPCD206_PEC_PToma"

data_inter_all = pd.concat([HP_PEC_NToma_inter, HP_PEC_PToma_inter, IL4_PEC_NToma_inter, IL4_PEC_PToma_inter, SS_PEC_NToma_inter, SS_PEC_PToma_inter, ThioIL4PCD206_PEC_PToma_inter, ThioPCD206_PEC_PToma_inter])

import numpy as np
X = data_inter_all.loc[:,data_inter_all.columns[:-1]]
y = data_inter_all['condition']

from sklearn.model_selection import train_test_split
import random
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 0)

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
rfc = RandomForestClassifier(n_estimators=200, random_state=0)
rfc.fit(X_train, y_train)
y_pred = rfc.predict(X_test)
print('Model accuracy score with 200 decision-trees : {0:0.4f}'. format(accuracy_score(y_test, y_pred)))

feature_scores = pd.Series(rfc.feature_importances_, index=X_train.columns).sort_values(ascending=False)

feature_scores_top40 = feature_scores.loc[feature_scores.index[:30]]
import seaborn as sns
import matplotlib.pyplot as plt
f, ax = plt.subplots(figsize=(5, 8))
ax = sns.barplot(x=feature_scores_top40, y=feature_scores_top40.index)
ax.set_yticklabels(feature_scores_top40.index)
ax.set_xlabel("Gini importance score")
ax.set_ylabel("Top 30 genes")
plt.show()

data_inter_all_expressed = pd.DataFrame(np.array(data_inter_all[feature_scores_top40.index.tolist()] >0, dtype=int), 
                                        columns=feature_scores_top40.index.tolist(), index=data_inter_all.index)
data_inter_all_expressed['condition'] = data_inter_all["condition"]
data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)

high_exp_gene = np.array(feature_scores_top40.index)[(data_inter_all[feature_scores_top40.index.tolist() + ["condition"]].groupby("condition").mean(0)).mean() >= 1.5]
data_inter_all_expressed_percent = data_inter_all_expressed_percent[high_exp_gene]
data_inter_dotplot = pd.melt(data_inter_all_expressed_percent*100, var_name='Top 40 genes', value_name="Percent Expressed")
data_inter_dotplot['Conditions'] = np.array(data_inter_all_expressed_percent.index.tolist()*len(high_exp_gene))
data_inter_dotplot['Average Expression'] = pd.melt(data_inter_all[list(high_exp_gene) + ["condition"]].groupby("condition").mean(0))['value']
data_inter_dotplot["Average Expression"] = [5 if AE >5 else AE for AE in data_inter_dotplot["Average Expression"]]
data_inter_dotplot["Average Expression"] = [1 if AE <1 else AE for AE in data_inter_dotplot["Average Expression"]]

import altair as alt
range_ = ['#990000', '#CC0000', '#FF0000', '#FF3333', '#FF6666', '#FF9999', '#FFCCCC']
alt.Chart(data_inter_dotplot).mark_circle().encode(
    x=alt.X('Top 40 genes:O', sort=high_exp_gene),
    y = 'Conditions:O',
    size = 'Percent Expressed:Q',
    color = alt.Color(
            'Average Expression:Q',
            scale = alt.Scale(domain= [max(data_inter_dotplot['Average Expression']),min(data_inter_dotplot['Average Expression'])], range = range_))
)

# 1. 
intersect_genes = list(set(HP_PEC_NToma.index) & set(HP_PEC_PToma.index) & set(IL4_PEC_NToma.index) & set(IL4_PEC_PToma.index) & 
                       set(SS_PEC_NToma.index) & set(SS_PEC_PToma.index) & set(ThioIL4PCD206_PEC_PToma.index) & set(ThioPCD206_PEC_PToma.index))

HP_PEC_NToma_inter = HP_PEC_NToma.loc[intersect_genes].T
HP_PEC_NToma_inter["condition"] = 'TRM'
HP_PEC_PToma_inter = HP_PEC_PToma.loc[intersect_genes].T
HP_PEC_PToma_inter["condition"] = 'TRM'

IL4_PEC_NToma_inter = IL4_PEC_NToma.loc[intersect_genes].T
IL4_PEC_NToma_inter["condition"] = 'TRM'
IL4_PEC_PToma_inter = IL4_PEC_PToma.loc[intersect_genes].T
IL4_PEC_PToma_inter["condition"] = 'TRM'

SS_PEC_NToma_inter = SS_PEC_NToma.loc[intersect_genes].T
SS_PEC_NToma_inter["condition"] = 'TRM'
SS_PEC_PToma_inter = SS_PEC_PToma.loc[intersect_genes].T
SS_PEC_PToma_inter["condition"] = 'TRM'

ThioIL4PCD206_PEC_PToma_inter = ThioIL4PCD206_PEC_PToma.loc[intersect_genes].T
ThioIL4PCD206_PEC_PToma_inter["condition"] = 'M\u03C6mo'
ThioPCD206_PEC_PToma_inter = ThioPCD206_PEC_PToma.loc[intersect_genes].T
ThioPCD206_PEC_PToma_inter["condition"] = 'M\u03C6mo'

data_inter_all = pd.concat([HP_PEC_NToma_inter, HP_PEC_PToma_inter, IL4_PEC_NToma_inter, IL4_PEC_PToma_inter, SS_PEC_NToma_inter, SS_PEC_PToma_inter, ThioIL4PCD206_PEC_PToma_inter, ThioPCD206_PEC_PToma_inter])

X = data_inter_all.loc[:,data_inter_all.columns[:-1]]
y = data_inter_all['condition']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 0)


rfc = RandomForestClassifier(n_estimators=200, random_state = 0, class_weight={'TRM': 1, 'M\u03C6mo': 4})
rfc.fit(X_train, y_train)
y_pred = rfc.predict(X_test)
print('Model accuracy score with 200 decision-trees : {0:0.4f}'. format(accuracy_score(y_test, y_pred)))

feature_scores = pd.Series(rfc.feature_importances_, index=X_train.columns).sort_values(ascending=False)
feature_scores_top40 = feature_scores.loc[feature_scores.index[:30]]

f, ax = plt.subplots(figsize=(5, 8))
ax = sns.barplot(x=feature_scores_top40, y=feature_scores_top40.index)
ax.set_yticklabels(feature_scores_top40.index)
ax.set_xlabel("Gini importance score")
ax.set_ylabel("Top 30 genes")
plt.show()

data_inter_all_expressed = pd.DataFrame(np.array(data_inter_all[feature_scores_top40.index.tolist()] >0, dtype=int), 
                                        columns=feature_scores_top40.index.tolist(), index=data_inter_all.index)
data_inter_all_expressed['condition'] = data_inter_all["condition"]
data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)

data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)
high_exp_gene = np.array(feature_scores_top40.index)[(data_inter_all[feature_scores_top40.index.tolist() + ["condition"]].groupby("condition").mean(0)).mean() >= 1.5]
data_inter_all_expressed_percent = data_inter_all_expressed_percent[high_exp_gene]

data_inter_dotplot = pd.melt(data_inter_all_expressed_percent*100, var_name='Top 40 genes', value_name="Percent Expressed")
data_inter_dotplot['Conditions'] = np.array(data_inter_all_expressed_percent.index.tolist()*len(high_exp_gene))
data_inter_dotplot['Average Expression'] = pd.melt(data_inter_all[list(high_exp_gene) + ["condition"]].groupby("condition").mean(0))['value']
data_inter_dotplot["Average Expression"] = [5 if AE >5 else AE for AE in data_inter_dotplot["Average Expression"]]
data_inter_dotplot["Average Expression"] = [1 if AE <1 else AE for AE in data_inter_dotplot["Average Expression"]]

y_axis = ['TRM', 'M\u03C6mo']
range_ = ['#990000', '#CC0000', '#FF0000', '#FF3333', '#FF6666', '#FF9999', '#FFCCCC']
alt.Chart(data_inter_dotplot).mark_circle().encode(
    x=alt.X('Top 40 genes:O', sort=high_exp_gene),
    y = alt.Y('Conditions:O', sort=y_axis),
    size = 'Percent Expressed:Q',
    color = alt.Color(
            'Average Expression:Q',
            scale = alt.Scale(domain= [max(data_inter_dotplot['Average Expression']),min(data_inter_dotplot['Average Expression'])], range = range_))
)

# 2. 
intersect_genes = list(set(HP_PEC_NToma.index) & set(HP_PEC_PToma.index) & set(IL4_PEC_NToma.index) & set(IL4_PEC_PToma.index) & 
                       set(SS_PEC_NToma.index) & set(SS_PEC_PToma.index))

HP_PEC_NToma_inter = HP_PEC_NToma.loc[intersect_genes].T
HP_PEC_NToma_inter["condition"] = 'TRM_em'
HP_PEC_PToma_inter = HP_PEC_PToma.loc[intersect_genes].T
HP_PEC_PToma_inter["condition"] = 'TRM_mo'

IL4_PEC_NToma_inter = IL4_PEC_NToma.loc[intersect_genes].T
IL4_PEC_NToma_inter["condition"] = 'TRM_em'
IL4_PEC_PToma_inter = IL4_PEC_PToma.loc[intersect_genes].T
IL4_PEC_PToma_inter["condition"] = 'TRM_mo'

SS_PEC_NToma_inter = SS_PEC_NToma.loc[intersect_genes].T
SS_PEC_NToma_inter["condition"] = 'TRM_em'
SS_PEC_PToma_inter = SS_PEC_PToma.loc[intersect_genes].T
SS_PEC_PToma_inter["condition"] = 'TRM_mo'


data_inter_all = pd.concat([HP_PEC_NToma_inter, HP_PEC_PToma_inter, IL4_PEC_NToma_inter, IL4_PEC_PToma_inter, SS_PEC_NToma_inter, SS_PEC_PToma_inter])
data_inter_all['condition'].value_counts()/len(data_inter_all)

X = data_inter_all.loc[:,data_inter_all.columns[:-1]]
y = data_inter_all['condition']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 0)
rfc = RandomForestClassifier(n_estimators=200, random_state=0, class_weight={'TRM_em': 1, 'TRM_mo': 2})
rfc.fit(X_train, y_train)
y_pred = rfc.predict(X_test)
print('Model accuracy score with 200 decision-trees : {0:0.4f}'. format(accuracy_score(y_test, y_pred)))

feature_scores = pd.Series(rfc.feature_importances_, index=X_train.columns).sort_values(ascending=False)
feature_scores_top40 = feature_scores.loc[feature_scores.index[:30]]

f, ax = plt.subplots(figsize=(5, 8))
ax = sns.barplot(x=feature_scores_top40, y=feature_scores_top40.index)
ax.set_yticklabels(feature_scores_top40.index)
ax.set_xlabel("Gini importance score")
ax.set_ylabel("Top 30 genes")
plt.show()

data_inter_all_expressed = pd.DataFrame(np.array(data_inter_all[feature_scores_top40.index.tolist()] >0, dtype=int), 
                                        columns=feature_scores_top40.index.tolist(), index=data_inter_all.index)
data_inter_all_expressed['condition'] = data_inter_all["condition"]
data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)
data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)
high_exp_gene = np.array(feature_scores_top40.index)[(data_inter_all[feature_scores_top40.index.tolist() + ["condition"]].groupby("condition").mean(0)).mean() >= 1.5]
data_inter_all_expressed_percent = data_inter_all_expressed_percent[high_exp_gene]

data_inter_dotplot = pd.melt(data_inter_all_expressed_percent*100, var_name='Top 40 genes', value_name="Percent Expressed")
data_inter_dotplot['Conditions'] = np.array(data_inter_all_expressed_percent.index.tolist()*len(high_exp_gene))
data_inter_dotplot['Average Expression'] = pd.melt(data_inter_all[feature_scores_top40.index.tolist() + ["condition"]].groupby("condition").mean(0))['value']
data_inter_dotplot['Average Expression'] = pd.melt(data_inter_all[list(high_exp_gene) + ["condition"]].groupby("condition").mean(0))['value']
data_inter_dotplot["Average Expression"] = [5 if AE >5 else AE for AE in data_inter_dotplot["Average Expression"]]
data_inter_dotplot["Average Expression"] = [1 if AE <1 else AE for AE in data_inter_dotplot["Average Expression"]]

y_axis = ["HP_PEC_NToma",  "IL4_PEC_NToma", "SS_PEC_NToma", "HP_PEC_PToma","IL4_PEC_PToma", "SS_PEC_PToma"]
range_ = ['#990000', '#CC0000', '#FF0000', '#FF3333', '#FF6666', '#FF9999', '#FFCCCC']
alt.Chart(data_inter_dotplot).mark_circle().encode(
    x=alt.X('Top 40 genes:O', sort=high_exp_gene),
    y = alt.Y('Conditions:O', sort=y_axis),
    size = 'Percent Expressed:Q',
    color = alt.Color(
            'Average Expression:Q',
            scale = alt.Scale(domain= [max(data_inter_dotplot['Average Expression']),min(data_inter_dotplot['Average Expression'])], range = range_))
)

# 3. 
intersect_genes = list(set(HP_PEC_NToma.index) & set(HP_PEC_PToma.index) & set(IL4_PEC_NToma.index) & set(IL4_PEC_PToma.index))

HP_PEC_NToma_inter = HP_PEC_NToma.loc[intersect_genes].T
HP_PEC_NToma_inter["condition"] = 'HP'
HP_PEC_PToma_inter = HP_PEC_PToma.loc[intersect_genes].T
HP_PEC_PToma_inter["condition"] = 'HP'

IL4_PEC_NToma_inter = IL4_PEC_NToma.loc[intersect_genes].T
IL4_PEC_NToma_inter["condition"] = 'IL4'
IL4_PEC_PToma_inter = IL4_PEC_PToma.loc[intersect_genes].T
IL4_PEC_PToma_inter["condition"] = 'IL4'

data_inter_all = pd.concat([HP_PEC_NToma_inter, HP_PEC_PToma_inter, IL4_PEC_NToma_inter, IL4_PEC_PToma_inter])

X = data_inter_all.loc[:,data_inter_all.columns[:-1]]
y = data_inter_all['condition']


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.3, random_state = 0)
rfc = RandomForestClassifier(n_estimators=200, random_state=0)
rfc.fit(X_train, y_train)
y_pred = rfc.predict(X_test)
print('Model accuracy score with 200 decision-trees : {0:0.4f}'. format(accuracy_score(y_test, y_pred)))

feature_scores = pd.Series(rfc.feature_importances_, index=X_train.columns).sort_values(ascending=False)
feature_scores_top40 = feature_scores.loc[feature_scores.index[:30]]
feature_scores

f, ax = plt.subplots(figsize=(5, 8))
ax = sns.barplot(x=feature_scores_top40, y=feature_scores_top40.index)
ax.set_yticklabels(feature_scores_top40.index)
ax.set_xlabel("Gini importance score")
ax.set_ylabel("Top 30 genes")
plt.show()

data_inter_all_expressed = pd.DataFrame(np.array(data_inter_all[feature_scores_top40.index.tolist()] >0, dtype=int), 
                                        columns=feature_scores_top40.index.tolist(), index=data_inter_all.index)
data_inter_all_expressed['condition'] = data_inter_all["condition"]
data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)

data_inter_all_expressed_percent = pd.DataFrame((np.array(data_inter_all_expressed.groupby("condition").sum(0)).T/np.array(data_inter_all['condition'].value_counts().loc[data_inter_all_expressed.groupby("condition").sum(0).index])).T,
                                                index = data_inter_all_expressed.groupby("condition").sum(0).index, columns=data_inter_all_expressed.groupby("condition").sum(0).columns)

high_exp_gene = np.array(feature_scores_top40.index)[(data_inter_all[feature_scores_top40.index.tolist() + ["condition"]].groupby("condition").mean(0)).mean() >= 1.5]
data_inter_all_expressed_percent = data_inter_all_expressed_percent[high_exp_gene]

data_inter_dotplot = pd.melt(data_inter_all_expressed_percent*100, var_name='Top 40 genes', value_name="Percent Expressed")
data_inter_dotplot['Conditions'] = np.array(data_inter_all_expressed_percent.index.tolist()*len(high_exp_gene))
data_inter_dotplot['Average Expression'] = pd.melt(data_inter_all[list(high_exp_gene) + ["condition"]].groupby("condition").mean(0))['value']
data_inter_dotplot["Average Expression"] = [5 if AE >5 else AE for AE in data_inter_dotplot["Average Expression"]]
data_inter_dotplot["Average Expression"] = [1 if AE <1 else AE for AE in data_inter_dotplot["Average Expression"]]


y_axis = [ "IL4_PEC_NToma","IL4_PEC_PToma", "HP_PEC_NToma","HP_PEC_PToma"]
range_ = ['#990000', '#CC0000', '#FF0000', '#FF3333', '#FF6666', '#FF9999', '#FFCCCC']
alt.Chart(data_inter_dotplot).mark_circle().encode(
    x=alt.X('Top 40 genes:O', sort=high_exp_gene),
    y = alt.Y('Conditions:O', sort=y_axis),
    size = 'Percent Expressed:Q',
    color = alt.Color(
            'Average Expression:Q',
            scale = alt.Scale(domain= [max(data_inter_dotplot['Average Expression']),min(data_inter_dotplot['Average Expression'])], range = range_))
)

