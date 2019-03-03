# coding=utf-8
# Load libraries

import networkx as nx
import pandas as pd
import statsmodels.stats.multitest
from scipy.stats import hypergeom
import os
import sys
import copy
import functools
import matplotlib.pyplot as plt


def Geom_p_value(ourSet, n_refSet, targetSet):
    m = n_refSet
    n = len(ourSet)
    interSecTar = ourSet.intersection(targetSet)
    k = len(interSecTar)
    # if k < 1:
    #    return(print('sets is disjoint'))
    l = len(targetSet)
    k_expect = l * n / m
    if k_expect > k:
        p_value = 0
        for i in range(0, k):
            p_value -= hypergeom.pmf(i, m, l, n)
    else:
        p_value = 0
        for i in range(k, n):
            p_value += hypergeom.pmf(i, m, l, n)
    if p_value == 0:
        return (1, 0)
    else:
        return (p_value, interSecTar)

def save_graph(nodes, path):
    color_map = []
    result_graph = G.subgraph(nodes)
    for node in result_graph:
        if node in compulsory_genes:
            color_map.append('red')
        else:
            if node in additional_genes:
                color_map.append('green')
            else:
                color_map.append('blue')
    plt.figure(figsize=(14, 14))
    nx.draw(result_graph, node_color=color_map, with_labels=True)
    plt.savefig(path)


PATH_TO_STRINGDB = "../data/combined_filtred_0.8.txt"
#PATH_TO_DISGENET = "../data/all_ather_0.08gda.tsv"
PATH_TO_DISGENET = sys.argv[1]

#PATH_TO_ORA_RESULTS = "../data/Result_ORA.csv"

#PATH_TO_DEG_FILE = "../data/degs_lgFC_1.5_with_ENSP.csv"
PATH_TO_DEG_FILE = sys.argv[2]

# Load STRING data

df = pd.read_csv(PATH_TO_STRINGDB, sep='\t')

G = nx.convert_matrix.from_pandas_edgelist(df, df.columns[0], df.columns[1])

# subset conencted component
# взять гены, уложенные в максимальной компоненте связности

LCC_vertex_set = max(nx.connected_components(G), key=len)


# по этим генам построить подграф, являющийся компонентой связности

LCC_subgraph = G.subgraph(LCC_vertex_set)


# # load data from DEG or DisGeNET

DEG_set = pd.read_csv(PATH_TO_DEG_FILE, usecols=[7])['hgnc_symbol']
DisGeNET = set(PATH_TO_DISGENET.strip().split('\r\n'))
#DisGeNET = pd.read_csv(PATH_TO_DISGENET, sep='\t', usecols=[3])['symbol']

# # subset by DisGeNET genes
# 

# создать набор генов, которые общие для выделенной компоненты связности и генов с DisGeNET

subsetted_graph_vs = set(LCC_subgraph.nodes) & set(DisGeNET)

# связность не проверяем - она низкая до идиотизма (а кстати, это идея - посмотреть; авось у нас
# она наоборот, слишком низкая?)

# затаскиваем результаты генного обогащенич

KEGGdb = dict()
KEGGf = open('../data/reactome.v6.1.symbols.gmt.txt', 'r') #reactome из MsigDB
for line in KEGGf:
    preList = line.strip().split('\t')
    KEGGdb[preList[0]] = preList[2:]

ORA_results = pd.DataFrame(columns=['p_value','interSecTar','Path','PathMembr'])

count = 0
for i in list(KEGGdb.keys()):
    p_value, interSecTar = Geom_p_value(subsetted_graph_vs, 20000, KEGGdb[i])
    if interSecTar != 0: interSecTar = len(interSecTar)
    ORA_results.loc[count] = [p_value, interSecTar, i, "/".join(KEGGdb[i])]
    count += 1

# парсим данные p-values с помощью поправки на множественные сравнения

passed_index = ORA_results.index[statsmodels.stats.multitest.multipletests(ORA_results['p_value'], 
                                                                           method='fdr_bh')[0]]




# извлекаем гены, которыми мы обогащаем сеть

network_enriching_genes = functools.reduce(lambda x, y: x + y, \
                 ORA_results.loc[passed_index, :]['PathMembr'].apply(lambda x:x.split('/')).values)


# избавляемся от дубликатов

network_enriching_genes = list(set(network_enriching_genes))


# добавляем ближайших соседей - они по идее должны сделать гены, лежащие в перечне внутри компоненты
# связности, связным графом.

LCC_complementing_genes = list()

for i in subsetted_graph_vs:
    LCC_complementing_genes.extend(list(LCC_subgraph.neighbors(i)))

# теперь объединяем все свои сеты. Это должно быть пересечение обогащающих и дополнящих генов +
# объединение с предполагаемыми генами, ассоциированными с атеросклерозом.

updated_vertices_set = ((set(LCC_complementing_genes) & 
                        set(network_enriching_genes)) | 
                        subsetted_graph_vs)


# На полученном перечне генов строим новый индуцированный подграф.
# Он должен иметь большую по размеру компоненту связности, чем просто subsetted_graph_vs. И большую,
# чем нечто другое(!)

# len(max(nx.connected_components(LCC_subgraph.subgraph(updated_vertices_set)), key=len))


# Теперь валидация. Надо итеративно ходить по графу
# На каждой итерации надо находить элемент, добавленный во время обогащения, с наименьшей
# центральностью, и его выкидывать.
# Если при этом граф не разваливается на несколько компонент связности - продолжать работу.
# Если граф разваивается - оставить ту компоненту, в которой сохранены все гены, взятые изначально.


ITERATE_FLAG = True

working_set = max(nx.connected_components(LCC_subgraph.subgraph(updated_vertices_set)), key=len)


compulsory_genes = set(subsetted_graph_vs) & working_set
additional_genes = (set(LCC_complementing_genes) & set(network_enriching_genes) &
                    set(working_set)).difference(compulsory_genes)

working_graph = LCC_subgraph.subgraph(working_set)

temp_additional_genes = copy.deepcopy(additional_genes)

while ITERATE_FLAG:

    temp_additional_genes = copy.deepcopy(additional_genes)
    
    # граф надо постоянно обновлять
    range_list = nx.betweenness_centrality_subset(G=working_graph.subgraph(additional_genes |
                                                                           compulsory_genes),
                                                  sources=additional_genes,
                                                  targets=compulsory_genes,)
    
    for gene in compulsory_genes:
        if gene in range_list.keys():
            range_list.pop(gene)
    
    
    for key, value in range_list.items():

    # если мера центральности равна нулю - сразу кидаем, без проверок
        if value == 0:

            temp_additional_genes.remove(key)

                # сразу сохраняем изменения
            additional_genes = copy.deepcopy(temp_additional_genes)


        # в противном случае сравниваем с минимумом
        elif value == min(range_list.values()):


            # если равна минимуму - выкидываем её
            if key not in compulsory_genes:
                temp_additional_genes.remove(key)

                # проверить, что при удалении минимального по значению узла компонента связности
                # не разваливается для это делаем перечень компонент связности
                temp_list = list(nx.connected_components(
                    working_graph.subgraph(temp_additional_genes | compulsory_genes)))
                # если компонента связности осталась одна - сохраняем изменения
                if len(temp_list) == 1:
                    additional_genes = copy.deepcopy(temp_additional_genes)

                # учесть вариант, когда компонента таки разваливается
                elif len(temp_list) > 1:

                    # создать список компонент с обязательными генами
                    lists_with_compulsory = list(filter(lambda x: len(x & compulsory_genes) > 0,
                                                        temp_list))

                    # если такой список состоит из 1 элемента - всё в порядке.
                    # Просто чистим рабочий список
                    if len(lists_with_compulsory) == 1:
                        additional_genes.union(set(lists_with_compulsory[0]))

                    # если список больше 1 - значит, компонента разваливается, и
                    # мы пришли к обязательным генам. не сохраняем
                    else:
                        ITERATE_FLAG = False
                        break

print("These are additional genes identified by algorithm:")
print(" ".join(additional_genes))

# # MARLEZONE BALETTE PART 2
balette_part_2_vertices = set(DEG_set).union(compulsory_genes).union(additional_genes).union(
    network_enriching_genes).intersection(G.nodes())

# len(balette_part_2_vertices)
# len(G.nodes())
# len(DEG_set)
# len(DEG_set)
genes_being_tested = set(DEG_set).difference(compulsory_genes).difference(
    additional_genes).difference(network_enriching_genes).intersection(G.nodes())
network_balette_part_2 = G.subgraph(balette_part_2_vertices)
# len(network_balette_part_2.nodes())
# len(genes_being_tested)
result_dict = {}

for gene in genes_being_tested:
    result_dict[gene] = {}
    min_length = 50000
    min_path = None
    for target in compulsory_genes.union(additional_genes):
        try:
            path = nx.bidirectional_shortest_path(network_balette_part_2, gene, target)
            if len(path) < min_length:
                result_dict[gene]['min_length'] = len(path)
                result_dict[gene]['min_path'] = path
                min_length = len(path)
                min_path = path
        except:
            pass


additional_genes.intersection(set(DEG_set))

work_dir = os.path.dirname(PATH_TO_DEG_FILE)
save_graph(compulsory_genes, work_dir + '/compulsory.png')
save_graph(working_graph, work_dir + '/all.png')
save_graph(compulsory_genes.union(additional_genes), work_dir + '/compulsory_additional.png')
