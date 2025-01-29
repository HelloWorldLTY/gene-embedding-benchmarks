import pandas as pd
import numpy as np
import os


def load_kegg_go_data(path, selected_row_terms, selected_column_terms):
    best_avg_matrices = []
    filename_order = []

    for file_name in os.listdir(path):
        if file_name.endswith(".csv"):
            filename_order.append(file_name.replace(".csv", ""))
            file_path = os.path.join(path, file_name)

            best_avg = pd.read_csv(file_path, delimiter=",", index_col=0)

            best_avg_matrix = best_avg.loc[selected_row_terms][
                selected_column_terms
            ].values

            best_avg_matrices.append(best_avg_matrix)

    return best_avg_matrices, filename_order


def generate_kegg_go_result(matrices, label_matrix):
    result_lists = [[] for i in range(len(matrices))]
    for i in range(label_matrix.shape[0]):
        for j, m in enumerate(matrices):
            label_indices = np.nonzero(label_matrix[i])

            order = matrices[j][i].argsort()
            ranks = order.argsort()
            result_lists[j].append(np.mean(ranks[list(label_indices)]))

    return result_lists


def generate_kegg_go_result2(matrices, label_matrix):
    label_matrix = label_matrix.T
    result_lists = [[] for i in range(len(matrices))]
    for i in range(label_matrix.shape[0]):
        for j, m in enumerate(matrices):
            label_indices = np.nonzero(label_matrix[i])

            order = matrices[j].T[i].argsort()
            ranks = order.argsort()
            result_lists[j].append(np.mean(ranks[list(label_indices)]))

    return result_lists
