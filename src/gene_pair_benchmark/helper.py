import os
import glob
import pandas as pd
import random
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import roc_auc_score, average_precision_score, precision_score
from sklearn.svm import SVC
import numpy as np
import time


def load_embeddings(folder_path="data/embeddings/intersect"):
    subfolders = [f.path for f in os.scandir(folder_path) if f.is_dir()]
    # dict of embeddings and gene list
    embeddings = {}
    gene_lists = {}

    for subfolder in subfolders:
        print(f"Processing subfolder: {subfolder}")
        gene_txt_files = glob.glob(os.path.join(subfolder, "*.txt"))
        if not gene_txt_files:
            print(f"No txt file found in {subfolder}")
            continue
        gene_file = gene_txt_files[0]
        with open(gene_file, "r") as f:
            genes = [line.strip() for line in f]
        gene_lists[subfolder] = genes

        csv_files = glob.glob(os.path.join(subfolder, "*.csv"))
        if not csv_files:
            print(f"No csv file found in {subfolder}")
            continue
        csv_file = csv_files[0]
        embedding = pd.read_csv(csv_file, header=None)
        embeddings[subfolder] = embedding

    reference_subfolder = list(gene_lists.keys())[0]
    reference_genes = gene_lists[
        reference_subfolder
    ]  # get an order of genes that the other embedigns will now follow
    reference_node2index = {x: i for i, x in enumerate(reference_genes)}

    for subfolder in embeddings:
        embeddings[subfolder] = embeddings[subfolder].values

    return embeddings, reference_node2index


def get_data(pairs, emb, mapping):
    mapped_pairs = [(mapping[x], mapping[y]) for x, y in pairs]
    data = [emb[x] + emb[y] for x, y in mapped_pairs]
    return data


def fold_split(used_nodes, num_folds=3, holdout_fraction=0.2):
    random.seed(42)
    holdout_size = int(len(used_nodes) * holdout_fraction)

    random.shuffle(used_nodes)

    holdout_nodes = used_nodes[:holdout_size]
    cv_nodes = used_nodes[holdout_size:]

    split_size = len(cv_nodes) // num_folds
    splits = [
        cv_nodes[:split_size],
        cv_nodes[split_size : 2 * split_size],
        cv_nodes[2 * split_size :],
    ]
    return holdout_nodes, splits, cv_nodes


def setup_data(positive_pairs, splits, holdout_nodes, cv_nodes, num_folds=3):
    random.seed(42)

    # create the data for models to run on/test on
    fold_splits = {}
    for i in range(num_folds):
        test_nodes = splits[i]
        train_nodes = []
        for j in range(num_folds):
            if j != i:
                train_nodes += splits[j]

        train_positive_pairs = [
            (x, y) for x, y in positive_pairs if x in train_nodes and y in train_nodes
        ]
        test_positive_pairs = [
            (x, y) for x, y in positive_pairs if x in test_nodes and y in test_nodes
        ]

        train_all_pairs = [(x, y) for x in train_nodes for y in train_nodes if x < y]
        test_all_pairs = [(x, y) for x in test_nodes for y in test_nodes if x < y]

        train_negative_pairs = list(set(train_all_pairs) - set(train_positive_pairs))
        test_negative_pairs = list(set(test_all_pairs) - set(test_positive_pairs))

        random.shuffle(train_negative_pairs)
        random.shuffle(test_negative_pairs)
        train_negative_pairs = train_negative_pairs[: len(train_positive_pairs) * 10]
        test_negative_pairs = test_negative_pairs[: len(test_positive_pairs) * 10]

        fold_splits[i] = {
            "train_pairs": train_positive_pairs + train_negative_pairs,
            "test_pairs": test_positive_pairs + test_negative_pairs,
            "train_labels": [1] * len(train_positive_pairs)
            + [0] * len(train_negative_pairs),
            "test_labels": [1] * len(test_positive_pairs)
            + [0] * len(test_negative_pairs),
        }

    holdout_positive_pairs = [
        (x, y) for x, y in positive_pairs if x in holdout_nodes and y in holdout_nodes
    ]
    holdout_all_pairs = [(x, y) for x in holdout_nodes for y in holdout_nodes if x < y]
    holdout_negative_pairs = list(set(holdout_all_pairs) - set(holdout_positive_pairs))
    random.shuffle(holdout_negative_pairs)
    holdout_negative_pairs = holdout_negative_pairs[: len(holdout_positive_pairs) * 10]

    holdout_pairs = holdout_positive_pairs + holdout_negative_pairs
    holdout_labels = [1] * len(holdout_positive_pairs) + [0] * len(
        holdout_negative_pairs
    )

    final_train_nodes = cv_nodes
    final_train_positive_pairs = [
        (x, y)
        for x, y in positive_pairs
        if x in final_train_nodes and y in final_train_nodes
    ]
    final_train_all_pairs = [
        (x, y) for x in final_train_nodes for y in final_train_nodes if x < y
    ]
    final_train_negative_pairs = list(
        set(final_train_all_pairs) - set(final_train_positive_pairs)
    )
    random.shuffle(final_train_negative_pairs)
    final_train_negative_pairs = final_train_negative_pairs[
        : len(final_train_positive_pairs) * 10
    ]

    final_train_pairs = final_train_positive_pairs + final_train_negative_pairs
    final_train_labels = [1] * len(final_train_positive_pairs) + [0] * len(
        final_train_negative_pairs
    )

    return (
        fold_splits,
        holdout_pairs,
        holdout_labels,
        final_train_pairs,
        final_train_labels,
    )


def run_SVM(
    embeddings,
    reference_node2index,
    splits,
    fold_splits,
    holdout_pairs,
    holdout_labels,
    final_train_pairs,
    final_train_labels,
    num_folds=3,
    log_file="genepair_output_log.txt",
):
    fold_results_records = []
    holdout_results_records = []

    def log_print(msg, log_file=log_file):
        print(msg)
        with open(log_file, "a") as f:
            f.write(msg + "\n")

    def precision_at_k(y_true, y_scores, k=10):
        """Calculate Precision at K."""
        if len(y_scores) > k:
            top_k_indices = np.argsort(y_scores)[-k:][::-1]
        else:
            top_k_indices = np.argsort(y_scores)[::-1]
        top_k_true = np.array(y_true)[top_k_indices]
        return np.mean(top_k_true)

    C_values = [0.1, 1, 10, 100, 1000]

    for subfolder in embeddings:
        emb = embeddings[subfolder]
        log_print(f"Evaluating embedding from subfolder: {subfolder}")

        c_performance = []

        for C in C_values:
            fold_results = []
            fold_times = []
            log_print(f"Testing C={C}")
            for i in range(num_folds):
                fold_start_time = time.time()

                test_nodes = splits[i]
                train_nodes = []
                for j in range(num_folds):
                    if j != i:
                        train_nodes += splits[j]

                fold_data = fold_splits[i]
                train_pairs = fold_data["train_pairs"]
                test_pairs = fold_data["test_pairs"]
                train_label = fold_data["train_labels"]
                test_label = fold_data["test_labels"]

                train_data = get_data(train_pairs, emb, reference_node2index)
                test_data = get_data(test_pairs, emb, reference_node2index)

                clf = SVC(class_weight="balanced", C=C)
                clf.fit(train_data, train_label)
                probabilities = clf.decision_function(test_data)

                auc = roc_auc_score(test_label, probabilities)
                auprc = average_precision_score(test_label, probabilities)
                pr10 = precision_at_k(test_label, probabilities, k=10)

                fold_end_time = time.time()
                fold_duration = fold_end_time - fold_start_time
                fold_times.append(fold_duration)

                log_print(
                    f"Fold {i + 1}: AUC: {auc:.4f}, AUPRC: {auprc:.4f}, PR@10: {pr10:.4f}, Time: {fold_duration:.2f}s"
                )
                fold_results.append((auc, auprc, pr10))

            avg_auc = np.mean([r[0] for r in fold_results])
            avg_auprc = np.mean([r[1] for r in fold_results])
            avg_pr10 = np.mean([r[2] for r in fold_results])
            avg_fold_time = np.mean(fold_times)
            log_print(
                f"Cross-val average for C={C}: AUC={avg_auc:.4f}, AUPRC={avg_auprc:.4f}, PR@10={avg_pr10:.4f}, Avg Time/Fold: {avg_fold_time:.2f}s"
            )

            c_performance.append(
                {
                    "C": C,
                    "AUC": avg_auc,
                    "AUPRC": avg_auprc,
                    "PR@10": avg_pr10,
                    "avg_fold_time": avg_fold_time,
                }
            )

        best_C_entry = max(c_performance, key=lambda x: x["AUC"])
        best_C = best_C_entry["C"]
        best_avg_auc = best_C_entry["AUC"]
        best_avg_auprc = best_C_entry["AUPRC"]
        best_avg_pr10 = best_C_entry["PR@10"]
        best_avg_fold_time = best_C_entry["avg_fold_time"]

        log_print(
            f"Best C for {subfolder}: C={best_C} with AUC={best_avg_auc:.4f}, AUPRC={best_avg_auprc:.4f}, PR@10={best_avg_pr10:.4f}, Avg Time/Fold: {best_avg_fold_time:.2f}s"
        )

        fold_results_records.append(
            {
                "subfolder": subfolder,
                "C": best_C,
                "AUC": best_avg_auc,
                "AUPRC": best_avg_auprc,
                "PR@10": best_avg_pr10,
                "avg_fold_time": best_avg_fold_time,
            }
        )

        holdout_start_time = time.time()

        final_train_data = get_data(final_train_pairs, emb, reference_node2index)
        final_clf = SVC(class_weight="balanced", C=best_C)
        final_clf.fit(final_train_data, final_train_labels)

        holdout_data = get_data(holdout_pairs, emb, reference_node2index)
        holdout_prob = final_clf.decision_function(holdout_data)

        holdout_auc = roc_auc_score(holdout_labels, holdout_prob)
        holdout_auprc = average_precision_score(holdout_labels, holdout_prob)
        holdout_pr10 = precision_at_k(holdout_labels, holdout_prob, k=10)

        holdout_end_time = time.time()
        holdout_duration = holdout_end_time - holdout_start_time

        log_print(f"Holdout performance for {subfolder} with best C={best_C}:")
        log_print(f"Holdout AUC: {holdout_auc:.4f}")
        log_print(f"Holdout AUPRC: {holdout_auprc:.4f}")
        log_print(f"Holdout PR@10: {holdout_pr10:.4f}")
        log_print(f"Holdout Duration: {holdout_duration:.2f}s")

        holdout_results_records.append(
            {
                "subfolder": subfolder,
                "C": best_C,
                "AUC": holdout_auc,
                "AUPRC": holdout_auprc,
                "PR@10": holdout_pr10,
                "holdout_time": holdout_duration,
            }
        )

    fold_results_df = pd.DataFrame(fold_results_records).set_index("subfolder")
    holdout_results_df = pd.DataFrame(holdout_results_records).set_index("subfolder")

    sorted_fold_results_df = fold_results_df.sort_values(by="AUPRC", ascending=False)
    sorted_holdout_results_df = holdout_results_df.sort_values(
        by="AUPRC", ascending=False
    )
    sorted_fold_results_df["dim"] = sorted_fold_results_df.index.map(
        lambda sf: embeddings[sf].shape[1]
    )
    sorted_holdout_results_df["dim"] = sorted_holdout_results_df.index.map(
        lambda sf: embeddings[sf].shape[1]
    )

    return sorted_fold_results_df, sorted_holdout_results_df
