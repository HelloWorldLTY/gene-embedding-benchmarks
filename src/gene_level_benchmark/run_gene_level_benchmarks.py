import helper
import pickle


if __name__ == "__main__":
    (
        reindexed_embeddings,
        gene_lists,
        reference_node2index,
        reference_genes,
    ) = helper.get_embeddings()

    # OMIM
    (
        doid_prop_use,
        doid_slim_use,
        doid_to_slim,
        graph,
        filtered_doid_to_slim,
    ) = helper.load_annotations(
        reference_genes,
        gmt_direct="data/gmt/omim.20231030.direct.gmt",
        gmt_prop="data/gmt/omim.20231030.prop.gmt",
        slim="data/slim_sets/doid_agr_slim.txt",
        obo="data/obo/doid.obo",
    )
    holdout_dict, cv_fold1_dict, cv_fold2_dict, cv_fold3_dict = helper.fold_split(
        doid_prop_use, doid_slim_use, doid_to_slim
    )

    file_names = [
        "bin/omim_cv_fold1_dict.pkl",
        "bin/omim_cv_fold2_dict.pkl",
        "bin/omim_cv_fold3_dict.pkl",
        "bin/omim_holdout_dict.pkl",
    ]
    data_dicts = [cv_fold1_dict, cv_fold2_dict, cv_fold3_dict, holdout_dict]

    for file_name, data_dict in zip(file_names, data_dicts):
        with open(file_name, "wb") as f:
            pickle.dump(data_dict, f)

    file_names

    all_fold_results, all_holdout_results = helper.run_benchmark(
        reindexed_embeddings, cv_fold1_dict, cv_fold2_dict, cv_fold3_dict, holdout_dict
    )

    with open("bin/omim_all_fold_results.pkl", "wb") as f:
        pickle.dump(all_fold_results, f)

    with open("bin/omim_all_holdout_results.pkl", "wb") as f:
        pickle.dump(all_holdout_results, f)

    (
        df_summary_with_meta,
        fold_auc_df,
        fold_auprc_df,
        fold_pr10_df,
        holdout_auc_df,
        holdout_auprc_df,
        holdout_pr10_df,
        holdout_time_df,
    ) = helper.reformat_results(all_fold_results, all_holdout_results)

    helper.plot_scatter(holdout_auc_df)
    helper.anova(df_summary_with_meta)

    g = helper.plot_slim_clustermap(
        df=fold_auc_df.loc[fold_auc_df.mean(axis=1).sort_values(ascending=False).index],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="OMIM Fold AUC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=fold_auprc_df.loc[
            fold_auprc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="OMIM Fold AUPRC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=fold_pr10_df.loc[
            fold_pr10_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="OMIM Fold PR@10",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_auc_df.loc[
            holdout_auc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="OMIM Holdout AUC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_auprc_df.loc[
            holdout_auprc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="OMIM Holdout AUPRC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_pr10_df.loc[
            holdout_pr10_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="OMIM Holdout PR@10",
        legend=True,
    )

    # GO All
    (
        go_prop_use,
        go_slim_use,
        go_to_slim,
        graph,
        filtered_go_to_slim,
    ) = helper.load_annotations(
        reference_genes,
        gmt_direct="data/gmt/hsa_EXP_ALL_BP_direct.gmt",
        gmt_prop="data/gmt/hsa_EXP_ALL_BP_propagated.gmt",
        slim="data/slim_sets/goslim_agr.tsv",
        obo="data/obo/go.obo",
    )
    holdout_dict, cv_fold1_dict, cv_fold2_dict, cv_fold3_dict = helper.fold_split(
        doid_prop_use, doid_slim_use, doid_to_slim
    )

    file_names = [
        "bin/go_cv_fold1_dict.pkl",
        "bin/go_cv_fold2_dict.pkl",
        "bin/go_cv_fold3_dict.pkl",
        "bin/go_holdout_dict.pkl",
    ]
    data_dicts = [cv_fold1_dict, cv_fold2_dict, cv_fold3_dict, holdout_dict]

    for file_name, data_dict in zip(file_names, data_dicts):
        with open(file_name, "wb") as f:
            pickle.dump(data_dict, f)

    file_names

    all_fold_results, all_holdout_results = helper.run_benchmark(
        reindexed_embeddings, cv_fold1_dict, cv_fold2_dict, cv_fold3_dict, holdout_dict
    )

    with open("bin/go_all_fold_results.pkl", "wb") as f:
        pickle.dump(all_fold_results, f)

    with open("bin/go_all_holdout_results.pkl", "wb") as f:
        pickle.dump(all_holdout_results, f)

    (
        df_summary_with_meta,
        fold_auc_df,
        fold_auprc_df,
        fold_pr10_df,
        holdout_auc_df,
        holdout_auprc_df,
        holdout_pr10_df,
        holdout_time_df,
    ) = helper.reformat_results(all_fold_results, all_holdout_results)

    helper.plot_scatter(holdout_auc_df)
    helper.anova(df_summary_with_meta)

    g = helper.plot_slim_clustermap(
        df=fold_auc_df.loc[fold_auc_df.mean(axis=1).sort_values(ascending=False).index],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO Fold AUC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=fold_auprc_df.loc[
            fold_auprc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO Fold AUPRC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=fold_pr10_df.loc[
            fold_pr10_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO Fold PR@10",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_auc_df.loc[
            holdout_auc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO Holdout AUC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_auprc_df.loc[
            holdout_auprc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO Holdout AUPRC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_pr10_df.loc[
            holdout_pr10_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO Holdout PR@10",
        legend=True,
    )

    # GO After 2020
    (
        go_prop_use,
        go_slim_use,
        go_to_slim,
        graph,
        filtered_go_to_slim,
    ) = helper.load_annotations(
        reference_genes,
        gmt_direct='data/gmt/2020_hsa_ALL_BP_direct.gmt"',
        gmt_prop='data/gmt/2020_hsa_ALL_BP_propagated.gmt"',
        slim="data/slim_sets/goslim_agr.tsv",
        obo="data/obo/go.obo",
    )
    holdout_dict, cv_fold1_dict, cv_fold2_dict, cv_fold3_dict = helper.fold_split(
        doid_prop_use, doid_slim_use, doid_to_slim
    )

    file_names = [
        "bin/go2020_cv_fold1_dict.pkl",
        "bin/go2020_cv_fold2_dict.pkl",
        "bin/go2020_cv_fold3_dict.pkl",
        "bin/go2020_holdout_dict.pkl",
    ]
    data_dicts = [cv_fold1_dict, cv_fold2_dict, cv_fold3_dict, holdout_dict]

    for file_name, data_dict in zip(file_names, data_dicts):
        with open(file_name, "wb") as f:
            pickle.dump(data_dict, f)

    file_names

    all_fold_results, all_holdout_results = helper.run_benchmark(
        reindexed_embeddings, cv_fold1_dict, cv_fold2_dict, cv_fold3_dict, holdout_dict
    )

    with open("bin/go2020_all_fold_results.pkl", "wb") as f:
        pickle.dump(all_fold_results, f)

    with open("bin/go2020_all_holdout_results.pkl", "wb") as f:
        pickle.dump(all_holdout_results, f)

    (
        df_summary_with_meta,
        fold_auc_df,
        fold_auprc_df,
        fold_pr10_df,
        holdout_auc_df,
        holdout_auprc_df,
        holdout_pr10_df,
        holdout_time_df,
    ) = helper.reformat_results(all_fold_results, all_holdout_results)

    helper.plot_scatter(holdout_auc_df)
    helper.anova(df_summary_with_meta)

    g = helper.plot_slim_clustermap(
        df=fold_auc_df.loc[fold_auc_df.mean(axis=1).sort_values(ascending=False).index],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO 2020 Fold AUC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=fold_auprc_df.loc[
            fold_auprc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO 2020 Fold AUPRC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=fold_pr10_df.loc[
            fold_pr10_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO 2020 Fold PR@10",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_auc_df.loc[
            holdout_auc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO 2020 Holdout AUC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_auprc_df.loc[
            holdout_auprc_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO 2020 Holdout AUPRC",
        legend=True,
    )
    g = helper.plot_slim_clustermap(
        df=holdout_pr10_df.loc[
            holdout_pr10_df.mean(axis=1).sort_values(ascending=False).index
        ],
        graph=graph,
        filtered_annot_to_slim=filtered_doid_to_slim,
        title="GO 2020 Holdout PR@10",
        legend=True,
    )
