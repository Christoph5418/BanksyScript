#!/usr/bin/env python
# coding: utf-8
import os


datasetID, sampleID, Nclusters = input().split()


def banksy(datasetID, sampleID, Nclusters):
    import sys
    print(sys.version)

    import os
    import numpy as np
    import pandas as pd
    import scipy.sparse as sparse

    import scanpy as sc
    from anndata import AnnData
    import anndata

    import matplotlib.pyplot as plt
    import seaborn as sns
    sc.logging.print_header()
    sc.set_figure_params(facecolor="white", figsize=(8, 8))
    sc.settings.verbosity = 3  # errors (0), warnings (1), info (2), hints (3)

    sns.set_style("white")

    np.random.seed(0)
    sampleFile = pd.read_csv('samples.csv')
    index = 0
    rowNum = -1
    headers = list(sampleFile.columns.values)
    sampleIDColumn = sampleFile['sample_id']
    for value in sampleIDColumn:
        if str(value) == sampleID:
            rowNum = index
        index = index + 1

    files_prefix = sampleFile.loc[rowNum, "files_prefix"]


    adata_filename = files_prefix + '_' + sampleID + '.h5ad'
    print(adata_filename)

    if os.path.isfile(adata_filename):

        adata = anndata.read_h5ad(adata_filename)

    else:
        counts = pd.read_csv(files_prefix + '_' + sampleID + '_UMImatrix.csv', header=None)
        coords = pd.read_csv(files_prefix + '_' + sampleID + '_coords.csv', header=None,
                             names=["xcoord", "ycoord"])

        counts = counts.transpose()

        sparse_X = sparse.csc_matrix(counts.values.T)

        adata = AnnData(X=sparse_X, obs=coords, var=pd.DataFrame(index=counts.index))

        adata.write(os.path(adata_filename))

    print(f"Matrix sparsity: {adata.X.nnz} filled elements "
          f"({adata.X.nnz / adata.X.shape[0] / adata.X.shape[1]:0.2f}) "
          f"out of {adata.X.shape[0] * adata.X.shape[1]}\n"
          f"max: {np.amax(adata.X.data)}, min: {np.amin(adata.X.data)}")

    raw_y, raw_x = adata.obs["xcoord"], adata.obs["ycoord"]
    adata.obsm["coord_xy"] = np.vstack((adata.obs["xcoord"].values, adata.obs["ycoord"].values)).T

    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 300], kde=False, bins=80, ax=axs[1])
    sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=80, ax=axs[2])
    sns.histplot(adata.obs["n_genes_by_counts"][
                     (adata.obs["n_genes_by_counts"] < 1000) & (adata.obs["n_genes_by_counts"] > 0)
                     ], kde=False, bins=100, ax=axs[3])
    fig.tight_layout()

    #THE FOLLOWING CODE WAS COMMENTED OUT BECAUSE IT FILTERS THE DATAPOINTS AND REDUCES ARRAY LENGTH SIGNFICANTLY

    # In[34]:

    # sc.pp.filter_cells(adata, min_counts=20)
    # sc.pp.filter_cells(adata, max_counts=1000)
    # print(f"#cells after count filter: {adata.n_obs}")

    # adata = adata[adata.obs["pct_counts_mt"] < 20]
    # print(f"#cells after MT filter: {adata.n_obs}")

    # sc.pp.filter_genes(adata, min_cells=10)
    # print(f"#genes after minimum cells per gene filter: {adata.n_vars}")

    # fig, axs = plt.subplots(1, 4, figsize=(15, 4))
    # sns.histplot(adata.obs["total_counts"], kde=False, ax=axs[0])
    # sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 2000], kde=False, bins=100, ax=axs[1])
    # sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axs[2])
    # sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 2000], kde=False, bins=100, ax=axs[3])
    # fig.tight_layout()

    # Filter out beads outside the puck
    # puck_center = (3330,3180) # (x,y)
    # puck_radius = 2550
    # puck_mask = np.sqrt((adata.obs["ycoord"] - puck_center[1])**2 + (adata.obs["xcoord"] - puck_center[0])**2) < puck_radius
    # adata = adata[puck_mask,:]
    # circle = plt.Circle(puck_center, puck_radius, color='g', fill=False)

    # fig, ax = plt.subplots(figsize = (8,8))
    # ax.scatter(raw_x, raw_y, s = 0.2, c = "red")
    # ax.scatter(adata.obs["xcoord"], adata.obs["ycoord"], s = 0.2, c = "slateblue")
    # ax.add_patch(circle)

    # adata

    # In[35]:

    # from skmisc.loess import loess
    # import skmisc.loess.loess as loess
    # import skmisc.loess

    def print_max_min(anndata: anndata.AnnData):
        print(f"\nMax: {anndata.X.max()}, Min: {anndata.X.min()}\n")

    n_top_genes = 2000
    filter_by_HVG = True
    # hvg_type = "seurat_v3"
    hvg_type = "seurat"

    # create a copy that will be log-transformed for old seurat version of HVG
    adata_log1p = adata.copy()

    # HVG detection on raw count data (seurat v3)
    # ----------------------------------------------

    # print("---- Raw data ----")
    # print_max_min(adata)
    # sc.pp.highly_variable_genes(adata,
    #                             flavor="seurat_v3",
    #                             n_top_genes=n_top_genes)

    # # normalize by total cell count
    sc.pp.normalize_total(adata, inplace=True)
    print("--- Normalized data -----")
    print_max_min(adata)

    # HVG detection on log-transformed data (old seurat)
    # --------------------------------------------------

    sc.pp.normalize_total(adata_log1p, inplace=True)
    sc.pp.log1p(adata_log1p)
    print("--- Normalized and log-transformed data -----")
    print_max_min(adata_log1p)

    sc.pp.highly_variable_genes(adata_log1p,
                                flavor="seurat",
                                n_top_genes=n_top_genes)

    # Compare
    # -------

    # print(adata_log1p.var["highly_variable"], )
    # print(f"\nOverlap between seurat method (on log transformed data)"
    #       f" and seurat v3 method (on raw count data):"
    #       f"{(adata.var['highly_variable'] & adata_log1p.var['highly_variable']).sum()}")

    # display(adata)
    # display(adata_log1p)

    # print("----- Top Genes, Seurat v3 -----")
    # display(adata.var[adata.var["highly_variable_rank"]<10])
    # print("----- Top Genes, Seurat -----")
    # display(adata_log1p.var[adata_log1p.var["highly_variable_rank"]<10])

    # Filter out non-highly variable genes
    # ------------------------------------

    if filter_by_HVG:

        if hvg_type == "seurat_v3":
            hvg_filter = adata.var["highly_variable"]
        elif hvg_type == "seurat":
            hvg_filter = adata_log1p.var["highly_variable"]
        else:
            raise ValueError("hvg_type not recognised")

        # save a copy of the unfiltered
        adata_allgenes = adata.copy()

        adata = adata[:, hvg_filter]

    # In[36]:

    # show annotation dataframes
    # --------------------------

    print("\n  Genes Info\n" + "=" * 20)
    print("\n  Cells Info\n" + "=" * 20)

    # Check data matrix
    # -----------------

    print(f"Matrix sparsity: {adata.X.nnz} filled elements "
          f"({adata.X.nnz / adata.X.shape[0] / adata.X.shape[1]:0.2f}) "
          f"out of {adata.X.shape[0] * adata.X.shape[1]}")

    # In[37]:

    from sklearn.neighbors import NearestNeighbors
    from banksycsr_operations import remove_greater_than, row_normalize
    from banksymain import p_equiv_radius, generate_spatial_weights_fixed_radius, generate_spatial_weights_fixed_nbrs
    from utilsplotting import plot_edge_histogram

    #
    # ---------------------------------------------------------------
    # set params
    # ---------------------------------------------------------------

    visualize_weights = True
    num_neighbours = 10  # only for fixed type
    p_outside = 0.05
    sigmas = (-1,)

    # ---------------------------------------------------------------
    #

    # Find median distance to closest neighbours
    # ------------------------------------------

    nbrs = NearestNeighbors(algorithm='ball_tree').fit(adata.obsm["coord_xy"])
    distances, indices = nbrs.kneighbors(n_neighbors=1)
    median_cell_distance = np.median(distances)
    print(f"\nMedian distance to closest cell = {median_cell_distance}\n")

    # Generate the spatial weights sparse adjacency matrix
    # ----------------------------------------------------

    processing_dict = {}

    for sigma in sigmas:

        if sigma < 0:

            if sigma == -1:
                decay_type = "reciprocal"
            elif sigma == -2:
                decay_type = "ranked"
            elif sigma == -3:
                decay_type = "uniform"
            else:
                raise ValueError("Type not recognised")

            weights_graph, distance_graph = generate_spatial_weights_fixed_nbrs(
                adata.obsm["coord_xy"],
                num_neighbours=num_neighbours,
                decay_type=decay_type,
                nbr_object=nbrs,
                verbose=False,
            )

        else:

            weights_graph, distance_graph = generate_spatial_weights_fixed_radius(
                adata.obsm["coord_xy"],
                p=p_outside,
                nbr_object=nbrs,
                sigma=sigma * median_cell_distance,
                decay_type="gaussian",
                max_num_neighbours=None,
                verbose=False,
            )

        processing_dict[sigma] = {"weights": weights_graph}

        # plot histograms of distances/weights between cells
        fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
        plot_edge_histogram(distance_graph, ax[0], title="distances between cells")
        plot_edge_histogram(weights_graph, ax[1], title="weights between cells")
        #     fig.suptitle(f"weights and distances for sigma = {sigma}")
        fig.tight_layout()

    print(processing_dict)

    # Visualize weights
    # -----------------

    # print(f"weights_graph for sigma {sigma}:\n{weights_graph}\n")

    if visualize_weights:

        from utilsplotting import plot_graph_weights

        for sigma in sigmas:
            plot_graph_weights(
                adata.obsm["coord_xy"],
                processing_dict[sigma]["weights"],
                # max_weight=.5,
                max_weight=1,
                # markersize=1.5,
                markersize=1,
                figsize=(8, 8),
                title=f"Sigma = {sigma}"
            )
            ax = plt.gca()
            ax.axis("equal")

            # Crop spatial graph for readme img
            # ---------------------------------
            # ax.set_xlim(1000,4000)
            # ax.set_ylim(1000,4000)

    # In[38]:

    from banksymain import zscore, weighted_concatenate, banksy_matrix_to_adata

    from datetime import datetime
    time_str = datetime.now().strftime("%b_%d_%Y_%H_%M")
    print(time_str)

    # List Genes
    # ----------

    gene_list = adata.var.index
    # print(f"\n{len(gene_list)} genes to be analysed:\n{gene_list}")

    #
    # -------------------------------------------------------------------------------------
    #                             Neighbourhood settings (lambda)
    # -------------------------------------------------------------------------------------

    nbrhood_contributions = [0.2, ]  # fraction of total dissimilarity from neighbouring cells

    # --------------------------------------------------------------------------------------
    #

    for sigma in processing_dict:

        print(sigma, processing_dict[sigma])

        # Compute neighbour averaged feature matrix
        # -----------------------------------------

        weights = processing_dict[sigma]["weights"]
        neighbour_agg_matrix = weights @ adata.X

        # Save the concatenated un-weighted data (for visualization later)
        # ----------------------------------------------------------------

        if sparse.issparse(adata.X):
            concatenated = sparse.hstack((adata.X, neighbour_agg_matrix), )
        else:
            concatenated = np.concatenate((adata.X, neighbour_agg_matrix), axis=1, )

        processing_dict[sigma]["norm_counts_concatenated"] = concatenated

        print(
            f"\nCell by gene matrix ({type(adata.X)}) has shape {adata.shape}\n"
            f"Normalized spatial weights graph adjacency matrix has shape {weights.shape}\n"
            f"Aggregated Neighbours matrix has shape {neighbour_agg_matrix.shape}\n"
        )

        # Concatenate self / neighbour data
        # ---------------------------------

        for nbrhood_contribution in nbrhood_contributions:

            banksy_matrix = weighted_concatenate(
                zscore(adata.X, axis=0),
                zscore(neighbour_agg_matrix, axis=0),
                nbrhood_contribution
            )

            # plot standard deviations per gene / nbr gene
            # --------------------------------------------

            if sparse.issparse(banksy_matrix):
                st_dev_pergene = banksy_matrix.todense().std(axis=0)
            else:
                st_dev_pergene = banksy_matrix.std(axis=0)

            fig, ax = plt.subplots(figsize=(8, 2))
            ax.bar(np.arange(len(st_dev_pergene)), st_dev_pergene,
                   width=1, color='slateblue', linewidth=0)
            ax.set_title(f"Standard deviations for neighbourhood contribution = {nbrhood_contribution}")

            # save as a new AnnData object
            # ----------------------------

            new_adata = banksy_matrix_to_adata(banksy_matrix, adata)
            # display(new_adata)

            processing_dict[sigma][nbrhood_contribution] = {
                "adata": new_adata,
            }

            # save the banksy matrix as a csv
            # -------------------------------

            save = False
            subfolder = os.path.join("slide_seq", "v2")
            save_folder = os.path.join("data", subfolder)
            save_name = f"adata_{sigma}_l{nbrhood_contribution}_{time_str}.csv"

            if save:

                try:

                    if not os.path.exists(save_folder):
                        os.makedirs(save_folder)

                    new_adata.write(filename=save_name)

                except PermissionError:
                    print("\nWARNING: Permission denied to save file. Not saving adata.\n")

    # set up an entry in the nested dictionary for non-spatial results
    # ----------------------------------------------------------------

    print_max_min(adata)

    # first we need to z-score the data matrix
    new_adata = adata.copy()
    new_adata.X = zscore(new_adata.X, axis=0)

    # here, sigma = 0 and nbr-contributions are also 0, original adata is assigned here as well
    processing_dict[0] = {0.0: {"adata": new_adata, }}

    # In[39]:

    from sklearn.decomposition import PCA
    import umap

    from utilspca import plot_remaining_variance, plot_singular_values

    #
    # -------------------------------------------------------------------------------------
    #                               PCA settings
    # -------------------------------------------------------------------------------------

    pca_dims = [16, ]  # number of PCA dimensions to reduce to, if autodetect_pcs is False

    # --------------------------------------------------------------------------------------
    #
    #

    print(list(processing_dict.keys()))

    for sigma in processing_dict:

        for nbrhood_contribution in processing_dict[sigma]:

            if isinstance(nbrhood_contribution, str):
                continue  # skip weights matrices

            print(
                f"\nReducing sigma = {sigma}, nbrhood contribution = {nbrhood_contribution}\n"
                + "=" * 50 + "\n"
            )

            # Retrieve anndata object
            # -----------------------

            adata_temp = processing_dict[sigma][nbrhood_contribution]["adata"]

            if sparse.issparse(adata_temp.X):
                X = adata_temp.X.todense()
                print("Converting sparse data matrix to dense")
            else:
                X = adata_temp.X

            print(X.shape)

            # Reduce dimensions by PCA and then UMAP
            # --------------------------------------

            for pca_dim in pca_dims:
                pca = PCA(n_components=pca_dim)
                reduced = pca.fit_transform(X)

                print(
                    f"\nReduced {reduced.shape} matrix:\n" + "-" * 60
                    + f"\nmin_value = {reduced.min()}, max = {reduced.max()}\n\n"
                )

                # Plot variance contribution for each component (elbow plot)
                # ----------------------------------------------------------

                plot_remaining_variance(
                    pca, figsize=(5, 2),
                    title=f"sigma = {sigma}, nbrhood contrib = {nbrhood_contribution}"
                )

                # UMAP
                # ----

                reducer = umap.UMAP(transform_seed=42)
                umap_embedding = reducer.fit_transform(reduced)
                print(f"shape of UMAP embedding: {umap_embedding.shape}\n\n")

                # save in dictionary
                # ------------------

                adata_temp.obsm[f"reduced_pc_{pca_dim}"] = reduced
                adata_temp.obsm[f"reduced_pc_{pca_dim}_umap"] = umap_embedding

    # processing_dict, processing_dict.keys(), processing_dict[0].keys()

    # In[40]:

    from banksymain import LeidenPartition
    import leidenalg
    from banksylabels import Label
    import re

    resolutions = [
        #     0.4,
        0.6
    ]
    num_nn = 30
    num_iterations = 100
    partition_seed = 1234  # match vipul's seed
    match_labels = True  # whether to keep labels consistent across different hyperparameter setting

    results = {}

    print(list(processing_dict.keys()))

    for sigma in processing_dict:

        print("Sigma: ", sigma)

        for nbrhood_contribution in processing_dict[sigma]:

            if isinstance(nbrhood_contribution, str):
                continue  # skip weights matrices

            print("Neighbourhood Contribution: ", nbrhood_contribution)

            adata_temp = processing_dict[sigma][nbrhood_contribution]["adata"]
            print("adata temp", adata_temp)
            pca_dims = []

            for key in adata_temp.obsm_keys():

                print(key, "\n")
                match = re.search(r"reduced_pc_([0-9]+)$", key)
                if match:
                    pca_dims.append(int(match.group(1)))

            print("PCA dims to analyse:", pca_dims)

            for pca_dim in pca_dims:

                if isinstance(pca_dim, str):
                    continue  # skip full concatenated matrices

                print("\n" + "=" * 100 +
                      f"\nSetting up partitioner for (sigma = {sigma}, "
                      f"Neighbourhood contribution = {nbrhood_contribution}, "
                      f"PCA dimensions = {pca_dim})\n" + "=" * 100 + "\n")

                banksy_reduced = adata_temp.obsm[f"reduced_pc_{pca_dim}"]
                #             banksy_reduced = processing_dict[sigma][nbrhood_contribution][pca_dim]["reduced"]

                partitioner = LeidenPartition(
                    banksy_reduced,
                    num_nn=num_nn,
                    nns_have_weights=True,
                    compute_shared_nn=True,
                    filter_shared_nn=True,
                    shared_nn_max_rank=3,
                    shared_nn_min_shared_nbrs=5,
                    verbose=True,
                )

                for resolution in resolutions:
                    print(f"\nResolution: {resolution}\n" + "-" * 30 + "\n")

                    # partition
                    # ---------

                    label, modularity = partitioner.partition(
                        resolution=resolution,
                        partition_metric=leidenalg.RBConfigurationVertexPartition,
                        # partition_metric=leidenalg.CPMVertexPartition,
                        n_iterations=num_iterations,
                        seed=partition_seed,
                    )

                    # store results in dictionary
                    # ---------------------------

                    param_str = f"s{sigma}_pc{pca_dim}_nc{nbrhood_contribution:0.2f}_r{resolution:0.2f}"
                    results[param_str] = {
                        "sigma": sigma,
                        "nbrhood_contribution": nbrhood_contribution,
                        "num_pcs": pca_dim,
                        "resolution": resolution,
                        "num_labels": label.num_labels,
                        "labels": label,
                        "adata": processing_dict[sigma][nbrhood_contribution]["adata"]
                    }

    # convert results dictionary into dataframe
    # -----------------------------------------

    results_df = pd.DataFrame.from_dict(results, orient='index')
    # display(results_df)

    results_df.sort_values(by=["sigma", "nbrhood_contribution", "num_pcs", "resolution"],
                           ascending=True, inplace=True)

    print("\nSorted Dataframe:\n")

    print(f"Shape of dataframe: {results_df.shape}")

    # find the maximum number of labels across all parameter sets
    # -----------------------------------------------------------

    max_num_labels = results_df["num_labels"].max()
    print(f"Maximum number of labels: {max_num_labels}\n")

    # Match the labels
    # ----------------

    if match_labels:
        from banksylabels import match_label_series

        results_df["relabeled"], max_num_labels = match_label_series(
            results_df["labels"],
            extra_labels_assignment="greedy",
            verbose=False,
        )


    # In[41]:

    from utilsplotting import plot_2d_embeddings, plot_labels_seperately, plot_label_subset
    from banksylabels import plot_connections

    import matplotlib as mpl
    import matplotlib.gridspec as gridspec

    use_sc_plot = False  # Plot another location plot using scanpy's embedding function
    save_all_h5ad = True  # save the adata objects

    # Generate discrete colormap
    # --------------------------
    # colours borrowed from spagcn tutorial, not sure where they got them from originally
    # can also use standard matplotlib colormaps e.g. spectral or tab10 or tab20

    plot_color = [
        "#F56867", "#FEB915", "#C798EE", "#59BE86", "#7495D3",
        "#D1D1D1", "#6D1A9C", "#15821E", "#3A84E6", "#997273",
        "#787878", "#DB4C6C", "#9E7A7A", "#554236", "#AF5F3C",
        "#93796C", "#F9BD3F", "#DAB370", "#877F6C", "#268785"
    ]
    cmap_discrete = mpl.colors.ListedColormap(plot_color[:max_num_labels])

    for params_name in results_df.index:

        if match_labels:
            labels = results_df.loc[params_name, "relabeled"]
        else:
            labels = results_df.loc[params_name, "labels"]

        adata_temp = results_df.loc[params_name, "adata"]
        num_pcs = results_df.loc[params_name, "num_pcs"]
        pc_temp = adata_temp.obsm[f"reduced_pc_{num_pcs}"]
        umap_temp = adata_temp.obsm[f"reduced_pc_{num_pcs}_umap"]

        label_name = f"labels_{params_name}"
        # adata_temp.obs[label_name] = labels.dense
        adata_temp.obs[label_name] = np.char.mod('%d', labels.dense)
        adata_temp.obs[label_name] = adata_temp.obs[label_name].astype('category')

        # labels.dense to tsv here

        from numpy import savetxt

        output = labels.dense.astype(int)
        savetxt('clustering.tsv', output, delimiter=',', fmt = '%d')




        # Main Figure
        # -----------

        adata_temp.obsm["coord_xy"] = np.vstack(
            (adata_temp.obs["xcoord"].values, adata_temp.obs["ycoord"].values)
        ).T

        if use_sc_plot:
            sc.pl.embedding(adata_temp, basis="coord_xy", color=label_name, size=5)

        fig = plt.figure(figsize=(15, 9), constrained_layout=True)
        grid = fig.add_gridspec(ncols=4, nrows=3,
                                width_ratios=[2, 0.1, 0.5, 0.5],
                                height_ratios=[1, 0.3, 1])

        # Plot Labels
        # -----------

        # cmap_name = 'tab20c'
        cmap_name = cmap_discrete

        ax_locs = fig.add_subplot(grid[:, 0])

        scatterplot = ax_locs.scatter(adata_temp.obsm["coord_xy"][:, 0],
                                      adata_temp.obsm["coord_xy"][:, 1],
                                      c=labels.dense,
                                      cmap=cmap_name,
                                      vmin=0, vmax=max_num_labels - 1,
                                      s=1.5, alpha=1.0)
        ax_locs.set_aspect('equal', 'datalim')
        ax_locs.set_title(f'BANKSY Labels ({params_name})', fontsize=20, fontweight="bold")
        # ax_locs.set_ylim(ax_locs.get_ylim()[::-1])

        ax_cbar = fig.add_subplot(grid[:, 1])
        cbar = fig.colorbar(
            scatterplot,
            boundaries=np.arange(max_num_labels + 1) - 0.5,
            cax=ax_cbar,
        )
        cbar.set_ticks(labels.ids)
        cbar.set_ticklabels(labels.ids)

        # Seperate Location plots
        # -----------------------

        plot_labels_seperately(
            labels, adata_temp.obsm["coord_xy"],
            embeddings=umap_temp,
            cmap_name=cmap_name,
            # cmap_name = None,
            default_colour="tab:red",
            plots_per_row=3,
            spot_size=0.3,
            subplot_size=(3, 3),
            flip_axes=False,
            verbose=False,
        )

        # Plot UMAP (again but with labels)
        # ---------------------------------

        ax_umap = fig.add_subplot(grid[0, -2:])
        plot_2d_embeddings(umap_temp, labels.dense,
                           method_str="UMAP",
                           space_str="",
                           xlabel="UMAP 1", ylabel="UMAP 2",
                           ax=ax_umap,
                           cmap_name=cmap_name,
                           plot_cmap=False,
                           )

        # Plot 1st 2 dimensions of PCA
        # ----------------------------

        dim_sets = (pc_temp[:, :2], pc_temp[:, 1:3])
        dims1 = (0, 1,)
        dims2 = (1, 2,)
        axes = [fig.add_subplot(grid[1, 2 + axnum]) for axnum in range(2)]

        for dim_set, dim1, dim2, ax in zip(dim_sets, dims1, dims2, axes):
            plot_2d_embeddings(dim_set, labels.dense,
                               method_str=f"PCA {dim1 + 1} / {dim2 + 1}",
                               space_str="",
                               xlabel=f"PCA {dim1 + 1}", ylabel=f"PCA {dim2 + 1}",
                               ax=ax,
                               cmap_name=cmap_name,
                               plot_cmap=False,
                               title_fontsize=9)

        # Plot connectivity between labels
        # --------------------------------

        ax_connections = fig.add_subplot(grid[-1, -2:])
        plot_connections(
            labels, weights_graph,
            ax_connections,
            zero_self_connections=True,
            title_str="connections between label",
            colormap_name=cmap_name,
        )


banksy(datasetID, sampleID, Nclusters)





