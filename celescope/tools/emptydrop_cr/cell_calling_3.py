import random
from collections import namedtuple

import numpy as np
import numpy.ma as ma

import celescope.tools.emptydrop_cr.sgt as cr_sgt  # # modified sgt.py
import celescope.tools.emptydrop_cr.stats as cr_stats  # # modified stats.py
from celescope.tools.matrix import CountMatrix


# Set random seed
random.seed(0)
np.random.seed(0)

# Number of additional barcodes to consider after the initial cell calling
N_CANDIDATE_BARCODES = 20000

# Number of partitions (max number of barcodes to consider for ambient estimation)
N_PARTITIONS = 90000

# Drop this top fraction of the barcodes when estimating ambient.
MAX_OCCUPIED_PARTITIONS_FRAC = 0.5

# Minimum number of UMIS per barcode to consider after the initial cell calling
MIN_UMIS = 500

# Minimum ratio of UMIs to the median (initial cell call UMI) to consider after the initial cell calling
MIN_UMI_FRAC_OF_MEDIAN = 0.01

# Maximum adjusted p-value to call a barcode as non-ambient
MAX_ADJ_PVALUE = 0.01


def adjust_pvalue_bh(p):
    """ Multiple testing correction of p-values using the Benjamini-Hochberg procedure """
    descending = np.argsort(p)[::-1]
    # q = p * N / k where p = p-value, N = # tests, k = p-value rank
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(scale * p[descending]))

    # Return to original order
    return q[np.argsort(descending)]


def estimate_profile_sgt(matrix, barcode_indices, nz_feat):
    """ Estimate a gene expression profile by Simple Good Turing.
    Args:
      raw_mat (sparse matrix): Sparse matrix of all counts
      barcode_indices (np.array(int)): Barcode indices to use
      nz_feat (np.array(int)): Indices of features that are non-zero at least once
    Returns:
      profile (np.array(float)): Estimated probabilities of length len(nz_feat).
    """
    # Initial profile estimate
    prof_mat = matrix[:, barcode_indices]

    profile = np.ravel(prof_mat[nz_feat, :].sum(axis=1))
    zero_feat = np.flatnonzero(profile == 0)

    # Simple Good Turing estimate
    p_smoothed, p0 = cr_sgt.sgt_proportions(profile[np.flatnonzero(profile)])

    # Distribute p0 equally among the zero elements.
    p0_i = p0/len(zero_feat)

    profile_p = np.repeat(p0_i, len(nz_feat))
    profile_p[np.flatnonzero(profile)] = p_smoothed

    assert np.isclose(profile_p.sum(), 1.0)
    return profile_p


# Construct a background expression profile from barcodes with <= T UMIs
def est_background_profile_sgt(matrix, use_bcs):
    """ Estimate a gene expression profile on a given subset of barcodes.
         Use Good-Turing to smooth the estimated profile.
    Args:
      matrix (scipy.sparse.csc_matrix): Sparse matrix of all counts
      use_bcs (np.array(int)): Indices of barcodes to use (col indices into matrix)
    Returns:
      profile (use_features, np.array(float)): Estimated probabilities of length use_features.
    """
    # Use features that are nonzero anywhere in the data
    use_feats = np.flatnonzero(np.asarray(matrix.sum(1)))

    # Estimate background profile
    bg_profile_p = estimate_profile_sgt(matrix, use_bcs, use_feats)

    return (use_feats, bg_profile_p)


def find_nonambient_barcodes(raw_mat, recovered_cells,
                             min_umi_frac_of_median=MIN_UMI_FRAC_OF_MEDIAN,
                             min_umis_nonambient=MIN_UMIS,
                             max_adj_pvalue=MAX_ADJ_PVALUE,):
    """ Call barcodes as being sufficiently distinct from the ambient profile
    Args:
      raw_mat: raw matrix of UMI counts
      recovered_cells: expected number of recovered cells
    Returns:
    TBD
    """
    NonAmbientBarcodeResult = namedtuple('NonAmbientBarcodeResult',
                                         ['eval_bcs',      # Candidate barcode indices (n)
                                          'log_likelihood',  # Ambient log likelihoods (n)
                                          'pvalues',       # pvalues (n)
                                          'pvalues_adj',   # B-H adjusted pvalues (n)
                                          'is_nonambient',  # Boolean nonambient calls (n)
                                          ])

    # Estimate an ambient RNA profile
    umis_per_bc = np.squeeze(np.asarray(raw_mat.sum(axis=0)))
    # get the index of sorted umis_per_bc (ascending, bc_order[0] is the index of the smallest element in umis_per_bc)
    bc_order = np.argsort(umis_per_bc)

    # Require non-zero barcodes
    nz_bcs = np.flatnonzero(umis_per_bc)

    # Take what we expect to be the barcodes associated w/ empty partitions.
    # Drop the top MAX_OCCUPIED_PARTITIONS_FRAC fraction of barcodes
    empty_bcs = bc_order[::-1][int(N_PARTITIONS*MAX_OCCUPIED_PARTITIONS_FRAC):N_PARTITIONS]
    empty_bcs.sort()

    # Get the non-zero barcodes overlapped with barcodes assumed to be associated w/ empty partitions
    use_bcs = np.intersect1d(empty_bcs, nz_bcs, assume_unique=True)

    if len(use_bcs) > 0:
        try:
            # Get used "Gene" features (eval_features)
            # and the smoothed prob profile per "Gene" (ambient_profile_p)
            eval_features, ambient_profile_p = est_background_profile_sgt(raw_mat.tocsc(), use_bcs)
        except cr_sgt.SimpleGoodTuringError as e:
            print(str(e))
    else:
        eval_features = np.zeros(0, dtype=int)
        ambient_profile_p = np.zeros(0)

    # Choose candidate cell barcodes
    # Regular ordmag filter
    gg_filtered_indices, gg_filtered_metrics, _msg = cr_stats.filter_cellular_barcodes_ordmag(
        umis_per_bc, recovered_cells=recovered_cells)

    print('Cell-called barcodes metrics:')
    print('\n'.join(list(map(lambda x: '{}: {}'.format(*x), list(gg_filtered_metrics.items())))))
    print('==============================')

    orig_cell_bc_set = set(gg_filtered_indices)
    orig_cells = np.flatnonzero(np.fromiter((bc in orig_cell_bc_set for bc in range(raw_mat.shape[1])), dtype=bool))

    # No good incoming cell calls
    if orig_cells.sum() == 0:
        print('Error: No original cells are selected!')
        return None, None, None

    # Look at non-cell barcodes above a minimum UMI count
    eval_bcs = np.ma.array(np.arange(raw_mat.shape[1]))
    eval_bcs[orig_cells] = ma.masked

    median_initial_umis = np.median(umis_per_bc[orig_cells])

    min_umis = int(max(min_umis_nonambient, round(np.ceil(median_initial_umis * min_umi_frac_of_median))))

    print('Median UMIs of initial cell calls: {}'.format(median_initial_umis))
    print('Min UMIs: {}'.format(min_umis))

    eval_bcs[umis_per_bc < min_umis] = ma.masked
    n_unmasked_bcs = len(eval_bcs) - eval_bcs.mask.sum()

    # Take the top N_CANDIDATE_BARCODES by UMI count, of barcodes that pass the above criteria
    # For evaluation of non-ambient bcs using background info estimated from SGT
    eval_bcs = np.argsort(ma.masked_array(umis_per_bc, mask=eval_bcs.mask))[:n_unmasked_bcs][-N_CANDIDATE_BARCODES:]

    if len(eval_bcs) == 0:
        print('Warning: no eval bcs are selected to evaluate non-empty bcs from SGT results!')
        print('Output bcs from 1st round cell calling ONLY.')
        return orig_cells, gg_filtered_metrics, None
    else:
        assert not np.any(np.isin(eval_bcs, orig_cells))
        print('Number of candidate bcs: {}'.format(len(eval_bcs)))
        print('Range candidate bc umis: {}, {}'.format(umis_per_bc[eval_bcs].min(), umis_per_bc[eval_bcs].max()))

        eval_mat = raw_mat.tocsc()[eval_features, :][:, eval_bcs]

        if len(ambient_profile_p) == 0:
            obs_loglk = np.repeat(np.nan, len(eval_bcs))
            pvalues = np.repeat(1, len(eval_bcs))
            sim_loglk = np.repeat(np.nan, len(eval_bcs))

        # Compute observed log-likelihood of barcodes being generated from ambient RNA
        obs_loglk = cr_stats.eval_multinomial_loglikelihoods(eval_mat, ambient_profile_p)

        # Simulate log likelihoods
        distinct_ns, sim_loglk = cr_stats.simulate_multinomial_loglikelihoods(ambient_profile_p, umis_per_bc[eval_bcs],
                                                                              num_sims=10000, verbose=True)

        # Compute p-values
        pvalues = cr_stats.compute_ambient_pvalues(umis_per_bc[eval_bcs], obs_loglk, distinct_ns, sim_loglk)

        pvalues_adj = adjust_pvalue_bh(pvalues)
        max_adj_pvalue = 0.01
        is_nonambient = pvalues_adj <= max_adj_pvalue

        print('Number of non-ambient barcodes from SGT:', len(eval_bcs[is_nonambient]))

        # Runxi's filtering
        print('Identify {} cell-associated barcodes'.format(len(orig_cells)+len(eval_bcs[is_nonambient])))

        # of barcodes overlapped w/ the cellranger results
        filtered_bc_indices = np.concatenate((orig_cells, eval_bcs[is_nonambient]), axis=None)

        return filtered_bc_indices, gg_filtered_metrics, NonAmbientBarcodeResult(
            eval_bcs=eval_bcs,
            log_likelihood=obs_loglk,
            pvalues=pvalues,
            pvalues_adj=pvalues_adj,
            is_nonambient=is_nonambient,
        )


def cell_calling_3(all_matrix_10X_dir, expected_cell_num):

    count_matrix = CountMatrix.from_matrix_dir(matrix_dir=all_matrix_10X_dir)

    # Run cell calling
    filtered_bc_indices, round_1_filtered_metrics, _non_ambient_barcode_result = find_nonambient_barcodes(
        raw_mat=count_matrix.get_matrix(), recovered_cells=expected_cell_num)

    raw_barcodes = np.array(count_matrix.get_barcodes())
    cell_bc = raw_barcodes[filtered_bc_indices]
    initial_cell_num = round_1_filtered_metrics['filtered_bcs']
    return cell_bc, initial_cell_num
