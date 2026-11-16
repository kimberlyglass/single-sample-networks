"""
Step 1: Generate synthetic expression matrices for the scalability benchmark.

For every (n_samples, n_genes) setting in the grid defined in config.py,
and for each of N_REPS repetitions, this script generates a synthetic gene
expression matrix with a block-correlated structure (modules of
co-expressed genes) and saves it as a .npy file (genes x samples, matching
the LocExp convention expected by ssn_methods.py) under ../data/.

"""

import csv

import numpy as np

import config


def make_block_correlated_expression(n_genes, n_samples, n_blocks, rho, rng):
    """
    Generate a synthetic gene expression matrix (n_genes x n_samples) with
    a block/modular correlation structure: genes are split into n_blocks
    contiguous modules, and genes within the same module share a latent
    factor so that their pairwise correlation is approximately `rho`.
    Genes in different modules are approximately independent.

    """
    block_ids = np.array_split(np.arange(n_genes), n_blocks)

    expr = np.empty((n_genes, n_samples), dtype=np.float64)
    for genes_in_block in block_ids:
        if len(genes_in_block) == 0:
            continue
        block_factor = rng.standard_normal(n_samples)
        noise = rng.standard_normal((len(genes_in_block), n_samples))
        expr[genes_in_block, :] = (
            np.sqrt(rho) * block_factor + np.sqrt(1 - rho) * noise
        )

    return expr


def build_grid():
    """Build the list of (n_samples, n_genes, sweep_name) settings to generate."""
    settings = []

    # Sweep over n_samples at fixed n_genes
    for n_samples in config.N_SAMPLES_GRID:
        settings.append((n_samples, config.FIXED_P_FOR_N_SWEEP, "vary_n_samples"))

    # Sweep over n_genes at fixed n_samples
    for n_genes in config.N_GENES_GRID:
        settings.append((config.FIXED_N_FOR_P_SWEEP, n_genes, "vary_n_genes"))

    # De-duplicate identical (n_samples, n_genes) pairs that may appear in
    # both sweeps (e.g. the shared default point), keeping first occurrence.
    seen = set()
    unique_settings = []
    for n_samples, n_genes, sweep in settings:
        key = (n_samples, n_genes)
        if key in seen:
            continue
        seen.add(key)
        unique_settings.append((n_samples, n_genes, sweep))

    return unique_settings


def main():
    config.DATA_DIR.mkdir(parents=True, exist_ok=True)

    settings = build_grid()

    manifest_rows = []
    for n_samples, n_genes, sweep in settings:
        n_edges = n_genes * (n_genes - 1) // 2
        for rep in range(1, config.N_REPS + 1):
            seed = config.BASE_SEED + 1000 * rep + n_samples + n_genes
            rng = np.random.default_rng(seed)

            expr = make_block_correlated_expression(
                n_genes=n_genes,
                n_samples=n_samples,
                n_blocks=config.N_BLOCKS,
                rho=config.WITHIN_BLOCK_RHO,
                rng=rng,
            )

            fname = f"expr_n{n_samples}_p{n_genes}_rep{rep}.npy"
            fpath = config.DATA_DIR / fname
            np.save(fpath, expr)

            manifest_rows.append(
                {
                    "file": fname,
                    "n_samples": n_samples,
                    "n_genes": n_genes,
                    "n_edges": n_edges,
                    "rep": rep,
                    "sweep": sweep,
                    "seed": seed,
                }
            )

            print(
                f"[generated] {fname}  "
                f"(n_samples={n_samples}, n_genes={n_genes}, sweep={sweep}, rep={rep})"
            )

    with open(config.MANIFEST_PATH, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "file",
                "n_samples",
                "n_genes",
                "n_edges",
                "rep",
                "sweep",
                "seed",
            ],
        )
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"\nWrote manifest with {len(manifest_rows)} entries to {config.MANIFEST_PATH}")


if __name__ == "__main__":
    main()
