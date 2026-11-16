"""
Step 3: Aggregate raw per-run timings into summary statistics.

Reads ../results/raw_timings.csv (one row per method x setting x rep) and
computes, for each (method, n_samples, n_genes, sweep) combination, the
mean, standard deviation, min, and max of both elapsed_sec and
peak_mem_mb across the N_REPS repetitions. Writes the result to
../results/summary_timings.csv.

"""

import csv
from collections import defaultdict
from statistics import mean, pstdev

import config


def load_raw_timings():
    if not config.RAW_TIMINGS_PATH.exists():
        raise FileNotFoundError(
            f"{config.RAW_TIMINGS_PATH} not found. Run 02_run_benchmark.py first."
        )
    with open(config.RAW_TIMINGS_PATH, newline="") as f:
        return list(csv.DictReader(f))


def main():
    rows = load_raw_timings()

    groups = defaultdict(list)
    for row in rows:
        key = (row["method"], int(row["n_samples"]), int(row["n_genes"]), row["sweep"])
        groups[key].append(row)

    summary_rows = []
    for (method, n_samples, n_genes, sweep), group_rows in sorted(groups.items()):
        elapsed = [float(r["elapsed_sec"]) for r in group_rows]
        mem = [float(r["peak_mem_mb"]) for r in group_rows]
        n_edges = int(group_rows[0]["n_edges"])

        summary_rows.append(
            {
                "method": method,
                "sweep": sweep,
                "n_samples": n_samples,
                "n_genes": n_genes,
                "n_edges": n_edges,
                "n_reps": len(group_rows),
                "elapsed_sec_mean": f"{mean(elapsed):.6f}",
                "elapsed_sec_sd": f"{pstdev(elapsed) if len(elapsed) > 1 else 0.0:.6f}",
                "elapsed_sec_min": f"{min(elapsed):.6f}",
                "elapsed_sec_max": f"{max(elapsed):.6f}",
                "peak_mem_mb_mean": f"{mean(mem):.3f}",
                "peak_mem_mb_sd": f"{pstdev(mem) if len(mem) > 1 else 0.0:.3f}",
                "peak_mem_mb_min": f"{min(mem):.3f}",
                "peak_mem_mb_max": f"{max(mem):.3f}",
            }
        )

    fieldnames = [
        "method",
        "sweep",
        "n_samples",
        "n_genes",
        "n_edges",
        "n_reps",
        "elapsed_sec_mean",
        "elapsed_sec_sd",
        "elapsed_sec_min",
        "elapsed_sec_max",
        "peak_mem_mb_mean",
        "peak_mem_mb_sd",
        "peak_mem_mb_min",
        "peak_mem_mb_max",
    ]
    with open(config.SUMMARY_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote summary ({len(summary_rows)} rows) to {config.SUMMARY_PATH}")


if __name__ == "__main__":
    main()
