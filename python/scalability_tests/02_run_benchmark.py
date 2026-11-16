"""
Step 2: Run the runtime/memory benchmark.

For every synthetic expression matrix listed in ../data/manifest.csv, and
for every method in config.METHODS (LIONESS, SSN, SWEET, BONOBO), this
script:
  1. loads the matrix,
  2. calls ssn_methods.run_SSN_method(expr, method_ID) once,
  3. records wall-clock runtime (seconds) and peak memory (MB) for that
     single call,
  4. appends one row of results to ../results/raw_timings.csv.

Each (n_samples, n_genes, rep) combination in the manifest already
corresponds to one of the N_REPS repetitions configured in config.py, so
running this script once covers all repetitions.

Peak memory is measured with resource.getrusage(), which reports the
process's maximum resident set size since the Python interpreter started.
Because memory is monotonically non-decreasing across calls within one
process, this script launches a FRESH PYTHON SUBPROCESS for every single
(method, file) run, so that the reported peak memory reflects only that
one run and is not inflated by previous runs accumulated in the same
process.


This can take a while depending on the grid size; progress is printed to
stdout as it goes, and results are appended incrementally to
raw_timings.csv so a partially completed run is not lost.
"""

import csv
import json
import subprocess
import sys

import config

# Template for the subprocess that times and measures memory for exactly
# one (method, data file) run. Runs in its own process so RSS reflects
# only this run.
_WORKER_TEMPLATE = """
import resource
import sys
import time
import numpy as np

sys.path.insert(0, {code_dir!r})
from ssn_methods import run_SSN_method

expr = np.load({data_path!r})
method_id = {method_id!r}

start = time.perf_counter()
result = run_SSN_method(expr, method_id)
elapsed = time.perf_counter() - start

usage = resource.getrusage(resource.RUSAGE_SELF)
# ru_maxrss is in KB on Linux, bytes on macOS.
peak_rss_kb = usage.ru_maxrss if sys.platform != "darwin" else usage.ru_maxrss / 1024.0
peak_mem_mb = peak_rss_kb / 1024.0

print("RESULT_JSON:" + __import__("json").dumps(
    {{"elapsed_sec": elapsed, "peak_mem_mb": peak_mem_mb}}
))
"""


def run_one(data_path, method_id):
    """Run a single (method, data file) benchmark in an isolated subprocess."""
    script = _WORKER_TEMPLATE.format(
        code_dir=str(config.CODE_DIR),
        data_path=str(data_path),
        method_id=method_id,
    )
    proc = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Benchmark subprocess failed for {data_path} (method_id={method_id}):\n"
            f"{proc.stderr}"
        )

    for line in proc.stdout.splitlines():
        if line.startswith("RESULT_JSON:"):
            payload = json.loads(line[len("RESULT_JSON:"):])
            return payload["elapsed_sec"], payload["peak_mem_mb"]

    raise RuntimeError(
        f"No RESULT_JSON found in subprocess output for {data_path} "
        f"(method_id={method_id}). Stdout was:\n{proc.stdout}"
    )


def load_manifest():
    if not config.MANIFEST_PATH.exists():
        raise FileNotFoundError(
            f"Manifest not found at {config.MANIFEST_PATH}. "
            "Run 01_generate_data.py first."
        )
    with open(config.MANIFEST_PATH, newline="") as f:
        return list(csv.DictReader(f))


def main():
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest()

    fieldnames = [
        "method",
        "n_samples",
        "n_genes",
        "n_edges",
        "rep",
        "sweep",
        "elapsed_sec",
        "peak_mem_mb",
        "file",
    ]

    write_header = not config.RAW_TIMINGS_PATH.exists()
    with open(config.RAW_TIMINGS_PATH, "a", newline="") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()

        total_runs = len(manifest) * len(config.METHODS)
        run_idx = 0

        for row in manifest:
            data_path = config.DATA_DIR / row["file"]

            for method_name, method_id in config.METHODS.items():
                run_idx += 1
                print(
                    f"[{run_idx}/{total_runs}] {method_name} on {row['file']} "
                    f"(n_samples={row['n_samples']}, n_genes={row['n_genes']}, "
                    f"rep={row['rep']}) ...",
                    end=" ",
                    flush=True,
                )

                elapsed_sec, peak_mem_mb = run_one(data_path, method_id)

                writer.writerow(
                    {
                        "method": method_name,
                        "n_samples": row["n_samples"],
                        "n_genes": row["n_genes"],
                        "n_edges": row["n_edges"],
                        "rep": row["rep"],
                        "sweep": row["sweep"],
                        "elapsed_sec": f"{elapsed_sec:.6f}",
                        "peak_mem_mb": f"{peak_mem_mb:.3f}",
                        "file": row["file"],
                    }
                )
                out_f.flush()

                print(f"{elapsed_sec:.3f}s, {peak_mem_mb:.1f} MB")

    print(f"\nWrote raw timings to {config.RAW_TIMINGS_PATH}")


if __name__ == "__main__":
    main()
