#!/usr/bin/env bash
# Run every metapredict benchmark on an Apple-silicon Mac (CPU + MPS), in sequence.
# Run on an otherwise idle machine: timings and peak-RAM numbers are sensitive to
# other load. Results land in ../data/benchmarks/<benchmark>/ as TSVs.
#
#   bash run_local_mac.sh                      # standard set (~1-2 h)
#   bash run_local_mac.sh --with-full-stream   # ... plus the full stream sweep (~4 h more)
#
# Optional environment variables:
#   MACHINE_LABEL           label for this machine in the figures (default: CPU name)
#   BATCH_SWEEP_COOLDOWN    seconds to rest before each batch-size measurement
#                           (default 0; 600 avoids thermal throttling but makes the
#                           sweep ~6 h)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"

label="${MACHINE_LABEL:-}"
cooldown="${BATCH_SWEEP_COOLDOWN:-0}"
full_stream=false
if [[ "${1:-}" == "--with-full-stream" ]]; then
    full_stream=true
fi

python -c "import metapredict; print('metapredict', metapredict.__version__, 'from', metapredict.__file__)"

echo "== 1/4 proteome benchmark (Figure 1E, S5): V1-V3 x CPU, MPS"
for version in 1 2 3; do
    for device in cpu mps; do
        python proteome/benchmark_metapredict.py --version "$version" --device "$device" --machine-label "$label"
    done
done

echo "== 2/4 batch-size sweep (Figure S6): V1-V3 x CPU, MPS"
python batch_size/run_batch_size_sweep.py --devices cpu mps --cooldown "$cooldown" --machine-label "$label"

echo "== 3/4 batch vs sequential (Figure S7): V1-V3"
python batch_vs_sequential/run_batch_vs_sequential.py --machine-label "$label"

echo "== 4/4 streaming memory, quick sweep (Figure S8)"
python stream/benchmark_stream_memory.py --quick --machine-label "$label"

if [[ "$full_stream" == true ]]; then
    # OPTIONAL, slow: 7 datasets up to ~1e9 residues (~4 h on an M3 Max, ~2 GB of
    # synthetic FASTA in ../data/benchmarks/stream_synthetic_data/). Resumable: re-run to continue.
    echo "== optional: streaming memory, FULL sweep (Figure S8)"
    python stream/benchmark_stream_memory.py --machine-label "$label"
fi

echo "done; results are in ../data/benchmarks/*/"
