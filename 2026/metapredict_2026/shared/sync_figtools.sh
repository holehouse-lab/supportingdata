#!/usr/bin/env bash
# Copy the master figtools.py into every figure directory that has notebooks, so
# each directory stays self-contained but all of them share one style/colour set.
# Run from anywhere: bash shared/sync_figtools.sh
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(dirname "$here")"
for dir in "$root"/maintext_figures/figure_* "$root"/supplementary_figures/figure_S*; do
    # Figure 3 (tree of life) is copied verbatim from the original analysis and is left untouched.
    if [[ "$(basename "$dir")" == "figure_3" ]]; then
        continue
    fi
    cp "$here/figtools.py" "$dir/figtools.py"
    echo "synced $dir/figtools.py"
done
