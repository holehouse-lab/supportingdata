"""Shared provenance and TSV helpers for the metapredict benchmark scripts.

Every benchmark row records exactly which metapredict was run (version string,
install path and git commit) and on what hardware, so results collected on
different machines can be pooled and checked before plotting. Each benchmark
script imports this module from the parent ``benchmarks/`` directory.
"""

from __future__ import annotations

import csv
import os
import platform
import socket
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

import psutil

# Benchmark inputs and results live in the repository's data directory, so the figure
# notebooks read them from the same place as every other input: each benchmark writes
# its TSVs to data/benchmarks/<benchmark>/ (see benchmark_results_dir).
BENCHMARK_DIR = Path(__file__).resolve().parent
BENCHMARK_DATA_DIR = BENCHMARK_DIR.parent / "data" / "benchmarks"

# The human reference proteome shared by the proteome, batch-size and
# batch-vs-sequential benchmarks (UniProt UP000005640, one sequence per gene).
HUMAN_PROTEOME_FASTA = BENCHMARK_DATA_DIR / "UP000005640_9606.fasta"

BYTES_PER_GIB = 1024.0**3

# Lower-cased fragments of the RuntimeError messages PyTorch raises when an
# allocation fails: "MPS backend out of memory" (MPS) and "DefaultCPUAllocator:
# can't allocate memory" (CPU). CUDA raises the dedicated torch.OutOfMemoryError.
OUT_OF_MEMORY_MESSAGES = ("out of memory", "can't allocate memory")


def benchmark_results_dir(benchmark: str) -> Path:
    """Default directory a benchmark writes its result TSVs to.

    Parameters
    ----------
    benchmark : str
        Benchmark name, which is also its script's subdirectory of ``benchmarks/``:
        ``"proteome"``, ``"batch_size"``, ``"batch_vs_sequential"`` or ``"stream"``.

    Returns
    -------
    pathlib.Path
        ``data/benchmarks/<benchmark>``, where the figure notebooks look for the TSVs
        (``figtools.benchmark_results_path``). The directory is not created here;
        each script creates it before writing.
    """
    return BENCHMARK_DATA_DIR / benchmark


def utc_now() -> datetime:
    """Return the current time in UTC."""
    return datetime.now(timezone.utc)


def utc_stamp(now: Optional[datetime] = None) -> str:
    """Return a compact UTC timestamp for filenames, e.g. ``20261029T171500Z``."""
    return (now or utc_now()).strftime("%Y%m%dT%H%M%SZ")


def short_hostname() -> str:
    """Return the hostname without any domain suffix."""
    return socket.gethostname().split(".")[0]


def _git(repo_dir: str, *args: str) -> str:
    """Run a git command in ``repo_dir``; return stripped stdout or '' if git fails."""
    try:
        completed = subprocess.run(
            ["git", "-C", repo_dir, *args], capture_output=True, text=True, timeout=10
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return completed.stdout.strip() if completed.returncode == 0 else ""


def metapredict_provenance() -> dict[str, str]:
    """Describe the metapredict install being benchmarked.

    Returns
    -------
    dict of str to str
        ``metapredict_version`` (``metapredict.__version__``), ``metapredict_path``
        (the package directory actually imported) and ``metapredict_git_commit``
        (short hash of that directory's git checkout, with ``+dirty`` if tracked
        files are modified; empty when metapredict is not a git checkout, e.g. a
        pip install).
    """
    import metapredict

    package_dir = os.path.dirname(os.path.abspath(metapredict.__file__))
    commit = _git(package_dir, "rev-parse", "--short", "HEAD")
    if commit and _git(package_dir, "status", "--porcelain", "--untracked-files=no"):
        commit += "+dirty"
    return {
        "metapredict_version": str(getattr(metapredict, "__version__", "")),
        "metapredict_path": package_dir,
        "metapredict_git_commit": commit,
    }


def host_record(device: str, machine_label: str = "") -> dict[str, Any]:
    """Describe the host hardware and software stack.

    Parameters
    ----------
    device : str
        torch device used for the measurement (``cpu``, ``mps`` or ``cuda``);
        decides which GPU description is recorded.
    machine_label : str, optional
        Free-text label for this machine (e.g. ``"x86-64 workstation"``). Figures
        use it to name the hardware group; when empty they build a label from the
        CPU / GPU names instead.

    Returns
    -------
    dict
        Flat record of hostname, platform, CPU, RAM, GPU and library versions.
        Fields that cannot be determined are empty strings.
    """
    import torch

    record: dict[str, Any] = {
        "hostname": short_hostname(),
        "machine_label": machine_label,
        "platform": platform.platform(),
        "cpu_brand": "",
        "physical_cores": psutil.cpu_count(logical=False),
        "logical_cores": psutil.cpu_count(logical=True),
        "total_ram_gb": round(psutil.virtual_memory().total / BYTES_PER_GIB, 2),
        "gpu_name": "",
        "gpu_memory_gb": "",
        "python_version": platform.python_version(),
        "torch_version": torch.__version__,
    }
    # py-cpuinfo gives a stable CPU brand string across platforms (optional).
    try:
        import cpuinfo

        record["cpu_brand"] = cpuinfo.get_cpu_info().get("brand_raw", "")
    except ImportError:
        record["cpu_brand"] = platform.processor()

    if device == "cuda" and torch.cuda.is_available():
        idx = torch.cuda.current_device()
        record["gpu_name"] = torch.cuda.get_device_name(idx)
        record["gpu_memory_gb"] = round(
            torch.cuda.get_device_properties(idx).total_memory / BYTES_PER_GIB, 2
        )
    elif device == "mps" and torch.backends.mps.is_available():
        record["gpu_name"] = "Apple MPS"
        record["gpu_memory_gb"] = round(torch.mps.recommended_max_memory() / BYTES_PER_GIB, 2)
    record.update(metapredict_provenance())
    return record


def device_available(device: str) -> bool:
    """True if the torch device (``cpu``, ``mps`` or ``cuda``) can be used here."""
    import torch

    if device == "cpu":
        return True
    if device == "mps":
        return torch.backends.mps.is_available()
    if device == "cuda":
        return torch.cuda.is_available()
    raise ValueError(f"Unknown device {device!r}; expected cpu, mps or cuda")


def sync_device(device: str) -> None:
    """Block until queued work on ``device`` has finished, so timings are complete."""
    import torch

    if device == "cuda":
        torch.cuda.synchronize()
    elif device == "mps":
        torch.mps.synchronize()


def release_device_memory(device: str) -> None:
    """Hand memory PyTorch has cached on ``device`` back to the driver.

    PyTorch's CUDA and MPS allocators keep freed blocks cached for reuse, so over a
    batch-size sweep the cache fills with blocks sized for earlier batch sizes. The
    next run's first large allocation can then fail, and the allocator flushes its
    cache and retries inside the timed region (on an RTX A4000 this printed a
    CUDACachingAllocator OOM warning during the V3 batch-size-256 run even though
    the run succeeded). Calling this off the clock before each timed run starts
    every batch size from the same, empty cache.

    Parameters
    ----------
    device : str
        ``cpu``, ``mps`` or ``cuda``. On ``cpu`` only Python garbage collection runs.
    """
    import gc

    import torch

    gc.collect()
    if device == "cuda":
        torch.cuda.empty_cache()
    elif device == "mps":
        torch.mps.empty_cache()


def is_out_of_memory(error: BaseException) -> bool:
    """True if ``error`` means the device or host ran out of memory.

    Recognises CUDA (``torch.OutOfMemoryError``), MPS and PyTorch's CPU allocator
    (RuntimeErrors whose message matches ``OUT_OF_MEMORY_MESSAGES``) and Python's
    own ``MemoryError``. A process killed by the Linux OOM killer receives SIGKILL
    and never gets the chance to raise, so that case cannot be caught here.

    Parameters
    ----------
    error : BaseException
        The exception raised by a prediction call.

    Returns
    -------
    bool
        True for an out-of-memory failure, False for anything else.
    """
    import torch

    if isinstance(error, (torch.OutOfMemoryError, MemoryError)):
        return True
    if isinstance(error, RuntimeError):
        message = str(error).lower()
        return any(fragment in message for fragment in OUT_OF_MEMORY_MESSAGES)
    return False


def default_device(version: str) -> str:
    """The device metapredict 3.1 auto-selects for ``version`` when none is given.

    Parameters
    ----------
    version : str
        Network version, ``"V1"``, ``"V2"`` or ``"V3"``.

    Returns
    -------
    str
        ``cpu``, ``mps`` or ``cuda``.
    """
    from metapredict.backend.network_parameters import metapredict_networks
    from metapredict.backend.predictor import check_device

    params = metapredict_networks[version]["parameters"]
    return check_device(None, default_device=params.get("device_order", ["cuda", "mps", "cpu"]))


def default_batch_size(version: str, device: str) -> int:
    """The batch size metapredict uses for ``version`` on ``device`` by default.

    Uses metapredict's own resolver (3.1+), which applies the per-network,
    per-device defaults (e.g. 256 on CPU/CUDA, 512 on MPS).

    Parameters
    ----------
    version : str
        Network version, ``"V1"``, ``"V2"`` or ``"V3"``.
    device : str
        ``cpu``, ``mps`` or ``cuda``.

    Returns
    -------
    int
        The batch size a call without ``batch_size=`` resolves to.
    """
    from metapredict.backend.network_parameters import metapredict_networks
    from metapredict.backend.predictor import resolve_batch_size

    return resolve_batch_size(None, device, metapredict_networks[version]["parameters"])


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> Path:
    """Write rows (all with identical keys) to a tab-separated file with a header.

    Parameters
    ----------
    path : pathlib.Path
        Output file; parent directories are created.
    rows : list of dict
        Records to write, in order. Every row must have the same keys.

    Returns
    -------
    pathlib.Path
        The written path.
    """
    assert rows, "refusing to write an empty TSV"
    fields = list(rows[0].keys())
    for row in rows:
        assert list(row.keys()) == fields, f"row keys differ from header: {row.keys()}"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")
    return path


def read_fasta(path: Path) -> dict[str, str]:
    """Read a FASTA file leniently (duplicate headers kept, odd residues converted)."""
    import protfasta

    return protfasta.read_fasta(
        str(path), invalid_sequence_action="convert", expect_unique_header=False
    )
