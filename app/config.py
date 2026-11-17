"""Runtime configuration for the paired FASTQ service.

All limits can be overridden through environment variables, but the defaults
match the delivery requirement: each uploaded file <= 256 MiB and the whole
process is expected to live inside a 256 MiB memory cap (see
docker-compose.yml).
"""

from __future__ import annotations

import os
from pathlib import Path

MAX_FILE_BYTES = int(os.environ.get("FASTQ_MAX_FILE_BYTES", 256 * 1024 * 1024))
"""Maximum size of a single uploaded FASTQ file (256 MiB)."""

DATA_DIR = Path(os.environ.get("FASTQ_DATA_DIR", "/data"))
"""Root directory for the SQLite database and all task files."""

DB_PATH = DATA_DIR / "tasks.db"
TASKS_DIR = DATA_DIR / "tasks"

# Uploads are spooled to a temporary file once they exceed this size, so a big
# multipart upload never accumulates in RAM.
UPLOAD_SPOOL_THRESHOLD = int(
    os.environ.get("FASTQ_UPLOAD_SPOOL_THRESHOLD", 1 * 1024 * 1024)
)

# Records are inserted in batches; only one batch is held in memory at a time.
INSERT_BATCH_SIZE = int(os.environ.get("FASTQ_INSERT_BATCH_SIZE", 512))

# SQLite page cache kept deliberately small. The heavy data stays on disk.
SQLITE_CACHE_PAGES = int(os.environ.get("FASTQ_SQLITE_CACHE_PAGES", 2048))
SQLITE_BUSY_TIMEOUT_MS = 30_000

# How often the worker refreshes the progress percentage in the database.
PROGRESS_UPDATE_EVERY_BYTES = 8 * 1024 * 1024

STAGE_UPLOAD = "upload"
STAGE_PARSE = "parse"
STAGE_VALIDATE = "validate"
STAGE_PUBLISH = "publish"
STAGE_DONE = "done"

# ---- Pair coverage check (see app/coverage.py) ----
COVERAGE_MAX_REFERENCE_BASES = int(
    os.environ.get("FASTQ_COVERAGE_MAX_REFERENCE_BASES", 10_000)
)
"""Maximum length of the single ACGT reference sequence (10,000 bp)."""

COVERAGE_MAX_PAIRS = int(os.environ.get("FASTQ_COVERAGE_MAX_PAIRS", 2_000))
"""Maximum number of read pairs one coverage check may process."""

COVERAGE_MIN_READ_BASES = int(os.environ.get("FASTQ_COVERAGE_MIN_READ_BASES", 40))
COVERAGE_MAX_READ_BASES = int(os.environ.get("FASTQ_COVERAGE_MAX_READ_BASES", 100))
"""Accepted per-end read length range (40-100 bp)."""

COVERAGE_MAX_SPAN = int(os.environ.get("FASTQ_COVERAGE_MAX_SPAN", 500))
"""Maximum outermost covered span of a compliant pair placement (500 bp)."""

COVERAGE_MIN_DEPTH_THRESHOLD = 1
COVERAGE_MAX_DEPTH_THRESHOLD = 20
"""Accepted depth-threshold range (1-20)."""

COVERAGE_CACHE_POSITIONS = int(
    os.environ.get("FASTQ_COVERAGE_CACHE_POSITIONS", 1_000_000)
)
"""Budget of cached match positions (int32, 1M entries = 4 MiB)."""
