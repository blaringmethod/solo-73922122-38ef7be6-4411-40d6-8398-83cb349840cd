"""Self-tests for the paired coverage check API.

Run with the project virtualenv:
    .venv/bin/python -I tests/test_coverage.py
Each test run gets an isolated FASTQ_DATA_DIR under tests/output.
"""

from __future__ import annotations

import asyncio
import hashlib
import io
import os
import random
import shutil
import sys
import tempfile
import time
import uuid
from pathlib import Path

# Make the project root importable when executed as a script.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import httpx  # noqa: E402

from app import config, db, worker  # noqa: E402
from app.main import app  # noqa: E402

RESULTS: list[tuple[str, bool, str]] = []


def record(name: str, ok: bool, detail: str = "") -> None:
    RESULTS.append((name, ok, detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f" - {detail}" if detail else ""))


def make_record(rid: str, seq: bytes, quality: bytes | None = None,
                suffix: str = "") -> bytes:
    if quality is None:
        quality = b"~" * len(seq)
    return b"\n".join((b"@" + rid.encode() + suffix.encode(), seq, b"+", quality)) + b"\n"


def fastq(records: list[bytes]) -> io.BytesIO:
    return io.BytesIO(b"".join(records))


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


# Deterministic pseudo-random reference: 40-100 bp slices are unique in it.
_RNG = random.Random(20261117)
REF_RANDOM = "".join(_RNG.choice("ACGT") for _ in range(600))


def fasta_text(seq: str, name: str = "ref1") -> str:
    return f">{name} panel\n" + seq + "\n"


async def wait_terminal(client: httpx.AsyncClient, task_id: str,
                        timeout: float = 60.0) -> dict:
    deadline = time.time() + timeout
    while time.time() < deadline:
        body = (await client.get(f"/tasks/{task_id}")).json()
        if body["status"] in ("succeeded", "failed", "cancelled"):
            return body
        await asyncio.sleep(0.05)
    raise AssertionError(f"task {task_id} did not finish in {timeout}s")


async def upload_pairs(client: httpx.AsyncClient,
                       pairs: list[tuple[bytes, bytes]]) -> dict:
    r1 = fastq([make_record(f"pair-{i}", s1) for i, (s1, _) in enumerate(pairs)])
    r2 = fastq([make_record(f"pair-{i}", s2) for i, (_, s2) in enumerate(pairs)])
    resp = await client.post(
        "/tasks",
        files={"r1": ("r1.fastq", r1, "application/octet-stream"),
               "r2": ("r2.fastq", r2, "application/octet-stream")},
    )
    assert resp.status_code == 202, resp.text
    return await wait_terminal(client, resp.json()["id"])


async def coverage(client: httpx.AsyncClient, task_id: str, fasta: str,
                   threshold: int) -> httpx.Response:
    return await client.post(
        f"/tasks/{task_id}/coverage-check",
        json={"reference_fasta": fasta, "depth_threshold": threshold},
    )


def result_file(task_id: str, side: int) -> Path:
    return config.TASKS_DIR / task_id / "output" / f"R{side}.fastq"


async def main_async() -> int:
    tmp_root = Path(tempfile.mkdtemp(prefix="fastq-coverage-",
                                     dir=str(ROOT / "tests" / "output")))
    data_dir = tmp_root / "data"
    os.environ["FASTQ_DATA_DIR"] = str(data_dir)
    config.DATA_DIR = data_dir
    config.DB_PATH = data_dir / "tasks.db"
    config.TASKS_DIR = data_dir / "tasks"

    db._initialized = False
    db.init_db()
    worker.recover_on_startup()

    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://test",
                                 timeout=120) as client:
        # ---- 1. happy path: one uniquely placed pair, exact 1-based intervals
        pairs = [(REF_RANDOM[100:160].encode(),
                  revcomp(REF_RANDOM[300:360]).encode())]
        body = await upload_pairs(client, pairs)
        record("coverage: fixture task succeeds", body["status"] == "succeeded",
               body.get("message", ""))
        tid = body["id"]
        hash_before = [hashlib.sha256(result_file(tid, s).read_bytes()).hexdigest()
                       for s in (1, 2)]

        r = await coverage(client, tid, fasta_text(REF_RANDOM), 1)
        record("coverage check returns 200", r.status_code == 200, str(r.status_code))
        got = r.json()
        record("counts: one placed pair",
               got["pairs_placed"] == 1 and got["pairs_no_placement"] == 0
               and got["pairs_multiple_placements"] == 0
               and got["pairs_total"] == 1, str(got))
        record("reference metadata echoed",
               got["reference_name"] == "ref1"
               and got["reference_length"] == 600
               and got["depth_threshold"] == 1, str(got))
        record("1-based closed low-coverage intervals",
               got["low_coverage_intervals"] == [
                   {"start": 1, "end": 100},
                   {"start": 161, "end": 300},
                   {"start": 361, "end": 600},
               ], str(got["low_coverage_intervals"]))
        r = await coverage(client, tid, fasta_text(REF_RANDOM), 2)
        record("depth-1 regions are below threshold 2",
               r.json()["low_coverage_intervals"] == [{"start": 1, "end": 600}],
               str(r.json()["low_coverage_intervals"]))
        hash_after = [hashlib.sha256(result_file(tid, s).read_bytes()).hexdigest()
                      for s in (1, 2)]
        record("result files untouched by the check", hash_before == hash_after)
        still = (await client.get(f"/tasks/{tid}")).json()
        record("task still succeeded after the check",
               still["status"] == "succeeded")

        # ---- 2. overlapping ends count twice
        pairs = [(REF_RANDOM[100:180].encode(),
                  revcomp(REF_RANDOM[150:230]).encode())]
        tid = (await upload_pairs(client, pairs))["id"]
        r = await coverage(client, tid, fasta_text(REF_RANDOM), 2)
        got = r.json()
        record("overlap of both ends counts twice",
               got["low_coverage_intervals"] == [
                   {"start": 1, "end": 150},
                   {"start": 181, "end": 600},
               ], str(got["low_coverage_intervals"]))

        # ---- 3. span constraint: only one of two R1 hits is compliant
        ref_span = "A" * 40 + "G" * 360 + "A" * 40 + "G" * 60 + "C" * 50 + "G" * 50
        pairs = [(b"A" * 40, b"G" * 50)]  # rc(R2) = C*50 at position 500
        tid = (await upload_pairs(client, pairs))["id"]
        r = await coverage(client, tid, fasta_text(ref_span), 1)
        got = r.json()
        record("span > 500 placements are not compliant",
               got["pairs_placed"] == 1
               and got["pairs_multiple_placements"] == 0, str(got))
        record("span-constrained placement intervals",
               got["low_coverage_intervals"] == [
                   {"start": 1, "end": 400},
                   {"start": 441, "end": 500},
                   {"start": 551, "end": 600},
               ], str(got["low_coverage_intervals"]))

        # ---- 4. repetitive hits -> multiple placements, no depth
        ref_rep = "A" * 200 + "C" * 200 + "G" * 200
        pairs = [(b"A" * 40, b"T" * 60)]
        tid = (await upload_pairs(client, pairs))["id"]
        r = await coverage(client, tid, fasta_text(ref_rep), 1)
        got = r.json()
        record("ambiguous pair counted as multiple placements",
               got["pairs_multiple_placements"] == 1 and got["pairs_placed"] == 0,
               str(got))
        record("ambiguous pair adds no depth",
               got["low_coverage_intervals"] == [{"start": 1, "end": 600}],
               str(got["low_coverage_intervals"]))

        # ---- 4b. span boundary: exactly 500 bp is compliant, 501 is not
        ref_500 = "A" * 40 + "G" * 400 + "C" * 60 + "G" * 50   # span 440+60-0=500
        tid = (await upload_pairs(client, [(b"A" * 40, b"G" * 60)]))["id"]
        r = await coverage(client, tid, fasta_text(ref_500), 1)
        got = r.json()
        record("span of exactly 500 bp is compliant",
               got["pairs_placed"] == 1
               and got["low_coverage_intervals"] == [
                   {"start": 41, "end": 440},
                   {"start": 501, "end": 550},
               ], str(got))
        ref_501 = "A" * 40 + "G" * 401 + "C" * 60 + "G" * 49   # span 501
        tid = (await upload_pairs(client, [(b"A" * 40, b"G" * 60)]))["id"]
        r = await coverage(client, tid, fasta_text(ref_501), 1)
        got = r.json()
        record("span of 501 bp is not compliant",
               got["pairs_no_placement"] == 1 and got["pairs_placed"] == 0,
               str(got))

        # ---- 4c. identical start positions are not "increasing"
        ref_same = "AC" * 30 + "T" * 40
        pairs = [(b"AC" * 30, b"GT" * 30)]  # rc(R2) = (AC)*30, same start 0
        tid = (await upload_pairs(client, pairs))["id"]
        r = await coverage(client, tid, fasta_text(ref_same), 1)
        got = r.json()
        record("equal starts are not an increasing placement",
               got["pairs_no_placement"] == 1 and got["pairs_placed"] == 0,
               str(got))

        # ---- 5. no compliant placement
        pairs = [(b"A" * 40, b"T" * 40)]
        tid = (await upload_pairs(client, pairs))["id"]
        r = await coverage(client, tid, fasta_text(REF_RANDOM), 3)
        got = r.json()
        record("unmatched pair counted as no placement",
               got["pairs_no_placement"] == 1 and got["pairs_placed"] == 0
               and got["low_coverage_intervals"] == [{"start": 1, "end": 600}],
               str(got))

        # ---- 6. one pair of each kind in a single run
        ref_mix = REF_RANDOM + "A" * 200
        pairs = [
            (REF_RANDOM[100:160].encode(), revcomp(REF_RANDOM[300:360]).encode()),
            (b"C" * 40, b"G" * 40),
            (b"A" * 40, b"T" * 60),
        ]
        tid = (await upload_pairs(client, pairs))["id"]
        r = await coverage(client, tid, fasta_text(ref_mix), 1)
        got = r.json()
        record("three pair categories counted separately",
               (got["pairs_placed"], got["pairs_no_placement"],
                got["pairs_multiple_placements"], got["pairs_total"]) == (1, 1, 1, 3),
               str(got))
        record("only the placed pair contributes depth",
               got["low_coverage_intervals"] == [
                   {"start": 1, "end": 100},
                   {"start": 161, "end": 300},
                   {"start": 361, "end": 800},
               ], str(got["low_coverage_intervals"]))

        # ---- 7. depth threshold boundaries
        r = await coverage(client, tid, fasta_text(ref_mix), 0)
        record("threshold 0 -> 422", r.status_code == 422, str(r.status_code))
        r = await coverage(client, tid, fasta_text(ref_mix), 21)
        record("threshold 21 -> 422", r.status_code == 422, str(r.status_code))
        r = await coverage(client, tid, fasta_text(ref_mix), 20)
        got = r.json()
        record("threshold 20 accepted, whole reference below it",
               r.status_code == 200
               and got["low_coverage_intervals"] == [{"start": 1, "end": 800}],
               str(got["low_coverage_intervals"]))

        # ---- 8. reference FASTA validation
        bad_fastas = {
            "missing header": "ACGTACGT\n",
            "two sequences": ">r\nACGT\n>r2\nACGT\n",
            "non-ACGT base": ">r\nACGTNACGT\n",
            "over 10000 bp": ">r\n" + "A" * 10001 + "\n",
            "empty input": "",
            "no sequence": ">r\n\n\n",
        }
        all_422 = True
        details = []
        for label, bad in bad_fastas.items():
            rr = await coverage(client, tid, bad, 5)
            details.append(f"{label}={rr.status_code}")
            all_422 = all_422 and rr.status_code == 422
        record("invalid reference FASTA -> 422", all_422, ", ".join(details))

        r = await coverage(client, tid, fasta_text(REF_RANDOM.lower() + "a" * 200), 1)
        record("lowercase ACGT reference accepted",
               r.status_code == 200
               and r.json()["low_coverage_intervals"][-1] == {"start": 361, "end": 800},
               str(r.status_code))
        wrapped = fasta_text(REF_RANDOM[:300] + "\n" + REF_RANDOM[300:] + "\n" + "A" * 200)
        r = await coverage(client, tid, wrapped, 1)
        record("multi-line wrapped FASTA accepted",
               r.status_code == 200
               and r.json()["low_coverage_intervals"][0] == {"start": 1, "end": 100},
               str(r.status_code))

        # ---- 9. task state errors
        r = await coverage(client, uuid.uuid4().hex, fasta_text(REF_RANDOM), 5)
        record("unknown task -> 404", r.status_code == 404, str(r.status_code))

        bad_r2 = fastq([make_record("x", b"AAAA", b"~~~")])
        ok_r1 = fastq([make_record("x", b"AAAA")])
        resp = await client.post(
            "/tasks",
            files={"r1": ("r1.fastq", ok_r1, "application/octet-stream"),
                   "r2": ("r2.fastq", bad_r2, "application/octet-stream")})
        failed_body = await wait_terminal(client, resp.json()["id"])
        r = await coverage(client, failed_body["id"], fasta_text(REF_RANDOM), 5)
        record("failed task -> 409",
               failed_body["status"] == "failed" and r.status_code == 409,
               f"task={failed_body['status']} http={r.status_code}")

        tid = (await upload_pairs(client, [(REF_RANDOM[0:60].encode(),
                                            revcomp(REF_RANDOM[100:160]).encode())]))["id"]
        result_file(tid, 1).unlink()
        r = await coverage(client, tid, fasta_text(REF_RANDOM), 5)
        record("missing result file -> error, no partial intervals",
               r.status_code >= 500, str(r.status_code))

        # ---- 10. result-file content errors fail the whole check (422)
        async def corrupted_task(r1_bytes: bytes | None,
                                 r2_bytes: bytes | None) -> str:
            pairs = [(REF_RANDOM[0:60].encode(), revcomp(REF_RANDOM[100:160]).encode())]
            body = await upload_pairs(client, pairs)
            cid = body["id"]
            if r1_bytes is not None:
                result_file(cid, 1).write_bytes(r1_bytes)
            if r2_bytes is not None:
                result_file(cid, 2).write_bytes(r2_bytes)
            return cid

        good1 = make_record("pair-0", REF_RANDOM[0:60].encode())
        good2 = make_record("pair-0", revcomp(REF_RANDOM[100:160]).encode())
        cases = {
            "read id mismatch": (None, make_record("other-id", b"A" * 60)),
            "seq/qual length mismatch": (make_record("pair-0", b"A" * 60, b"~" * 59), None),
            "header without @": (b"pair-0\n" + b"A" * 60 + b"\n+\n" + b"~" * 60 + b"\n", None),
            "truncated record": (b"@pair-0\n" + b"A" * 60 + b"\n+\n", None),
            "record count mismatch": (None, good2 + make_record("extra", b"T" * 60)),
        }
        all_422 = True
        details = []
        for label, (r1b, r2b) in cases.items():
            cid = await corrupted_task(r1b, r2b)
            rr = await coverage(client, cid, fasta_text(REF_RANDOM), 5)
            details.append(f"{label}={rr.status_code}")
            all_422 = all_422 and rr.status_code == 422
        record("corrupt result files -> 422, no partial intervals",
               all_422, ", ".join(details))

        short = await upload_pairs(client, [(b"A" * 30, b"T" * 30)])
        r = await coverage(client, short["id"], fasta_text(REF_RANDOM), 5)
        record("read length outside 40-100 bp -> 422",
               r.status_code == 422, str(r.status_code))

        too_many = await upload_pairs(client, [(b"A" * 40, b"T" * 40)] * 2001)
        r = await coverage(client, too_many["id"], fasta_text(REF_RANDOM), 5)
        record("more than 2000 pairs -> 422", r.status_code == 422,
               str(r.status_code))
        exactly = await upload_pairs(client, [(b"A" * 40, b"T" * 40)] * 2000)
        r = await coverage(client, exactly["id"], fasta_text(REF_RANDOM), 5)
        got = r.json()
        record("exactly 2000 pairs accepted",
               r.status_code == 200 and got["pairs_no_placement"] == 2000,
               f"http={r.status_code} no_placement={got.get('pairs_no_placement')}")

        # ---- 11. repetitive reference stress: 10k homopolymer, 2000 pairs
        ref_homo = "A" * 10000
        stress_pairs = [
            (b"A" * (40 + i % 61), b"T" * (40 + (i * 13) % 61))
            for i in range(2000)
        ]
        body = await upload_pairs(client, stress_pairs)
        record("stress fixture task succeeds", body["status"] == "succeeded",
               body.get("message", ""))
        started = time.monotonic()
        r = await coverage(client, body["id"], fasta_text(ref_homo), 5)
        elapsed = time.monotonic() - started
        got = r.json()
        record("repetitive stress completes quickly",
               r.status_code == 200 and elapsed < 30, f"{elapsed:.2f}s")
        record("repetitive pairs all counted as multiple placements",
               got["pairs_multiple_placements"] == 2000
               and got["pairs_total"] == 2000, str(got))
        record("repetitive stress intervals",
               got["low_coverage_intervals"] == [{"start": 1, "end": 10000}],
               str(got["low_coverage_intervals"]))

        # ---- 12. fully covered reference -> no low-coverage intervals
        cover_pairs = [
            (REF_RANDOM[s:s + 60].encode(), revcomp(REF_RANDOM[s + 60:s + 120]).encode())
            for s in range(0, 600, 120)
        ]
        tid = (await upload_pairs(client, cover_pairs))["id"]
        r = await coverage(client, tid, fasta_text(REF_RANDOM), 1)
        got = r.json()
        record("fully covered reference -> empty interval list",
               r.status_code == 200 and got["pairs_placed"] == 5
               and got["low_coverage_intervals"] == [], str(got))

        # ---- 13. interactive docs expose the new endpoint
        r = await client.get("/docs")
        record("/docs reachable", r.status_code == 200, str(r.status_code))
        r = await client.get("/openapi.json")
        record("openapi lists the coverage endpoint",
               r.status_code == 200
               and "/tasks/{task_id}/coverage-check" in r.json().get("paths", {}),
               str(r.status_code))

    shutil.rmtree(tmp_root, ignore_errors=True)
    failed = [name for name, ok, _ in RESULTS if not ok]
    print("\n%d/%d checks passed" % (len(RESULTS) - len(failed), len(RESULTS)))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(asyncio.run(main_async()))
