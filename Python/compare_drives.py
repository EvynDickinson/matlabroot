#!/usr/bin/env python3
"""
compare_drives.py — Deep folder comparison for drive verification
Checks for missing files, size mismatches, and MD5/SHA256 corruption

** GENERALLY RUN THIS ONE: **
python compare_drives.py --hash md5 --workers 2

* Run this line to set a custom output location for the final report *
python compare_drives.py /source /dest --out ~/Desktop/my_report

* Run this way if you want to use parallel computing powers *
python compare_drives.py /source /dest --workers 4 --out ~/Desktop/my_report

* Run this way for most sophisticated check (includes integrity check) *
python compare_drives.py /source /dest --hash md5 --out ~/Desktop/my_report
python compare_drives.py --hash md5

"""
#!/usr/bin/env python3
"""
================================================================================
compare_drives.py — Deep Folder Comparison Tool for Drive Verification
================================================================================

PURPOSE:
    This script compares the complete contents of two directories (recursively,
    meaning it descends into every subfolder, and every subfolder of that, all
    the way down) and produces a detailed report of any differences it finds.

    It is designed to answer the question: "Did my files copy correctly?"

    Specifically, it detects:
      1. FILES MISSING IN DESTINATION — files that exist on the source drive
         but did not make it to the destination drive at all. This catches
         incomplete copies where the transfer was interrupted or skipped files.

      2. EXTRA FILES IN DESTINATION — files that exist on the destination but
         NOT on the source. This is usually harmless (e.g. the destination had
         pre-existing files) but worth knowing about.

      3. SIZE MISMATCHES — files that exist on both drives but have different
         byte counts. This almost always indicates a corrupted or truncated copy,
         where the transfer stopped partway through writing a file.

      4. CHECKSUM MISMATCHES — files that have the exact same size on both drives
         but whose contents differ. This is the sneakiest kind of corruption: a
         file can be the same number of bytes yet contain completely different
         (wrong) data due to bit-flips, write errors, or filesystem bugs. Only a
         cryptographic hash can catch this.

      5. UNREADABLE FILES — files that exist but could not be opened or read,
         usually because of permission errors or filesystem damage.

HOW IT WORKS (high level):
    Phase 1 — SCAN SOURCE: Walk every file in the source directory tree.
              For each file, record its relative path, size, and modification
              time. Optionally compute a cryptographic hash of its contents.

    Phase 2 — SCAN DESTINATION: Do exactly the same for the destination tree.

    Phase 3 — COMPARE: For every file in the source, check whether it exists
              in the destination and whether the size/hash match. Also flag any
              files that are only in the destination.

    Phase 4 — REPORT: Print a human-readable summary to the terminal and write
              detailed JSON and CSV report files to disk.

HASH ALGORITHMS SUPPORTED:
    md5     — Fast. Produces a 128-bit hash. More than sufficient for detecting
               accidental corruption (not for security purposes).
    sha1    — Slightly slower. 160-bit hash. Rarely needed over md5 for this use.
    sha256  — Slower but produces a 256-bit hash. Use if you want maximum
               confidence or are paranoid about hash collisions (extremely rare).
    sha512  — Slowest. 512-bit hash. Overkill for drive verification but here
               if you want it.
    none    — Skip hashing entirely. Only size and existence are checked. Very
               fast but will NOT catch same-size corruption.

PERFORMANCE NOTES:
    - Hashing is the bottleneck. A spinning hard drive can do ~100-150 MB/s;
      an SSD ~500 MB/s; NVMe ~3000 MB/s. A 1 TB drive at 150 MB/s takes ~2 hrs.
    - The --workers flag controls how many files are hashed in parallel. More
      workers help on SSDs/NVMe but can actually slow down spinning drives due
      to seek thrashing. Default of 4 is a safe middle ground.

OUTPUT FILES:
    drive_compare_YYYYMMDD_HHMMSS.json  — Full machine-readable report
    drive_compare_YYYYMMDD_HHMMSS.csv   — Spreadsheet-friendly list of issues

EXIT CODES:
    0 — No issues found. The copy is verified complete and intact.
    1 — One or more issues were detected. Check the report.

USAGE EXAMPLES:
    # Basic check with MD5 hashing (recommended default)
    python compare_drives.py /Volumes/SourceDrive /Volumes/BackupDrive

    # Maximum integrity check with SHA256
    python compare_drives.py /source /dest --hash sha256

    # Fast size-only check (no hash, will miss same-size corruption)
    python compare_drives.py /source /dest --hash none

    # Save reports to a specific location
    python compare_drives.py /source /dest --out ~/Desktop/verification_report

    # Use 8 parallel workers (good for NVMe drives)
    python compare_drives.py /source /dest --workers 8

REQUIREMENTS:
    Python 3.6 or newer. No third-party packages required — only standard
    library modules are used (os, sys, hashlib, argparse, json, csv, pathlib,
    datetime, concurrent.futures, dataclasses, typing).

AUTHOR NOTE:
    This script is read-only with respect to both drives. It will never write,
    move, delete, or modify any file on either the source or destination.
    It only reads files and writes its own report files to the current directory
    (or wherever --out points).
================================================================================
"""

# ── Standard library imports ──────────────────────────────────────────────────
# os: Used for os.walk(), which is the function that recursively traverses a
#     directory tree yielding (dirpath, subdirs, filenames) tuples for each
#     directory it visits.
import os

# sys: Used for sys.exit(), which terminates the script with a specific exit
#     code (0 = success, 1 = issues found). Also used for sys.argv implicitly
#     via argparse.
import sys

# hashlib: Python's built-in cryptographic hashing library. Provides MD5, SHA1,
#     SHA256, SHA512, and many others. We use hashlib.new(algorithm) to create
#     a hasher object dynamically based on the user's --hash argument.
import hashlib

# argparse: The standard library module for parsing command-line arguments.
#     It automatically generates --help text, validates argument types, and
#     makes positional vs optional arguments easy to define.
import argparse

# json: Used to serialize the comparison results into a human-readable and
#     machine-parseable JSON report file. json.dump() writes Python dicts/lists
#     to a file with optional pretty-printing (indent=2).
import json

# csv: Used to write a comma-separated values report that can be opened directly
#     in Excel, Numbers, Google Sheets, etc. csv.DictWriter makes it easy to
#     write rows from Python dictionaries.
import csv

# pathlib.Path: A modern, object-oriented way to work with filesystem paths.
#     Path objects know how to join paths, resolve symlinks, check existence,
#     call .stat() for file metadata, etc. Much cleaner than os.path string
#     manipulation.
from pathlib import Path

# datetime: Used to record when the comparison was run (for the report header)
#     and to calculate how long the whole process took (elapsed time).
from datetime import datetime

# concurrent.futures: Provides ThreadPoolExecutor, which manages a pool of
#     worker threads. We submit file-hashing tasks to the pool so multiple files
#     can be hashed concurrently. as_completed() lets us process results as each
#     thread finishes rather than waiting for all of them.
from concurrent.futures import ThreadPoolExecutor, as_completed

# dataclasses: The @dataclass decorator automatically generates __init__,
#     __repr__, and other boilerplate methods for classes that are primarily
#     just containers for data. field(default_factory=list) creates a fresh
#     empty list for each instance rather than sharing one mutable default.
from dataclasses import dataclass, field, asdict

# typing.Optional: Type hint meaning "this value is either the specified type
#     OR None". Used to indicate that checksum may not be computed (e.g. when
#     --hash none is passed, or if the file couldn't be read).
from typing import Optional

import sys
sys.stdout.reconfigure(encoding="utf-8")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1: TERMINAL OUTPUT HELPERS
# ══════════════════════════════════════════════════════════════════════════════
#
# ANSI escape codes are special character sequences that terminal emulators
# interpret as formatting instructions rather than printing them literally.
# The format is: ESC[ + code + m
#   ESC is the escape character, represented in Python as \033 (octal) or \x1b (hex)
#   The '[' is the "Control Sequence Introducer"
#   The number is the formatting code
#   'm' ends the sequence
#
# Code reference:
#   92 = bright green foreground
#   93 = bright yellow foreground
#   91 = bright red foreground
#   96 = bright cyan foreground
#    1 = bold text
#    0 = reset all formatting (return to default)
#
# We define these as module-level constants so they can be reused everywhere
# without repeating the magic strings. They're intentionally short names
# because they're used inline inside f-strings constantly.

GREEN  = "\033[92m"   # Bright green  — used for success/OK messages
YELLOW = "\033[93m"   # Bright yellow — used for warnings (non-fatal issues)
RED    = "\033[91m"   # Bright red    — used for errors (files missing/corrupted)
CYAN   = "\033[96m"   # Bright cyan   — used for informational/progress messages
BOLD   = "\033[1m"    # Bold text     — used for headers and section titles
RESET  = "\033[0m"    # Reset all     — MUST be placed after every coloured string
                      #                 to stop the colour bleeding into the next line


# ── Convenience print functions ───────────────────────────────────────────────
# These four one-liners wrap print() with consistent formatting so every
# category of message looks the same throughout the script. The icon prefix
# (✔ ⚠ ✘ ·) gives a quick visual scan of the output even without colour support.
#
# Usage:
#   ok("Everything matched")      →  green   ✔  Everything matched
#   warn("3 extra files found")   →  yellow  ⚠  3 extra files found
#   err("12 files missing")       →  red     ✘  12 files missing
#   info("Scanning source …")     →  cyan    ·  Scanning source …

def ok(msg):   print(f"{GREEN}  ✔  {msg}{RESET}")   # Success — all good
def warn(msg): print(f"{YELLOW}  ⚠  {msg}{RESET}")  # Warning — noteworthy but not fatal
def err(msg):  print(f"{RED}  ✘  {msg}{RESET}")     # Error — something is wrong
def info(msg): print(f"{CYAN}  ·  {msg}{RESET}")    # Info — general progress update

#def ok(msg):   print(f"{GREEN}  [OK]  {msg}{RESET}")
#def warn(msg): print(f"{YELLOW}  [!!]  {msg}{RESET}")
#def err(msg):  print(f"{RED}  [XX]  {msg}{RESET}")
#def info(msg): print(f"{CYAN}  [..]  {msg}{RESET}")

# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2: DATA STRUCTURES
# ══════════════════════════════════════════════════════════════════════════════
#
# We define two dataclasses to hold our data cleanly. Using dataclasses instead
# of plain dicts gives us:
#   - Type hints (self-documenting what each field contains)
#   - Auto-generated __init__ (no boilerplate constructor needed)
#   - Auto-generated __repr__ (useful for debugging)
#   - Immutability option if needed later


@dataclass
class FileStat:
    """
    Holds all the metadata we collect about a single file during the scan phase.

    This is created once per file, for both the source and destination scans.
    After both scans complete, we pair up FileStat objects by their rel_path
    and compare their fields.

    Fields:
        rel_path (str):
            The file's path RELATIVE to the scan root. For example, if we're
            scanning /Volumes/Drive and find /Volumes/Drive/Photos/2023/IMG001.jpg,
            the rel_path is "Photos/2023/IMG001.jpg". We use relative paths (not
            absolute) as the dictionary key so we can compare source and destination
            files that live at different absolute locations but the same relative
            position within their respective trees.

        size (int):
            The file's size in bytes, as reported by the filesystem. This is the
            first and fastest check — if sizes differ, the files are definitely
            different and we don't need to bother hashing.

        mtime (float):
            The file's last-modification timestamp as a Unix epoch float (seconds
            since Jan 1 1970 UTC). We collect this but don't currently use it in
            comparisons — it's stored in case you want to extend the script to
            detect files that were re-copied (same size/hash, different timestamp).

        checksum (Optional[str]):
            The hex-encoded hash digest of the file's full contents. For example,
            an MD5 checksum looks like "d41d8cd98f00b204e9800998ecf8427e".
            This is None if --hash none was specified OR if a read error occurred.

        checksum_error (Optional[str]):
            If we tried to hash the file but hit an exception (e.g. permission
            denied, file disappeared mid-scan), we store the error message here
            instead of crashing. This is None if hashing succeeded or was skipped.
    """
    rel_path:        str            # Relative path from the scan root (used as dict key)
    size:            int            # File size in bytes
    mtime:           float          # Last-modified timestamp (Unix epoch)
    checksum:        Optional[str] = None  # Hex hash digest, or None if not computed
    checksum_error:  Optional[str] = None  # Error message if hashing failed, else None


@dataclass
class CompareResult:
    """
    Accumulates all the findings from comparing the source and destination scans.

    This object is populated by the compare() function and then consumed by
    the print_report(), write_json_report(), and write_csv_report() functions.

    The lists contain either plain strings (for simple missing-file cases) or
    dicts with extra detail fields (for mismatches where we want to show both
    the source and destination values).

    Fields:
        missing_in_dest (list[str]):
            Relative paths of files that exist in the source scan but have NO
            corresponding entry in the destination scan. These are files that
            simply didn't get copied over. The list is sorted alphabetically.

        missing_in_src (list[str]):
            Relative paths of files that exist in the destination scan but NOT
            in the source scan. These are "extra" files — the destination has
            something the source doesn't. Usually harmless but reported for
            completeness. The list is sorted alphabetically.

        size_mismatches (list[dict]):
            Files present on both drives with DIFFERENT sizes. Each entry is a
            dict with keys: path, src_size, dst_size, diff_bytes.
            A size mismatch almost certainly means the copy was truncated or
            the file was modified after copying. We stop checking this file
            further (no point hashing if sizes differ).

        checksum_mismatches (list[dict]):
            Files present on both drives with the SAME size but DIFFERENT hash.
            Each entry is a dict with keys: path, src_checksum, dst_checksum, size.
            This is the most insidious kind of corruption: the byte count is
            identical but the actual bits differ. Could be caused by:
              - RAM bit-flip during copy (rare but real on servers)
              - Failing storage sector that reads as zeros
              - Filesystem driver bug
              - File was modified after the copy started

        checksum_errors (list[dict]):
            Files we couldn't hash due to I/O errors. Each entry is a dict with
            keys: path, src_error, dst_error. At least one of the error fields
            will be a non-None string describing what went wrong.

        ok_files (int):
            Count of files that passed ALL checks (exist on both sides, same size,
            same hash). This is the "good news" counter.

        total_src_files (int):
            Total number of files found in the source tree. Set at the start of
            compare() from len(src_stats).

        total_dst_files (int):
            Total number of files found in the destination tree.

        src_total_bytes (int):
            Sum of all file sizes in the source tree (used for display only).

        dst_total_bytes (int):
            Sum of all file sizes in the destination tree.

        elapsed_seconds (float):
            How long the entire operation took from start to finish, in seconds.
            Set in main() after all phases complete.
    """
    missing_in_dest:     list  = field(default_factory=list)  # In src, not in dst
    missing_in_src:      list  = field(default_factory=list)  # In dst, not in src (extras)
    size_mismatches:     list  = field(default_factory=list)  # Same path, different byte count
    checksum_mismatches: list  = field(default_factory=list)  # Same size, different hash digest
    checksum_errors:     list  = field(default_factory=list)  # Files that couldn't be read/hashed
    ok_files:            int   = 0      # Count of files that passed all checks
    total_src_files:     int   = 0      # Total files found in source tree
    total_dst_files:     int   = 0      # Total files found in destination tree
    src_total_bytes:     int   = 0      # Sum of all source file sizes in bytes
    dst_total_bytes:     int   = 0      # Sum of all destination file sizes in bytes
    elapsed_seconds:     float = 0.0   # Wall-clock time for the whole run


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3: FILE SCANNING
# ══════════════════════════════════════════════════════════════════════════════


def compute_checksum(path: Path, algorithm: str = "md5", chunk: int = 1 << 20) -> str:
    """
    Read a file and compute a cryptographic hash of its entire contents.

    This function is the heart of corruption detection. Rather than loading the
    whole file into memory at once (which would be catastrophic for a 50GB video
    file), we read it in chunks and feed each chunk to the hasher incrementally.
    The hasher internally accumulates a running state so the final digest is the
    same as if we'd fed it the whole file at once.

    WHY HASHING WORKS FOR CORRUPTION DETECTION:
        A good hash function has the "avalanche effect": changing even a single
        bit anywhere in the input produces a completely different hash output.
        So if even one byte of a file got corrupted during the copy, the hash
        of the copy will be wildly different from the hash of the original.
        The probability of two different files having the same MD5 hash is
        approximately 1 in 2^128 — astronomically unlikely by accident.

    Args:
        path (Path):
            Absolute path to the file to hash. Must be readable.

        algorithm (str):
            Name of the hash algorithm to use. Must be a name recognised by
            hashlib.new() — e.g. "md5", "sha1", "sha256", "sha512".
            Defaults to "md5" which is fast and collision-resistant enough for
            data integrity (not security) purposes.

        chunk (int):
            Number of bytes to read per iteration. Default is 1 << 20 = 1,048,576
            bytes = 1 MiB. This is a good balance between:
              - Too small: many read() syscalls, slow due to overhead
              - Too large: excessive RAM usage for large files
            You can increase this (e.g. 1 << 23 = 8 MiB) on systems with lots
            of RAM and fast storage for a small speed improvement.

    Returns:
        str: The hexadecimal string representation of the digest.
             MD5 produces 32 hex chars (128 bits).
             SHA256 produces 64 hex chars (256 bits).

    Raises:
        OSError: If the file cannot be opened or read (permission denied,
                 file disappeared, I/O error on the drive, etc.).
                 The caller (scan_tree's _hash inner function) catches this
                 and stores the error message in FileStat.checksum_error.
    """
    # hashlib.new() creates a hasher object for the named algorithm.
    # Using new() instead of e.g. hashlib.md5() allows the algorithm to be
    # chosen dynamically at runtime from the user's --hash argument.
    h = hashlib.new(algorithm)

    # Open the file in binary mode ("rb").
    # Binary mode is essential — we must read the raw bytes exactly as they
    # are stored on disk, with no newline translation or encoding applied.
    with open(path, "rb") as f:
        # Read the file in chunks until we've consumed it all.
        # f.read(chunk) returns an empty bytes object b"" when EOF is reached,
        # which is falsy in Python, so `if not data: break` exits cleanly.
        while True:
            data = f.read(chunk)
            if not data:
                break  # EOF reached — we've processed the entire file

            # Feed this chunk into the hasher. The hasher updates its internal
            # state to incorporate these bytes. We can call update() as many
            # times as we want; the result is identical to calling it once with
            # all the bytes concatenated.
            h.update(data)

    # hexdigest() returns the final hash as a lowercase hex string.
    # e.g. "d41d8cd98f00b204e9800998ecf8427e" for MD5 of an empty file.
    return h.hexdigest()


def scan_tree(root: Path, algorithm: str, workers: int, verbose: bool) -> dict[str, FileStat]:
    """
    Recursively scan an entire directory tree and return a dictionary of file metadata.

    This function has two phases:
      1. WALK PHASE: Use os.walk() to discover every file in the tree and collect
         basic filesystem metadata (size, mtime) via stat(). This is fast because
         it only reads directory entries, not file contents.

      2. HASH PHASE: Submit all files to a thread pool for parallel hashing.
         Hashing reads every byte of every file, which is I/O-intensive and
         benefits from parallelism on multi-drive or SSD setups.

    Args:
        root (Path):
            The top-level directory to scan. Can be relative or absolute —
            we call root.resolve() at the start to normalise it to an absolute
            path, which is necessary for computing correct relative paths later.

        algorithm (str):
            Hash algorithm to use, forwarded to compute_checksum(). If "none",
            hashing is skipped entirely and checksum fields remain None.

        workers (int):
            Maximum number of threads in the thread pool. More threads = more
            files hashed simultaneously. Optimal value depends on storage type:
              - Spinning HDD: 1-2 (seek thrashing kills performance with more)
              - SATA SSD: 4-8
              - NVMe SSD: 8-16
              - Network drive: depends on network latency and server capacity

        verbose (bool):
            If True, print a progress update after every single file rather than
            only at every 5% milestone. Useful for debugging slow scans.

    Returns:
        dict[str, FileStat]:
            A dictionary mapping each file's relative path string to its FileStat.
            Example key: "Documents/Work/report.pdf"
            Example value: FileStat(rel_path="Documents/Work/report.pdf",
                                    size=204800, mtime=1699123456.0,
                                    checksum="a3f5...", checksum_error=None)
            The relative path is used as the key because it's the only thing
            that's meaningful to compare between source and destination
            (absolute paths will differ between the two drives).
    """
    # Resolve the root to an absolute, symlink-free path.
    # This is necessary so that Path.relative_to(root) works correctly later
    # when computing each file's relative path.
    root = root.resolve()

    # ── Phase 1: Walk the directory tree ─────────────────────────────────────
    # os.walk() is a generator that yields a 3-tuple for each directory:
    #   dirpath:   the absolute path of the current directory (as a string)
    #   dirnames:  list of subdirectory names in dirpath (we ignore this — walk
    #              will recurse into them automatically on the next iteration)
    #   filenames: list of file names in dirpath (NOT full paths, just names)
    #
    # We iterate over this to build a flat list of all files with their metadata.
    # We collect everything into a list first (rather than hashing inline) so
    # we know the total file count before starting the thread pool, which lets
    # us show accurate percentage progress.
    paths = []  # Will hold tuples of (rel_path_str, full_path, size_int, mtime_float)

    for dirpath, _, filenames in os.walk(root):
        for fn in filenames:
            # Build the full absolute path by joining the directory path with
            # the filename. We use Path() to get cross-platform path handling.
            full = Path(dirpath) / fn

            # Compute the relative path by stripping the root prefix.
            # e.g. if root is /Volumes/Drive and full is /Volumes/Drive/a/b.txt
            # then rel is "a/b.txt"
            # We convert to str immediately because we use it as a dict key.
            rel = str(full.relative_to(root))

            try:
                # stat() makes a single syscall to get file metadata from the
                # filesystem without reading the file's contents. Very fast.
                # st_size: file size in bytes
                # st_mtime: last modification time as Unix timestamp float
                st = full.stat()
                paths.append((rel, full, st.st_size, st.st_mtime))
            except OSError as e:
                # stat() can fail if the file is a broken symlink, if we don't
                # have read permission on the parent directory, or if the file
                # disappeared between os.walk() finding it and us stat()-ing it
                # (e.g. another process deleted it). We warn and skip rather
                # than crashing the whole scan.
                warn(f"Cannot stat {rel}: {e}")

    # ── Phase 2: Hash files in parallel ──────────────────────────────────────
    stats: dict[str, FileStat] = {}  # The result dict we'll return
    total = len(paths)               # Total file count (for progress display)

    def _hash(item):
        """
        Inner function that hashes a single file. This is what each worker
        thread executes. It's defined as a closure inside scan_tree so it has
        access to the `algorithm` variable from the outer scope.

        Args:
            item (tuple): A (rel, full, size, mtime) tuple from the paths list.

        Returns:
            tuple: (rel_path_str, FileStat) — the key-value pair to store in stats.

        Note on thread safety:
            Each call to _hash works on a different file and creates its own
            FileStat object. There's no shared mutable state between threads,
            so no locks are needed. The stats dict is only written to in the
            main thread (in the as_completed loop below).
        """
        rel, full, size, mtime = item

        # Create a FileStat with the metadata we already have from the stat() call.
        # checksum and checksum_error start as None (the dataclass defaults).
        fs = FileStat(rel_path=rel, size=size, mtime=mtime)

        # Only attempt hashing if the user didn't pass --hash none.
        if algorithm != "none":
            try:
                # compute_checksum() reads the entire file and returns a hex digest.
                # This is the slow part — for a 1 GB file it might take several seconds.
                fs.checksum = compute_checksum(full, algorithm)
            except OSError as e:
                # Store the error message in the FileStat instead of crashing.
                # The compare() function will check for this and report it separately.
                fs.checksum_error = str(e)

        return rel, fs

    # ThreadPoolExecutor manages a pool of daemon threads.
    # max_workers controls how many threads run simultaneously.
    # The `with` statement ensures all threads are properly joined (waited for)
    # when the block exits, even if an exception occurs.
    done = 0  # Counter for progress display — incremented as each future completes
    with ThreadPoolExecutor(max_workers=workers) as ex:
        # Submit all files to the thread pool at once. ex.submit(fn, arg) schedules
        # fn(arg) to run on a worker thread and returns a Future object immediately
        # without blocking. We store the futures in a dict mapping future->item so
        # we could look up which item each future corresponds to (useful for error
        # messages, though we don't currently use it beyond what _hash returns).
        futs = {ex.submit(_hash, item): item for item in paths}

        # as_completed() yields futures one by one as they finish, in the order
        # they complete (which is generally NOT the order they were submitted).
        # This lets us process results and update progress as fast as threads finish,
        # rather than waiting for all threads to complete before doing anything.
        for fut in as_completed(futs):
            # fut.result() retrieves the return value of _hash(). If _hash raised
            # an uncaught exception, fut.result() would re-raise it here. We rely
            # on _hash catching its own OSErrors internally, so this should always
            # succeed.
            rel, fs = fut.result()

            # Store the completed FileStat in our result dict.
            # This is safe because only the main thread writes to `stats`;
            # worker threads only read their input and return values.
            stats[rel] = fs
            done += 1

            # Progress reporting: print a \r (carriage return without newline)
            # to overwrite the current line. We show progress either on every file
            # (verbose mode) or at every 5% milestone (total // 20 files).
            # max(1, ...) guards against division by zero when total is 0.
            if verbose or done % max(1, total // 20) == 0:
                pct = done * 100 // total if total else 100
                # end="" suppresses the automatic newline; flush=True forces the
                # output to appear immediately rather than being buffered.
                print(f"\r  Scanning … {done}/{total}  ({pct}%)", end="", flush=True)

    # Print a newline to end the \r progress line cleanly before the next output.
    print()
    return stats


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4: COMPARISON LOGIC
# ══════════════════════════════════════════════════════════════════════════════


def compare(src_stats: dict, dst_stats: dict) -> CompareResult:
    """
    Compare two scan result dictionaries and categorise every difference found.

    This function takes the output of two scan_tree() calls (one for source,
    one for destination) and produces a CompareResult describing all discrepancies.

    The comparison logic uses Python's set operations on the dictionary keys
    (which are relative file paths) to efficiently find which files are unique
    to each side, then iterates over the common files to check their metadata.

    COMPARISON PRIORITY (for files that exist on both sides):
        We check conditions in a specific order and `continue` (skip to the
        next file) as soon as we find a problem. This means:
          1. If there's a read error -> log it, stop checking this file.
          2. If sizes differ -> log it, stop checking (no point hashing).
          3. If checksums differ -> log it, stop checking.
          4. Otherwise -> increment ok_files counter.
        A file can only appear in ONE of the result lists (the first problem wins).

    Args:
        src_stats (dict[str, FileStat]):
            Output of scan_tree() for the source directory. Keys are relative
            paths, values are FileStat objects.

        dst_stats (dict[str, FileStat]):
            Output of scan_tree() for the destination directory. Same structure.

    Returns:
        CompareResult: A fully populated result object. See the CompareResult
                       docstring for a description of each field.
    """
    # Initialise the result object. We pre-populate the totals from the scan
    # results immediately, since these don't require any per-file comparison.
    r = CompareResult(
        total_src_files=len(src_stats),
        total_dst_files=len(dst_stats),
        # Generator expressions: sum file sizes across all FileStat values.
        # This gives us the total data volume on each drive for the report header.
        src_total_bytes=sum(f.size for f in src_stats.values()),
        dst_total_bytes=sum(f.size for f in dst_stats.values()),
    )

    # Convert dict keys to sets so we can use Python's set algebra operators.
    # Set operations are O(n) and extremely fast — finding the difference between
    # two sets of 500,000 file paths takes milliseconds.
    src_keys = set(src_stats)  # Set of all relative paths found in source
    dst_keys = set(dst_stats)  # Set of all relative paths found in destination

    # SET DIFFERENCE: src_keys - dst_keys gives us paths that are in the source
    # but NOT in the destination. These are the missing files.
    # sorted() makes the output deterministic and easier to read in the report.
    r.missing_in_dest = sorted(src_keys - dst_keys)

    # SET DIFFERENCE: dst_keys - src_keys gives us paths that are in the
    # destination but NOT in the source. These are "extra" files.
    r.missing_in_src = sorted(dst_keys - src_keys)

    # SET INTERSECTION: src_keys & dst_keys gives us paths that exist on BOTH
    # sides. These are the files we need to compare in detail.
    # We sort for deterministic output ordering.
    for rel in sorted(src_keys & dst_keys):
        # Look up the FileStat for this relative path in each scan.
        sf = src_stats[rel]  # Source FileStat
        df = dst_stats[rel]  # Destination FileStat

        # ── Check 1: Read errors ──────────────────────────────────────────────
        # If either side had a read error during hashing, we can't meaningfully
        # compare this file. Log the error and move on.
        # We check BEFORE the size comparison because a read error means the
        # checksum is None (not computed), so a checksum comparison would give
        # a false negative.
        if sf.checksum_error or df.checksum_error:
            r.checksum_errors.append({
                "path":      rel,
                "src_error": sf.checksum_error,  # None if source was fine
                "dst_error": df.checksum_error,  # None if destination was fine
            })
            continue  # Don't check this file further; go to the next file

        # ── Check 2: Size mismatch ────────────────────────────────────────────
        # Compare byte counts. If they differ, the files are definitely different.
        # We record the difference (diff_bytes) which can be negative (dst smaller
        # than src) or positive (dst larger than src, e.g. if the file grew).
        # We skip the checksum check if sizes differ — no point reading both files
        # fully when we already know they're different. Also, most hash algorithms
        # incorporate the data length, so a size mismatch would cause a hash
        # mismatch anyway.
        if sf.size != df.size:
            r.size_mismatches.append({
                "path":       rel,
                "src_size":   sf.size,
                "dst_size":   df.size,
                "diff_bytes": df.size - sf.size,  # Negative = dst is smaller (truncated)
            })
            continue  # Don't check checksums for this file; go to the next

        # ── Check 3: Checksum mismatch ────────────────────────────────────────
        # This is the critical check. Files have the same size but we hash them
        # to verify every single byte matches.
        #
        # The `sf.checksum and df.checksum` guard skips this check if either
        # checksum is None (which happens when --hash none was specified).
        # Without that guard, None != None would be False (they'd match), which
        # is actually fine, but the explicit guard makes the intent clear.
        if sf.checksum and df.checksum and sf.checksum != df.checksum:
            r.checksum_mismatches.append({
                "path":           rel,
                "src_checksum":   sf.checksum,  # Full hex digest from source
                "dst_checksum":   df.checksum,  # Full hex digest from destination
                "size":           sf.size,       # Size (same on both sides at this point)
            })
            continue  # This file failed; go to the next

        # ── All checks passed ─────────────────────────────────────────────────
        # If we reach here, the file exists on both sides, has the same size,
        # and (if hashing was enabled) has the same hash. It's verified good.
        r.ok_files += 1

    return r


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5: REPORTING
# ══════════════════════════════════════════════════════════════════════════════


def fmt_bytes(n: int) -> str:
    """
    Convert a raw byte count to a human-readable string with appropriate units.

    This function uses a simple iterative approach: keep dividing by 1024 and
    advancing to the next unit until the value is less than 1024.

    We use 1024-based units (kibibytes etc.) rather than 1000-based (kilobytes)
    because that's how filesystem sizes are traditionally reported. The labels
    use the informal "KB/MB/GB" names rather than the technically correct
    "KiB/MiB/GiB" to match what users are used to seeing.

    Args:
        n (int): File size in bytes. Can be 0 or any positive integer.

    Returns:
        str: Human-readable size string, e.g. "1.5 GB", "204.0 MB", "512.0 B"

    Examples:
        fmt_bytes(0)              -> "0.0 B"
        fmt_bytes(1023)           -> "1023.0 B"
        fmt_bytes(1024)           -> "1.0 KB"
        fmt_bytes(1_073_741_824)  -> "1.0 GB"
        fmt_bytes(1_500_000_000)  -> "1.4 GB"
    """
    # Iterate through units from smallest to largest.
    # The loop exits early via `return` when the value is small enough to display.
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if abs(n) < 1024:
            # Value fits in this unit — format with one decimal place.
            # abs() handles the (theoretically possible) negative case.
            return f"{n:.1f} {unit}"
        # Value is too large for this unit — divide by 1024 and try the next.
        n /= 1024
    # If we exhaust all units (> 1 PB), fall through to PB.
    return f"{n:.1f} PB"


def print_report(r: CompareResult, src: str, dst: str, algorithm: str):
    """
    Print a formatted human-readable comparison report to stdout.

    The report has three sections:
      1. HEADER: Shows the paths being compared, hash algorithm used, and timing.
      2. SUMMARY: One line per check type — either ✔ (all clear) or ✘ (issues found).
      3. DETAIL LISTS: For each problem category, lists up to 20 affected paths.
         If more than 20 exist, tells the user to check the report file.

    Args:
        r (CompareResult): The populated result object from compare().
        src (str): Absolute path of the source directory (for display).
        dst (str): Absolute path of the destination directory (for display).
        algorithm (str): Hash algorithm name (for display). "none" means no hashing.
    """
    W = 60  # Width of the horizontal divider lines (─ repeated W times)

    # ── Header ────────────────────────────────────────────────────────────────
    print()
    print(BOLD + "-" * W + RESET)
    print(BOLD + "  DRIVE COMPARISON REPORT" + RESET)
    print(f"  Source : {src}")
    print(f"  Dest   : {dst}")
    print(f"  Hash   : {algorithm.upper()}")
    print(f"  Time   : {r.elapsed_seconds:.1f}s")
    print("-" * W)
    # Show file counts and total data sizes for both sides.
    # The :>8, format spec right-aligns the number in an 8-character field
    # with comma-thousands separators, e.g. "  123,456".
    print(f"  Source files : {r.total_src_files:>8,}   ({fmt_bytes(r.src_total_bytes)})")
    print(f"  Dest   files : {r.total_dst_files:>8,}   ({fmt_bytes(r.dst_total_bytes)})")
    print("-" * W)

    # ── Summary section ───────────────────────────────────────────────────────
    # Count total issue categories (not total files — each category counts as 1
    # even if it has many affected files). This drives the final verdict message.
    issues = (len(r.missing_in_dest) + len(r.size_mismatches) +
              len(r.checksum_mismatches) + len(r.checksum_errors))

    # Print one line per check category. Each line is either green or red.
    if r.missing_in_dest:
        err(f"Missing in destination  : {len(r.missing_in_dest):,} file(s)")
    else:
        ok("No files missing in destination")

    if r.missing_in_src:
        # Use warn (yellow) rather than err (red) for extra files — they're
        # unexpected but don't mean the source data is lost or corrupted.
        warn(f"Extra files in destination (not in source): {len(r.missing_in_src):,}")
    else:
        ok("No extra files in destination")

    if r.size_mismatches:
        err(f"Size mismatches          : {len(r.size_mismatches):,} file(s)")
    else:
        ok("No size mismatches")

    # Only show the checksum result line if hashing was actually performed.
    # If --hash none was used, checksums are all None, so we'd get misleading
    # "All checksums match" messages when we never actually checked anything.
    if algorithm != "none":
        if r.checksum_mismatches:
            err(f"Checksum mismatches      : {len(r.checksum_mismatches):,} file(s)  <- likely corruption")
        else:
            ok("All checksums match — no corruption detected")

    if r.checksum_errors:
        # Use warn rather than err — these are unreadable files, which is
        # concerning, but not the same as a confirmed mismatch.
        warn(f"Checksum errors (unreadable): {len(r.checksum_errors):,} file(s)")

    # ── Verdict ───────────────────────────────────────────────────────────────
    print("-" * W)
    if issues == 0:
        print(GREEN + BOLD + "  ALL CLEAR — copy is verified complete and intact." + RESET)
    else:
        print(RED + BOLD + f"  {issues} issue category(ies) found — see details below or in the report file." + RESET)
    print("-" * W + "\n")

    # ── Detail lists ──────────────────────────────────────────────────────────
    # Show the first LIMIT items from each problem list. Showing all could
    # spam the terminal with thousands of lines for large problems.
    LIMIT = 20  # Maximum number of paths to print per category

    def show_list(title, items, key="path"):
        """
        Inner helper to print a titled list of problem items.

        Args:
            title (str): Section heading (e.g. "Missing in destination:")
            items (list): Either a list of strings or list of dicts.
            key (str): If items are dicts, which key to display as the path.
                       Defaults to "path" which all our dict items use.
        """
        if not items:
            return  # Nothing to show for this category — skip silently

        print(BOLD + f"  {title}" + RESET)
        # Print up to LIMIT items. For string items, print directly.
        # For dict items, extract the "path" (or specified key) field.
        for i in items[:LIMIT]:
            print(f"    {i if isinstance(i, str) else i[key]}")

        # If there are more items than we're showing, tell the user where to
        # find the full list.
        if len(items) > LIMIT:
            print(f"    … and {len(items) - LIMIT} more (see report file)")
        print()  # Blank line between categories

    # Call show_list for each problem category. Categories with no issues
    # are silently skipped by the `if not items: return` guard inside show_list.
    show_list("Missing in destination:", r.missing_in_dest)
    show_list("Extra in destination:", r.missing_in_src)
    show_list("Size mismatches:", r.size_mismatches)
    show_list("Checksum mismatches (CORRUPTED):", r.checksum_mismatches)
    show_list("Unreadable files:", r.checksum_errors)


def write_json_report(r: CompareResult, src: str, dst: str, algorithm: str, out_path: str):
    """
    Write a complete machine-readable JSON report to disk.

    The JSON report contains EVERYTHING — all problem lists in full (no LIMIT
    truncation like the terminal output), plus metadata about the run. It's
    useful for:
      - Scripting: pipe the output through `jq` to filter/query specific issues
      - Archiving: keep a record of when and how the verification was done
      - Programmatic follow-up: e.g. automatically re-copy the missing files

    The output uses indent=2 for human readability. For very large results
    (millions of files) this makes the file large; you could change indent=None
    for compact output if size matters.

    Args:
        r (CompareResult): The populated result object.
        src (str): Source directory absolute path.
        dst (str): Destination directory absolute path.
        algorithm (str): Hash algorithm name.
        out_path (str): Full path to write the .json file to.
    """
    # Build a plain Python dict. json.dump() can serialise dicts, lists, strings,
    # ints, floats, and None natively. Our CompareResult contains only these types.
    data = {
        # ISO 8601 timestamp: when this report was generated.
        # isoformat() produces e.g. "2024-01-15T14:32:07.123456"
        "generated":      datetime.now().isoformat(),
        "source":         src,
        "destination":    dst,
        "hash_algorithm": algorithm,
        "elapsed_seconds": r.elapsed_seconds,

        # A compact summary section — useful for at-a-glance status checks
        # without parsing the full problem lists.
        "summary": {
            "source_files":        r.total_src_files,
            "dest_files":          r.total_dst_files,
            "source_bytes":        r.src_total_bytes,
            "dest_bytes":          r.dst_total_bytes,
            "ok_files":            r.ok_files,
            "missing_in_dest":     len(r.missing_in_dest),
            "extra_in_dest":       len(r.missing_in_src),
            "size_mismatches":     len(r.size_mismatches),
            "checksum_mismatches": len(r.checksum_mismatches),
            "checksum_errors":     len(r.checksum_errors),
        },

        # Full problem lists — all items, not truncated.
        "missing_in_dest":     r.missing_in_dest,      # list of strings
        "extra_in_dest":       r.missing_in_src,       # list of strings
        "size_mismatches":     r.size_mismatches,      # list of dicts
        "checksum_mismatches": r.checksum_mismatches,  # list of dicts
        "checksum_errors":     r.checksum_errors,      # list of dicts
    }

    # Write the JSON file. We use UTF-8 encoding explicitly to handle filenames
    # that contain non-ASCII characters (accented chars, Chinese, emoji paths, etc.)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        # indent=2: pretty-print with 2-space indentation
        # ensure_ascii=False (default): allow non-ASCII characters in strings
        json.dump(data, f, indent=2)

    ok(f"JSON report saved -> {out_path}")


def write_csv_report(r: CompareResult, out_path: str):
    """
    Write a flat CSV report of all issues to disk.

    CSV (Comma-Separated Values) is the most universally compatible format for
    tabular data. It can be opened directly in Excel, Numbers, Google Sheets,
    LibreOffice Calc, or processed with standard Unix tools like `awk` and `grep`.

    The CSV has three columns:
        status  — one of: MISSING_IN_DEST, EXTRA_IN_DEST, SIZE_MISMATCH,
                          CHECKSUM_MISMATCH, READ_ERROR
        path    — the relative path of the affected file
        detail  — additional context (e.g. "src=1024 dst=1020" for size mismatches)

    All problem categories are flattened into a single list of rows. This makes
    it easy to sort/filter in a spreadsheet (e.g. "show me all SIZE_MISMATCH rows").

    Args:
        r (CompareResult): The populated result object.
        out_path (str): Full path to write the .csv file to.
    """
    # Build a flat list of dicts — one dict per problem file.
    rows = []

    # Missing files: just the path, no extra detail needed.
    for p in r.missing_in_dest:
        rows.append({"status": "MISSING_IN_DEST", "path": p, "detail": ""})

    # Extra files: just the path.
    for p in r.missing_in_src:
        rows.append({"status": "EXTRA_IN_DEST", "path": p, "detail": ""})

    # Size mismatches: include both sizes so the user can see how big the
    # discrepancy is. A 1-byte diff might be a line-ending issue;
    # a 500 MB diff is a clearly truncated copy.
    for m in r.size_mismatches:
        rows.append({
            "status": "SIZE_MISMATCH",
            "path":   m["path"],
            "detail": f"src={m['src_size']} dst={m['dst_size']}"
        })

    # Checksum mismatches: include truncated versions of both hashes. The full
    # hash is in the JSON report; here we just show enough to confirm they differ.
    # [:12] takes the first 12 hex characters (48 bits of the hash). The "…"
    # ellipsis signals that the hash is truncated, not that the file path is.
    for m in r.checksum_mismatches:
        rows.append({
            "status": "CHECKSUM_MISMATCH",
            "path":   m["path"],
            "detail": f"src={m['src_checksum'][:12]}... dst={m['dst_checksum'][:12]}..."
        })

    # Read errors: show the error message so the user knows WHY it failed
    # (e.g. "Permission denied" vs "No such file or directory").
    # We use `or` to pick whichever error field is non-None (at least one will be).
    for m in r.checksum_errors:
        rows.append({
            "status": "READ_ERROR",
            "path":   m["path"],
            "detail": str(m.get("src_error") or m.get("dst_error"))
        })

    # Write the CSV file.
    # newline="" is required by the csv module on Windows to prevent it from
    # adding an extra \r before each \n (the csv writer handles line endings itself).
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        # DictWriter writes rows from dictionaries. fieldnames specifies column
        # order — the dict keys can be in any order, but the CSV columns will
        # always appear in this sequence.
        writer = csv.DictWriter(f, fieldnames=["status", "path", "detail"])
        # writeheader() writes the column names as the first row of the file.
        writer.writeheader()
        # writerows() writes all the dicts in order. Each dict becomes one row.
        writer.writerows(rows)

    ok(f"CSV  report saved -> {out_path}")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6: ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════


def main():
    """
    Script entry point — parse arguments, orchestrate the full comparison pipeline.

    This function is the conductor: it sets up argument parsing, validates inputs,
    calls scan_tree() twice (once for each drive), calls compare(), then calls
    the reporting functions. It measures wall-clock time for the whole operation.

    Exit codes:
        0 — No issues found. Safe to treat the copy as verified.
        1 — One or more issues detected. The user should investigate.

    These exit codes make the script useful in shell pipelines and CI scripts:
        python compare_drives.py /src /dst && echo "VERIFIED OK" || echo "ISSUES FOUND"
    """
    # ── Optional GUI folder picker ─────────────────────────────────────────────
    # Uses tkinter's built-in directory dialog — no third-party packages needed.
    # The hidden root window is created just to host the dialogs, then destroyed.
    # If the user cancels either dialog, the script exits cleanly.
    import tkinter as tk
    from tkinter import filedialog, messagebox

    def pick_folders_with_gui():
        """
        Opens two folder-picker dialogs in sequence and returns the chosen paths.
        Call this INSTEAD of passing source/destination as command-line arguments.

        Returns:
            tuple[str, str]: (source_path, destination_path)
            Calls sys.exit() if the user cancels either dialog.
        """
        # Create a hidden root window. tkinter requires one to exist before
        # any dialogs can open, but we don't actually want a window visible.
        root = tk.Tk()
        root.withdraw()          # Hide the root window — we only want the dialogs
        root.attributes("-topmost", True)  # Ensure dialogs appear on top of other windows

        # ── Dialog 1: Pick source ──────────────────────────────────────────────
        messagebox.showinfo(
            "Step 1 of 2",
            "Select the SOURCE folder (the original drive you copied FROM).",
            parent=root
        )
        source = filedialog.askdirectory(
            title="Select SOURCE folder (original drive)",
            parent=root
        )
        if not source:
            # User hit Cancel — exit gracefully rather than crashing
            sys.exit("No source folder selected. Exiting.")

        # ── Dialog 2: Pick destination ─────────────────────────────────────────
        messagebox.showinfo(
            "Step 2 of 2",
            "Now select the DESTINATION folder (the copy you want to verify).",
            parent=root
        )
        destination = filedialog.askdirectory(
            title="Select DESTINATION folder (the copy to verify)",
            parent=root
        )
        if not destination:
            sys.exit("No destination folder selected. Exiting.")

        root.destroy()  # Clean up the hidden root window now that we're done with it
        return source, destination


    # ── Argument parser setup ─────────────────────────────────────────────────
    # argparse.ArgumentParser handles all command-line argument parsing.
    # description: shown in --help above the argument list.
    # formatter_class=RawDescriptionHelpFormatter: preserves newlines/indentation
    #   in the `epilog` string (which contains our usage examples).
    parser = argparse.ArgumentParser(
        description="Deep folder comparison — verify a drive copy is complete and uncorrupted.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic size+checksum check (MD5, fast enough for most drives)
  python compare_drives.py /Volumes/SourceDrive /Volumes/BackupDrive

  # Use SHA256 for higher integrity assurance
  python compare_drives.py /source /dest --hash sha256

  # Skip checksums — size-only, very fast
  python compare_drives.py /source /dest --hash none

  # Save results to custom path
  python compare_drives.py /source /dest --out ~/Desktop/report

  # Use 8 threads (good for fast NVMe drives)
  python compare_drives.py /source /dest --workers 8
""")

    # Positional arguments (required, no -- prefix needed):
    # Change nargs="?" makes the arguments optional (defaults to None if omitted)
    parser.add_argument("source",      nargs="?", default=None, help="Source folder / drive root")
    parser.add_argument("destination", nargs="?", default=None, help="Destination folder / drive root")

    # Optional arguments (keyword, with -- prefix):
    parser.add_argument("--hash", default="md5",
        choices=["md5", "sha1", "sha256", "sha512", "none"],
        help="Hash algorithm for content integrity check (default: md5). "
             "Use 'none' to skip hashing and only check file existence and size.")

    parser.add_argument("--workers", type=int, default=4,
        help="Number of parallel threads for hashing (default: 4). "
             "Increase for SSDs (8-16), decrease for spinning HDDs (1-2).")

    parser.add_argument("--out", default="",
        help="Base path for output report files, without extension "
             "(default: ./drive_compare_YYYYMMDD_HHMMSS). "
             "The script appends .json and .csv automatically.")

    parser.add_argument("--no-json", action="store_true",
        help="Skip writing the JSON report file.")

    parser.add_argument("--no-csv", action="store_true",
        help="Skip writing the CSV report file.")

    parser.add_argument("--verbose", action="store_true",
        help="Print a progress update after every single file instead of every 5%%.")

    # parse_args() reads sys.argv[1:], matches arguments against the definitions
    # above, and returns a Namespace object where each argument is an attribute.
    # If required arguments are missing or types are wrong, argparse prints an
    # error message and exits with code 2 automatically.
    args = parser.parse_args()

    # ── GUI folder picker ──────────────────────────────────────────────────────
    # If no paths were passed on the command line, open folder picker dialogs.
    # This lets you double-click the script and pick folders visually instead
    # of typing paths. Command-line arguments still work as normal if provided.
    if not args.source or not args.destination:
        args.source, args.destination = pick_folders_with_gui()

    # ── Input validation ──────────────────────────────────────────────────────
    # Convert the user-supplied path strings to Path objects for easier handling.
    src_root = Path(args.source)
    dst_root = Path(args.destination)

    # Validate that both paths exist and are directories.
    # We check both before starting the scan so we don't waste time scanning
    # the source only to then discover the destination doesn't exist.
    for p, label in [(src_root, "Source"), (dst_root, "Destination")]:
        if not p.exists():
            # sys.exit() with a string prints the string to stderr and exits
            # with code 1. This is the standard Python way to abort with an error.
            sys.exit(f"ERROR: {label} path does not exist: {p}")
        if not p.is_dir():
            sys.exit(f"ERROR: {label} path is not a directory: {p}")

    # ── Output path setup ─────────────────────────────────────────────────────
    # Create a timestamp string for the default output filename.
    # strftime format: %Y=4-digit year, %m=month, %d=day, %H=hour, %M=min, %S=sec
    # Example: "20240115_143207"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Use the user's --out value if provided, otherwise generate a timestamped name.
    # The `or` operator returns the right side if the left side is falsy (empty string).
    out_base = args.out or f"drive_compare_{timestamp}"
    # Report files will be: out_base + ".json" and out_base + ".csv"

    # ── Startup banner ────────────────────────────────────────────────────────
    print()
    print(BOLD + CYAN + "  DRIVE COMPARISON TOOL" + RESET)
    print(f"  Source      : {src_root.resolve()}")   # resolve() shows the absolute path
    print(f"  Destination : {dst_root.resolve()}")
    print(f"  Hash        : {args.hash.upper()}")
    print(f"  Workers     : {args.workers}")
    print()

    # ── Main pipeline ─────────────────────────────────────────────────────────
    # Record the start time so we can compute elapsed_seconds at the end.
    start = datetime.now()

    # Phase 1: Scan the source tree.
    # This walks all files and (unless --hash none) reads and hashes each one.
    info("Scanning source …")
    src_stats = scan_tree(src_root, args.hash, args.workers, args.verbose)
    ok(f"Source scanned: {len(src_stats):,} files")

    # Phase 2: Scan the destination tree.
    # Same process as above but for the destination drive.
    info("Scanning destination …")
    dst_stats = scan_tree(dst_root, args.hash, args.workers, args.verbose)
    ok(f"Destination scanned: {len(dst_stats):,} files")

    # Phase 3: Compare the two scan results.
    # This is CPU-bound set operations and dict lookups — very fast.
    info("Comparing …")
    result = compare(src_stats, dst_stats)

    # Record how long the whole thing took (scan + hash + compare).
    result.elapsed_seconds = (datetime.now() - start).total_seconds()

    # Phase 4: Output results.
    # Terminal report first — the user sees immediate feedback.
    print_report(result, str(src_root.resolve()), str(dst_root.resolve()), args.hash)

    # Write file reports (unless disabled by --no-json / --no-csv).
    if not args.no_json:
        write_json_report(result, str(src_root.resolve()), str(dst_root.resolve()),
                          args.hash, f"{out_base}.json")
    if not args.no_csv:
        write_csv_report(result, f"{out_base}.csv")

    # ── Exit code ─────────────────────────────────────────────────────────────
    # Count total issues across all categories. Note: we count the number of
    # CATEGORIES with issues (each non-empty list is +1), NOT the total number
    # of affected files. This is intentional — the exit code is binary (ok/not ok)
    # so the exact count doesn't matter, only whether anything failed.
    issues = (len(result.missing_in_dest) + len(result.size_mismatches) +
              len(result.checksum_mismatches) + len(result.checksum_errors))

    # Exit 0 = success (no issues). Exit 1 = failure (at least one issue).
    # Shell scripts can check: `if python compare_drives.py ...; then ...; fi`
    sys.exit(0 if issues == 0 else 1)


# ── Module guard ──────────────────────────────────────────────────────────────
# This block only executes when the script is run directly (e.g. `python compare_drives.py`)
# and NOT when it's imported as a module by another Python script.
# When Python imports this file as a module, __name__ is set to the module's name
# ("compare_drives"), not "__main__". The guard prevents main() from running
# automatically on import, which would be confusing and break anyone trying to
# import individual functions (e.g. compute_checksum or scan_tree) from this file.
if __name__ == "__main__":
    main()
