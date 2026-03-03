"""
tools/cache.py - Disk-Based Result Caching for ProteusIQ

Provides a SequenceCache that stores analysis results on disk using
diskcache. This prevents redundant API calls when the same sequence
is analyzed multiple times (by the same or different users).

Cache key = SHA-256(sequence + options hash)
Cache TTL = 7 days by default

Large binary data (PDB file content) is excluded from the cache to
save disk space — it is re-downloaded on cache hit if needed.
"""

import os
import json
import hashlib
import logging

logger = logging.getLogger(__name__)

# ─── Cache configuration ───
CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), ".proteusiq_cache")
CACHE_TTL = 7 * 24 * 3600  # 7 days in seconds

# Keys to exclude from cached data (too large or ephemeral)
EXCLUDE_KEYS = {"pdb_content", "pdb_path"}

# Try to import diskcache; fall back to no-op if unavailable
try:
    import diskcache
    HAS_DISKCACHE = True
except ImportError:
    HAS_DISKCACHE = False
    logger.info("diskcache not installed — caching disabled. Install with: pip install diskcache")


def _make_cache_key(
    sequence: str,
    skip_conservation: bool = False,
    compute_msa: bool = False,
    build_tree: bool = False,
    run_interpro: bool = False,
) -> str:
    """
    Generate a unique cache key from the sequence and analysis options.

    Returns:
        SHA-256 hex digest string.
    """
    key_parts = {
        "seq": sequence,
        "skip_cons": skip_conservation,
        "msa": compute_msa,
        "tree": build_tree,
        "interpro": run_interpro,
    }
    key_string = json.dumps(key_parts, sort_keys=True)
    return hashlib.sha256(key_string.encode("utf-8")).hexdigest()


def _strip_large_fields(data: dict) -> dict:
    """Remove large/ephemeral fields before storing in cache."""
    cleaned = {}
    for k, v in data.items():
        if k in EXCLUDE_KEYS:
            continue
        cleaned[k] = v
    return cleaned


class SequenceCache:
    """
    Disk-based cache for protein analysis results.

    Usage:
        cache = SequenceCache()
        key = cache.make_key(sequence, ...)
        if cache.has(key):
            return cache.get(key)
        # ... run analysis ...
        cache.set(key, results)
    """

    def __init__(self, cache_dir: str = CACHE_DIR, ttl: int = CACHE_TTL):
        self.enabled = HAS_DISKCACHE
        self.ttl = ttl

        if self.enabled:
            try:
                self._cache = diskcache.Cache(cache_dir, size_limit=500 * 1024 * 1024)  # 500MB max
                logger.info("Disk cache initialized at %s", cache_dir)
            except Exception as e:
                logger.warning("Failed to initialize disk cache: %s", e)
                self.enabled = False
                self._cache = None
        else:
            self._cache = None

    def make_key(self, sequence: str, **options) -> str:
        """Generate cache key from sequence and options."""
        return _make_cache_key(sequence, **options)

    def has(self, key: str) -> bool:
        """Check if a key exists in the cache."""
        if not self.enabled:
            return False
        try:
            return key in self._cache
        except Exception:
            return False

    def get(self, key: str) -> dict:
        """
        Retrieve cached results.

        Returns:
            Cached data dict, or None if not found.
        """
        if not self.enabled:
            return None
        try:
            result = self._cache.get(key)
            if result is not None:
                logger.info("Cache HIT for key %s...", key[:16])
            return result
        except Exception as e:
            logger.warning("Cache read error: %s", e)
            return None

    def set(self, key: str, data: dict) -> None:
        """
        Store analysis results in the cache.

        Large fields (pdb_content, pdb_path) are stripped automatically.
        """
        if not self.enabled:
            return
        try:
            cleaned = _strip_large_fields(data)
            self._cache.set(key, cleaned, expire=self.ttl)
            logger.info("Cache SET for key %s... (TTL=%ds)", key[:16], self.ttl)
        except Exception as e:
            logger.warning("Cache write error: %s", e)

    def clear(self) -> None:
        """Clear the entire cache."""
        if not self.enabled:
            return
        try:
            self._cache.clear()
            logger.info("Cache cleared")
        except Exception as e:
            logger.warning("Cache clear error: %s", e)

    def close(self) -> None:
        """Close the cache (release file handles)."""
        if self._cache is not None:
            try:
                self._cache.close()
            except Exception:
                pass
