"""
tools/api_utils.py - Resilient HTTP Client & Rate Limiting

Provides a shared resilient_request() function that wraps requests with:
  - Exponential backoff with jitter (3 retries)
  - HTTP 429 (Too Many Requests) handling with Retry-After support
  - Global concurrency semaphore (max 3 simultaneous outbound API calls)
  - Configurable timeouts
  - NCBI API key support via environment variable

All API-calling modules should use resilient_request() instead of raw
requests.get/post to prevent IP blocking under concurrent usage.
"""

import os
import time
import random
import logging
import threading
import requests

logger = logging.getLogger(__name__)

# ─── Global concurrency limiter ───
# Prevents more than MAX_CONCURRENT outbound API calls at once,
# which keeps us within NCBI/EBI fair-use policies.
MAX_CONCURRENT = 3
_semaphore = threading.Semaphore(MAX_CONCURRENT)

# ─── Retry configuration ───
MAX_RETRIES = 3
BASE_DELAY = 1.0        # seconds
MAX_DELAY = 30.0         # cap on backoff
JITTER_RANGE = 0.5       # ±0.5s random jitter

# ─── Default timeouts ───
DEFAULT_CONNECT_TIMEOUT = 10   # seconds
DEFAULT_READ_TIMEOUT = 30      # seconds


def get_ncbi_params() -> dict:
    """
    Return NCBI-required parameters for API calls.

    If NCBI_API_KEY is set in the environment, includes it for higher
    rate limits (10 req/s instead of 3 req/s).

    Returns:
        Dict with 'tool' and 'email', plus 'api_key' if available.
    """
    params = {
        "tool": "ProteusIQ",
        "email": "proteusiq_tool@example.com",
    }
    api_key = os.environ.get("NCBI_API_KEY", "")
    if api_key:
        params["api_key"] = api_key
        logger.debug("Using NCBI API key for elevated rate limits")
    return params


def resilient_request(
    method: str,
    url: str,
    max_retries: int = MAX_RETRIES,
    base_delay: float = BASE_DELAY,
    timeout: tuple = None,
    **kwargs,
) -> requests.Response:
    """
    Make an HTTP request with automatic retry, backoff, and rate limiting.

    This function acquires a global semaphore before making the request,
    ensuring at most MAX_CONCURRENT outbound calls are active at any time.

    Args:
        method: HTTP method ('get' or 'post').
        url: Request URL.
        max_retries: Maximum number of retry attempts.
        base_delay: Base delay for exponential backoff (seconds).
        timeout: Optional (connect, read) timeout tuple.
        **kwargs: Additional arguments passed to requests.request().

    Returns:
        requests.Response object.

    Raises:
        requests.RequestException: If all retries are exhausted.
    """
    if timeout is None:
        timeout = (DEFAULT_CONNECT_TIMEOUT, DEFAULT_READ_TIMEOUT)
    kwargs.setdefault("timeout", timeout)

    last_exception = None

    for attempt in range(1, max_retries + 1):
        # Acquire semaphore to limit concurrency
        _semaphore.acquire()
        try:
            response = requests.request(method, url, **kwargs)

            # Handle HTTP 429 Too Many Requests
            if response.status_code == 429:
                retry_after = response.headers.get("Retry-After")
                if retry_after:
                    try:
                        wait_time = int(retry_after)
                    except ValueError:
                        wait_time = base_delay * (2 ** (attempt - 1))
                else:
                    wait_time = base_delay * (2 ** (attempt - 1))

                wait_time = min(wait_time, MAX_DELAY)
                logger.warning(
                    "HTTP 429 from %s — retrying in %.1fs (attempt %d/%d)",
                    url[:80], wait_time, attempt, max_retries,
                )
                time.sleep(wait_time)
                continue

            # Handle server errors (500, 502, 503, 504)
            if response.status_code >= 500:
                delay = min(base_delay * (2 ** (attempt - 1)), MAX_DELAY)
                jitter = random.uniform(-JITTER_RANGE, JITTER_RANGE)
                wait_time = max(0.1, delay + jitter)

                logger.warning(
                    "HTTP %d from %s — retrying in %.1fs (attempt %d/%d)",
                    response.status_code, url[:80], wait_time, attempt, max_retries,
                )
                time.sleep(wait_time)
                continue

            # Success or client error (4xx other than 429) — return immediately
            return response

        except requests.exceptions.ConnectionError as e:
            last_exception = e
            delay = min(base_delay * (2 ** (attempt - 1)), MAX_DELAY)
            jitter = random.uniform(-JITTER_RANGE, JITTER_RANGE)
            wait_time = max(0.1, delay + jitter)

            logger.warning(
                "Connection error to %s — retrying in %.1fs (attempt %d/%d): %s",
                url[:80], wait_time, attempt, max_retries, str(e)[:100],
            )
            time.sleep(wait_time)

        except requests.exceptions.Timeout as e:
            last_exception = e
            delay = min(base_delay * (2 ** (attempt - 1)), MAX_DELAY)
            logger.warning(
                "Timeout to %s — retrying in %.1fs (attempt %d/%d)",
                url[:80], delay, attempt, max_retries,
            )
            time.sleep(delay)

        except requests.RequestException as e:
            last_exception = e
            logger.warning(
                "Request error to %s (attempt %d/%d): %s",
                url[:80], attempt, max_retries, str(e)[:100],
            )
            if attempt < max_retries:
                time.sleep(base_delay)

        finally:
            _semaphore.release()

    # All retries exhausted
    if last_exception:
        raise last_exception
    raise requests.RequestException(
        f"All {max_retries} attempts to {url[:80]} failed"
    )
