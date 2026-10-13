# Single-container deployment for the Quantum Drug Discovery Explorer.
#
# One process (`uvicorn app.main:app`) serves the JSON API *and* the frontend in
# `frontend/static/` (see `app/main.py`, and docs/DEPLOYMENT.md for why this is
# the recommended shape instead of a serverless platform).
#
# PySCF and qiskit-nature are deliberately NOT installed here. They are used
# offline only, to generate `app/data/hamiltonians/*.json`; the running API reads
# those files and never calls a quantum-chemistry package live. See
# docs/DEPLOYMENT.md.
FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PORT=8000

WORKDIR /srv/app

# Dependencies first, so editing application code does not re-install them.
COPY requirements.txt ./
RUN python -m pip install --upgrade pip \
    && python -m pip install --no-cache-dir -r requirements.txt

# Application code, the precomputed Hamiltonians, the curated manifests and the
# frontend this app serves itself.
COPY app ./app
COPY frontend ./frontend
COPY README.md SPEC.md LIMITATIONS.md RESULTS.md ./

# Run as a non-root user. Everything the app writes at runtime (the disk cache,
# the SQLite job store) lives under /srv/app/app/data, so that tree must be
# writable by this user and, on hosts with ephemeral disks, is disposable.
RUN useradd --create-home --uid 1000 appuser \
    && chown -R appuser:appuser /srv/app
USER appuser

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD python -c "import os, sys, urllib.request; port = os.environ.get('PORT', '8000'); sys.exit(0 if urllib.request.urlopen('http://127.0.0.1:' + port + '/health', timeout=4).status == 200 else 1)"

# Bind 0.0.0.0: the default HOST of 127.0.0.1 would be unreachable from outside
# the container. PORT is read from the environment because hosts such as Render
# assign it (Hugging Face Spaces expects 7860).
CMD ["sh", "-c", "uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}"]
