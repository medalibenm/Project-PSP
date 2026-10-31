#!/usr/bin/env bash
# install_cron.sh - Installe la tache cron du projet, sans dupliquer une ancienne entree.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
command -v crontab >/dev/null || { echo "crontab introuvable" >&2; exit 5; }

{
    crontab -l 2>/dev/null | grep -v '# projet_maladie$' || true
    grep -v '^#' "$ROOT/cron/pipeline.cron" | sed "s|@ROOT@|$ROOT|g"
} | crontab -

echo "Tache installee :"
crontab -l | grep '# projet_maladie$'
