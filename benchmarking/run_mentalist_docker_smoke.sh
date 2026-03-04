#!/usr/bin/env bash
set -euo pipefail

IMAGE="${IMAGE:-matnguyen/mentalist:latest}"

docker pull "${IMAGE}"
docker run --rm "${IMAGE}" mentalist -h
