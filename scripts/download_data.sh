#!/usr/bin/env bash
# Restore the data folders left out of this repo. Run from the repo root.
#
#   scripts/download_data.sh                    prediction data (about 1.6 GB)
#   scripts/download_data.sh --with-swissprot   also SwissProt (about 0.5 GB more)
#   scripts/download_data.sh --from /path/to/DeepPBI-KG   copy from a local checkout
#   scripts/download_data.sh --only all_phage_db          one folder, for testing
set -euo pipefail

UPSTREAM=https://github.com/Tongqing-Wei/DeepPBI-KG.git
COMMIT=956850cac597255b1aff869786f34d144cff583b
folders=(phage_raw_data host_raw_data all_phage_db all_host_db)
source_dir=""

while [ $# -gt 0 ]; do
  case "$1" in
    --with-swissprot) folders+=(code/swiss); shift ;;
    --from) source_dir="$2"; shift 2 ;;
    --only) folders=("$2"); shift 2 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

[ -f DeepPBI-KG.py ] || { echo "run this from the repo root" >&2; exit 2; }

expected_count() {
  case "$1" in
    phage_raw_data) echo 3513 ;;
    host_raw_data) echo 395 ;;
    all_phage_db|all_host_db) echo 8 ;;
    code/swiss) echo 2 ;;
  esac
}

todo=()
for folder in "${folders[@]}"; do
  if [ -e "$folder" ]; then echo "keeping existing $folder"; else todo+=("$folder"); fi
done

if [ ${#todo[@]} -gt 0 ]; then
  if [ -n "$source_dir" ]; then
    for folder in "${todo[@]}"; do
      mkdir -p "$(dirname "$folder")"
      cp -a "$source_dir/$folder" "$folder"
    done
  else
    command -v git-lfs >/dev/null || { echo "git-lfs is required" >&2; exit 1; }
    work=$(mktemp -d .download_data.XXXXXX)
    trap 'rm -rf "$work"' EXIT
    git -C "$work" init -q
    git -C "$work" remote add origin "$UPSTREAM"
    git -C "$work" sparse-checkout set "${todo[@]}"
    git -C "$work" fetch -q --depth 1 --filter=blob:none origin "$COMMIT"
    git -C "$work" checkout -q FETCH_HEAD
    for folder in "${todo[@]}"; do
      mkdir -p "$(dirname "$folder")"
      mv "$work/$folder" "$folder"
    done
  fi
fi

status=0
for folder in "${folders[@]}"; do
  found=$(find "$folder" -type f | wc -l)
  want=$(expected_count "$folder")
  if [ "$found" -ne "$want" ]; then
    echo "FAIL $folder: $found files, expected $want"; status=1
  elif grep -rlq --include='*.nsq' --include='swissprot*' '^version https://git-lfs' "$folder"; then
    echo "FAIL $folder: holds Git LFS pointers, not data"; status=1
  else
    echo "ok   $folder ($found files)"
  fi
done
exit $status
