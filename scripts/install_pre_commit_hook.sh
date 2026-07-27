#!/usr/bin/env bash
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/dev/install_pre_commit_hook.sh" "$@"
