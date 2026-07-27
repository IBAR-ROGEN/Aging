#!/usr/bin/env bash
exec "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/dev/ukbb_ci_compliance_audit.sh" "$@"
