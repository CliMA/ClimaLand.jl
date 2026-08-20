#!/bin/bash
# Start a Julia Distributed worker bound to the node's HSN interface.
#
# Why: `julia --worker` binds and advertises `getipaddr()`, which on Derecho
# compute nodes is the management network (10.169.0.0/16). Login nodes
# cannot reach that network, so a driver running in tmux on a login node
# never connects and every worker exits with "Master process (id 1) could
# not connect within 60.0 seconds". The HSN (10.14.0.0/18) is shared by
# login and compute nodes, so binding there works for drivers on either.
set -u
BIND=$(ip -o -4 addr show | awk '$2 ~ /^hsn/ {print $4; exit}' | cut -d/ -f1)
if [ -z "$BIND" ]; then
    echo "julia_worker_hsn_bind.sh: no hsn interface found, binding default" >&2
    exec julia "$@"
fi
exec julia --bind-to "$BIND" "$@"
