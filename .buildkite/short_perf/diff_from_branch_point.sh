#!/bin/bash

# Based on example from monorepo-diff-buildkite-plugin
set -ue

BRANCH_POINT_COMMIT=$(diff -u <(git rev-list --first-parent HEAD) <(git rev-list --first-parent main) | sed -ne 's/^ //p' | head -1)

diff=$(git diff --name-only "$BRANCH_POINT_COMMIT")
echo "$diff"
