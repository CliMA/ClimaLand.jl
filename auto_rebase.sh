#!/bin/bash
# Rebase er/uspac onto main while keeping ALL local changes

set -e

echo "================================"
echo "Rebasing er/uspac onto main"
echo "Keeping ALL local changes (ours)"
echo "================================"

# Store current branch
CURRENT_BRANCH=$(git branch --show-current)
echo "Current branch: $CURRENT_BRANCH"

if [ "$CURRENT_BRANCH" != "er/uspac" ]; then
    echo "ERROR: Not on er/uspac branch. Please checkout er/uspac first."
    exit 1
fi

# Fetch latest from remote
echo ""
echo "Fetching latest changes from remote..."
git fetch origin

# Get main branch
echo "Updating main branch..."
git fetch origin main:main

# Start rebase
echo ""
echo "Starting rebase onto main..."
echo "Strategy: Keep ALL changes from er/uspac (ours)"
echo ""

# Use strategy option to prefer our changes
git rebase -X ours origin/main

if [ $? -eq 0 ]; then
    echo ""
    echo "================================"
    echo "✓ Rebase completed successfully!"
    echo "================================"
    echo ""
    echo "Your er/uspac changes have been preserved."
    echo ""
    echo "To push to remote (force push required):"
    echo "  git push origin er/uspac --force-with-lease"
else
    echo ""
    echo "================================"
    echo "Rebase encountered conflicts"
    echo "================================"
    echo ""
    echo "For remaining conflicts, use:"
    echo "  git checkout --ours <file>    # Keep your version"
    echo "  git add <file>"
    echo "  git rebase --continue"
    echo ""
    echo "To abort the rebase:"
    echo "  git rebase --abort"
fi
EOF
chmod +x rebase_keep_ours.sh
./rebase_keep_ours.sh