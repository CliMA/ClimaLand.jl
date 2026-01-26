#!/bin/bash
# Rebase strategy to keep uSPAC, SMAP, and aridity changes
# Run this script to help preserve your work during rebase

set -e

echo "=========================================="
echo "Git Rebase Strategy Helper"
echo "=========================================="
echo ""

# Define files to keep (your versions)
KEEP_FILES=(
    # Entire calibration folder
    "experiments/calibration/"
    
    # Vegetation standalone models
    "src/standalone/Vegetation/Canopy.jl"
    "src/standalone/Vegetation/stomatalconductance.jl"
    "src/standalone/Vegetation/uspacmodel.jl"
    "src/standalone/Vegetation/photosynthesis.jl"
)

# Step 1: Create a backup branch
echo "Step 1: Creating backup branch..."
BACKUP_BRANCH="backup-before-rebase-$(date +%Y%m%d-%H%M%S)"
git branch "$BACKUP_BRANCH"
echo "✓ Created backup branch: $BACKUP_BRANCH"
echo ""

# Step 2: Stash any uncommitted changes
echo "Step 2: Checking for uncommitted changes..."
if ! git diff-index --quiet HEAD --; then
    echo "Found uncommitted changes. Stashing..."
    git stash push -m "Pre-rebase stash $(date +%Y%m%d-%H%M%S)"
    STASHED=1
else
    echo "✓ No uncommitted changes"
    STASHED=0
fi
echo ""

# Step 3: Show what files you've modified
echo "Step 3: Files modified in your branch vs main:"
echo "----------------------------------------------"
git diff --name-status origin/main...HEAD | head -50
echo ""

# Step 4: Instructions for interactive rebase
echo "=========================================="
echo "REBASE INSTRUCTIONS"
echo "=========================================="
echo ""
echo "To rebase while keeping your changes:"
echo ""
echo "1. Start interactive rebase:"
echo "   git rebase -i origin/main"
echo ""
echo "2. Git will open an editor with your commits. Keep all your commits as 'pick'"
echo ""
echo "3. When conflicts occur, for files you want to KEEP (your version):"
echo "   git checkout --ours <file>"
echo ""
echo "4. For files you want to ACCEPT from main:"
echo "   git checkout --theirs <file>"
echo ""
echo "5. Files to keep YOUR version (use --ours):"
for file in "${KEEP_FILES[@]}"; do
    echo "   - $file"
done
echo ""
echo "6. After resolving each conflict:"
echo "   git add <resolved-files>"
echo "   git rebase --continue"
echo ""
echo "7. If something goes wrong:"
echo "   git rebase --abort"
echo "   git checkout $BACKUP_BRANCH"
echo ""
echo "=========================================="
echo ""

# Step 5: Create helper script for conflict resolution
cat > /tmp/resolve_uspac_conflicts.sh << 'EOF'
#!/bin/bash
# Helper script to resolve conflicts by keeping your uSPAC/SMAP versions

echo "Resolving conflicts - keeping your versions of uSPAC/SMAP files..."

# Keep entire calibration folder
if git status --porcelain | grep -q "experiments/calibration/"; then
    git checkout --ours experiments/calibration/
    git add experiments/calibration/
    echo "✓ Kept your experiments/calibration/"
fi

# Keep vegetation files
for file in \
    "src/standalone/Vegetation/Canopy.jl" \
    "src/standalone/Vegetation/stomatalconductance.jl" \
    "src/standalone/Vegetation/uspacmodel.jl" \
    "src/standalone/Vegetation/photosynthesis.jl"
do
    if git status --porcelain | grep -q "$file"; then
        git checkout --ours "$file"
        git add "$file"
        echo "✓ Kept your $file"
    fi
done

echo ""
echo "Files resolved. Check 'git status' for any remaining conflicts."
echo "When ready, run: git rebase --continue"
EOF

chmod +x /tmp/resolve_uspac_conflicts.sh

echo "Created helper script: /tmp/resolve_uspac_conflicts.sh"
echo "During rebase, if you have conflicts in your uSPAC files, run:"
echo "  bash /tmp/resolve_uspac_conflicts.sh"
echo ""
echo "=========================================="
echo "Ready to rebase? Run: git rebase -i origin/main"
echo "=========================================="
