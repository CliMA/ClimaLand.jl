#!/bin/bash
# Sync ClimaLand.jl to tofu including git repository

SOURCE="/glade/derecho/scratch/reich/ClimaLand.jl/"
TARGET="egreich@tofu.gps.caltech.edu:/projects/egreich/ClimaLand.jl/"

echo "=========================================="
echo "Syncing ClimaLand.jl to tofu"
echo "=========================================="
echo "Source: $SOURCE"
echo "Target: $TARGET"
echo ""
echo "This will copy:"
echo "  - All source code"
echo "  - Git repository (.git/)"
echo "  - Calibration experiments and results"
echo "  - Compiled modules"
echo ""
read -p "Continue? (y/n) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 1
fi

# Sync everything including .git
rsync -avz --progress \
    --exclude='*.swp' \
    --exclude='*.swo' \
    --exclude='.DS_Store' \
    "$SOURCE" "$TARGET"

echo ""
echo "=========================================="
echo "✓ Sync complete!"
echo "=========================================="
echo ""
echo "To verify on tofu:"
echo "  ssh tofu.gps.caltech.edu"
echo "  cd /home/egreich/Projects/Climaexplore/ClimaLand.jl"
echo "  git status"
echo "  git log --oneline -5"
