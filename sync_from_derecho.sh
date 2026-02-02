#!/bin/bash
# Run this script ON TOFU to pull changes from derecho

SOURCE="reich@derecho.hpc.ucar.edu:/glade/derecho/scratch/reich/ClimaLand.jl/"
TARGET="/projects/egreich/ClimaLand.jl/"

echo "=========================================="
echo "Syncing ClimaLand.jl from derecho to tofu"
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

# Create target directory if it doesn't exist
mkdir -p "$TARGET"

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
echo "Next steps:"
echo "  cd $TARGET"
echo "  git status"
echo "  git push origin er/uspac"
