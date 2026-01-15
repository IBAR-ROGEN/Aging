#!/bin/bash
# Script to find R installation on macOS

echo "Searching for R installation..."
echo ""

# Check common locations
locations=(
  "/Library/Frameworks/R.framework/Resources/bin/R"
  "/usr/local/bin/R"
  "/opt/homebrew/bin/R"
  "/usr/bin/R"
  "$HOME/.local/bin/R"
)

found=false
for loc in "${locations[@]}"; do
  if [ -f "$loc" ] && [ -x "$loc" ]; then
    echo "✓ Found R at: $loc"
    echo "  Version: $($loc --version 2>&1 | head -1)"
    echo ""
    echo "Add this to your .vscode/settings.json:"
    echo "  \"r.rpath.mac\": \"$loc\""
    found=true
    break
  fi
done

if [ "$found" = false ]; then
  echo "R not found in common locations."
  echo ""
  echo "To install R on macOS:"
  echo "1. Download from CRAN: https://cran.r-project.org/bin/macosx/"
  echo "2. Or install via Homebrew: brew install r"
  echo ""
  echo "After installation, run this script again to find the path."
fi
