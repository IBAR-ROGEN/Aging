#!/bin/bash
# Helper script to install Graphviz on macOS

set -e

echo "Installing Graphviz..."

# Check for Homebrew
if command -v brew &> /dev/null; then
    echo "Found Homebrew, installing Graphviz..."
    brew install graphviz
    echo "Graphviz installed successfully!"
    exit 0
fi

# Check if Homebrew exists but not in PATH
if [ -f "/opt/homebrew/bin/brew" ]; then
    echo "Found Homebrew at /opt/homebrew/bin/brew"
    /opt/homebrew/bin/brew install graphviz
    echo "Graphviz installed successfully!"
    exit 0
fi

if [ -f "/usr/local/bin/brew" ]; then
    echo "Found Homebrew at /usr/local/bin/brew"
    /usr/local/bin/brew install graphviz
    echo "Graphviz installed successfully!"
    exit 0
fi

echo "Homebrew not found. Please install Homebrew first:"
echo "  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
echo ""
echo "Then run this script again, or install Graphviz manually:"
echo "  brew install graphviz"
exit 1
