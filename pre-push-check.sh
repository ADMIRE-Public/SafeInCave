#!/bin/bash
# Pre-push verification script
# Run this before pushing to GitHub to catch CI failures early

set -e

echo "🔍 Running ruff linter..."
ruff check .

echo ""
echo "🧪 Running tests..."
python -m pytest tests/ -q

echo ""
echo "✅ All checks passed! Safe to push."
