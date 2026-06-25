#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$PROJECT_DIR/bin"

echo "Installing ReaxTools 2.1"
echo "Project directory: $PROJECT_DIR"

echo ""
echo "[1/3] Building C++ analysis core"
cmake -S "$PROJECT_DIR" -B "$PROJECT_DIR/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$PROJECT_DIR/build" -j"$(nproc 2>/dev/null || echo 4)"

echo ""
echo "[2/3] Installing Python plotting dependencies"
if command -v python3 >/dev/null 2>&1; then
    python3 -m pip install --user -r "$BIN_DIR/requirements.txt"
else
    echo "Warning: python3 was not found. Analysis will work, but plotting commands need Python."
fi

echo ""
echo "[3/3] Registering command path"
current_shell="$(basename "${SHELL:-bash}")"
case "$current_shell" in
    bash) shell_rc="$HOME/.bashrc" ;;
    zsh) shell_rc="$HOME/.zshrc" ;;
    ksh) shell_rc="$HOME/.kshrc" ;;
    csh) shell_rc="$HOME/.cshrc" ;;
    tcsh) shell_rc="$HOME/.tcshrc" ;;
    *) shell_rc="$HOME/.profile" ;;
esac

touch "$shell_rc"
export_path_line="export PATH=\"\$PATH:$BIN_DIR\""
temp_file="$(mktemp)"
grep -v "ReaxTools installation" "$shell_rc" | grep -v "$BIN_DIR" > "$temp_file" 2>/dev/null || true
mv "$temp_file" "$shell_rc"

{
    echo ""
    echo "# ReaxTools installation"
    echo "$export_path_line"
} >> "$shell_rc"

chmod +x "$BIN_DIR/reax_tools" "$BIN_DIR/reax_tools_core" 2>/dev/null || true

echo ""
echo "ReaxTools installation completed."
echo "Command: reax_tools"
echo "Installed path: $BIN_DIR"
echo "Shell configuration updated: $shell_rc"
echo "Restart your terminal or run: source $shell_rc"
