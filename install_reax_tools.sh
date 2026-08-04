#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="$PROJECT_DIR/bin"
VENV_DIR="$PROJECT_DIR/venv"

echo "Installing ReaxTools 2.1"
echo "Project directory: $PROJECT_DIR"

echo ""
echo "[1/3] Checking build tools"
for tool in cmake g++ gcc; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Error: $tool not found. Please install cmake and a C++17 compiler first." >&2
        echo "  Debian/Ubuntu: sudo apt install cmake g++" >&2
        echo "  RHEL/Alibaba:  sudo dnf install cmake gcc-c++" >&2
        exit 1
    fi
done

echo ""
echo "[2/3] Building C++ analysis core"
cmake -S "$PROJECT_DIR" -B "$PROJECT_DIR/build" -DCMAKE_BUILD_TYPE=Release
cmake --build "$PROJECT_DIR/build" -j"$(nproc 2>/dev/null || echo 4)"

echo ""
echo "[3/3] Configuring Python environment"
find_python() {
    local candidate
    for candidate in python3.13 python3.12 python3.11 python3.10 python3.9 python3.8 python3; do
        if command -v "$candidate" >/dev/null 2>&1 && "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 8) else 1)' 2>/dev/null; then
            echo "$candidate"
            return 0
        fi
    done
    return 1
}

PY_BIN="$(find_python || true)"
if [ -z "$PY_BIN" ]; then
    echo "Error: Python 3.8 or newer is required but was not found." >&2
    echo "  Debian/Ubuntu: sudo apt install python3 python3-venv" >&2
    echo "  RHEL/Alibaba:  sudo dnf install python39 python39-pip" >&2
    echo "  Or activate a conda environment with Python >= 3.8 before running." >&2
    exit 1
fi
echo "Using Python: $("$PY_BIN" --version 2>&1) ($(command -v "$PY_BIN"))"

use_user_install=0
if [ ! -x "$VENV_DIR/bin/python" ]; then
    echo "Creating virtual environment at $VENV_DIR"
    if ! "$PY_BIN" -m venv "$VENV_DIR"; then
        echo "Warning: venv creation failed; falling back to --user install."
        use_user_install=1
    fi
fi

if [ "$use_user_install" -eq 0 ] && [ -x "$VENV_DIR/bin/python" ]; then
    VENV_PY="$VENV_DIR/bin/python"
    "$VENV_PY" -m pip install --upgrade pip >/dev/null 2>&1 || true
    "$VENV_PY" -m pip install -r "$BIN_DIR/requirements.txt"
    echo "Python dependencies installed into $VENV_DIR"
else
    "$PY_BIN" -m pip install --user -r "$BIN_DIR/requirements.txt"
    echo "Python dependencies installed with --user"
fi

echo ""
echo "[4/4] Registering command path"
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

chmod +x "$BIN_DIR/reax_tools" "$BIN_DIR/reax_tools.py" "$BIN_DIR/reax_tools_core" 2>/dev/null || true

echo ""
echo "ReaxTools installation completed."
echo "Command: reax_tools"
echo "Installed path: $BIN_DIR"
echo "Python environment: $VENV_DIR"
echo "Shell configuration updated: $shell_rc"
echo "Restart your terminal or run: source $shell_rc"
