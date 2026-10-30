@echo off
REM ============================================================
REM  Chemora Extract - one-click local run for Windows
REM  Double-click this file, or run it from the VS Code terminal.
REM  It creates a conda env the first time, then launches the app.
REM ============================================================
cd /d "%~dp0"

REM -- Does the 'chemora' conda env already exist? --
call conda env list | findstr /C:"chemora" >nul
if errorlevel 1 (
    echo.
    echo [setup] First run - creating conda env "chemora". This takes a few minutes...
    REM Python 3.9: OpenChemIE's pinned stack (torch 1.13, Pillow 9.5, opencv 4.5)
    REM has clean 3.9 wheels. Do NOT use 3.11+.
    call conda create -y -n chemora python=3.9
    call conda activate chemora
    echo [setup] Installing poppler + pdftotext + rdkit via conda-forge...
    REM pdftotext is an OpenChemIE dependency that compiles against poppler and
    REM fails to pip-build on Windows; conda-forge ships a prebuilt wheel.
    call conda install -y -c conda-forge poppler pdftotext rdkit
    echo [setup] Installing PyTorch 1.13.1 (CPU build) BEFORE OpenChemIE...
    REM OpenChemIE pins torch>=1.10,<2.0 — install 1.13.1 first so pip can't pull a 2.x.
    call pip install torch==1.13.1 --index-url https://download.pytorch.org/whl/cpu
    echo [setup] Installing OpenChemIE from GitHub source (not on PyPI)...
    REM Installed before the app deps so its Pillow==9.5.0 / opencv==4.5.5.64 pins win.
    call pip install "OpenChemIE @ git+https://github.com/CrystalEye42/OpenChemIE"
    echo [setup] Installing app deps (no pillow pin - let OpenChemIE's 9.5.0 stand)...
    call pip install pymupdf streamlit streamlit-ketcher
) else (
    call conda activate chemora
)

REM -- Engine: openchemie = molecules + reactions. Change to "manual" for a no-model UI demo. --
set CHEMORA_ENGINE=openchemie

echo.
echo [run] Launching Chemora Extract... your browser will open shortly.
streamlit run app.py
pause
