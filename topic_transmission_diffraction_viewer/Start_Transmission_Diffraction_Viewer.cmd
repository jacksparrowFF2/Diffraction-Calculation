@echo off
setlocal
title Transmission Diffraction Viewer
pushd "%~dp0"
if errorlevel 1 goto :directory_error

set "VIEWER_PYTHON_ARGS="
set "VIEWER_PYTHON_CHECK=import sys; import streamlit, h5py, numpy, matplotlib, PIL; from packaging.version import Version; assert sys.version_info >= (3, 10) and Version(streamlit.__version__) >= Version('1.65')"

for %%P in (".venv\Scripts\python.exe" "..\.venv\Scripts\python.exe") do (
    call :try_python "%%~P"
    if not errorlevel 1 goto :run
)

for /f "delims=" %%P in ('where.exe python.exe 2^>nul') do (
    call :try_python "%%P"
    if not errorlevel 1 goto :run
)

rem Explorer can start with a different PATH from an activated Conda shell.
for %%D in ("%ProgramData%" "%LocalAppData%" "%UserProfile%") do (
    for %%C in (Anaconda3 Miniconda3) do (
        call :try_python "%%~D\%%C\python.exe"
        if not errorlevel 1 goto :run
    )
)

set "VIEWER_PYTHON=py.exe"
set "VIEWER_PYTHON_ARGS=-3"
"%VIEWER_PYTHON%" %VIEWER_PYTHON_ARGS% -c "%VIEWER_PYTHON_CHECK%" >nul 2>&1
if not errorlevel 1 goto :run

echo Could not find a Python environment with the viewer dependencies.
echo Install the dependencies in requirements.txt.
popd
pause
exit /b 1

:run
echo Starting Transmission Diffraction Viewer. Close this window to stop the server.
echo Python: "%VIEWER_PYTHON%" %VIEWER_PYTHON_ARGS%
rem Leave the port unset so Streamlit can retry an occupied default port.
set "STREAMLIT_SERVER_PORT="
set "STREAMLIT_BROWSER_SERVER_PORT="
"%VIEWER_PYTHON%" %VIEWER_PYTHON_ARGS% -m streamlit run "%~dp0src\streamlit_app.py" --server.address 127.0.0.1 --server.headless false --server.showEmailPrompt false --browser.serverAddress 127.0.0.1 --browser.gatherUsageStats false
set "VIEWER_EXIT=%errorlevel%"
popd
if not "%VIEWER_EXIT%"=="0" pause
exit /b %VIEWER_EXIT%

:try_python
if /i "%~dp1"=="%LocalAppData%\Microsoft\WindowsApps\" exit /b 1
if not exist "%~1" exit /b 1
"%~1" -c "%VIEWER_PYTHON_CHECK%" >nul 2>&1
if errorlevel 1 exit /b 1
set "VIEWER_PYTHON=%~f1"
exit /b 0

:directory_error
echo Could not open the viewer directory.
pause
exit /b 1
