$InstallPath = $HOME\AppData\Local\Executables

if (-not (Test-Path $InstallPath)) { mkdir $InstallPath }

pyinstaller -wF pylinegui.py
