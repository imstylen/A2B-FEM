@echo off

echo %cd%
pause

C:\Octave\octave.vbs -p "%cd%" --no-gui --persist --eval run("'Dependencies\OCTAVE\RunAnalysis.m'")

