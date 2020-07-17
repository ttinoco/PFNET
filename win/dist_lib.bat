call clean.bat
set PFNET=pfnet-1.3.6
set DIR=%CD%
7z a -ttar %PFNET%.tar -xr!.git -xr!python %DIR%
7z a -tgzip %PFNET%.tar.gz %PFNET%.tar
del %PFNET%.tar
