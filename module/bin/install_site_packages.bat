@echo off
set fpath=%cd%
set python=%fpath%/python/python.exe
call %python% %fpath%/get-pip.py
call %python% -m pip install --upgrade pip
call %python% -m pip install --upgrade pip
call %python% -m pip install numpy
call %python% -m pip install matplotlib
call %python% -m pip install pandas
call %python% -m pip install h5py
call %python% -m pip install netcdf4
pause