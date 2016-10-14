%-------------------------%
% Startup file for octave % 
%-------------------------%

%--- list additional packages installed
disp(" additional packages installed ")
disp("-------------------------------")
disp("                               ")
pkg list
disp("                               ")
disp("                               ")

%--- load packages
disp(" load package| netcdf ")
pkg load netcdf
disp(" load package| odepkg ")
pkg load odepkg
