module constants 
use iso_fortran_env, only: real64
implicit none 

        ! real(kind=8), parameter :: Kb = 1.3807e-23 
        ! real(kind=8), parameter :: gyromagnetic_ratio = 1.76e11
        ! real(kind=8), parameter :: Bohr_magneton = 9.2740e-24

        real(kind=8), parameter :: Kb = 8.6173303e-5
        real(kind=8), parameter :: gyromagnetic_ratio = 1.76e11
        real(kind=8), parameter :: Bohr_magneton = 5.7883818060e-5 ! ev/ T
        real(kind=8), parameter :: pi = 3.14159265358979323846 
        integer, parameter :: dp = real64
        public :: Kb, gyromagnetic_ratio, Bohr_magneton, dp
end module constants 
