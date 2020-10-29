To repeat the yasso15 one month experiment
the following files are needed:

A)
1. awen_monthly.txt 
2. weather_monthly.txt

The monthly litter fall and monthly climate data.

B)
1. y15_subroutinen.f90
2. y15_subroutine_month.f90

yasso15  implementation in  fortran, the  latter one  has change  on a
single line for one month. The loop for average temperature dependence
is replaced with a single call to one month temperature dependence.

C)
1. load_yasso15.r
2. load_yasso15_month.r

Load compiled fortran yasso implementations.

D)
1. y15month_test.r

To run the test install all the files in the same directory. Compile
the fortran files (usually from command line: R CMD SHLIB <files>).
Finally in R, run y15month_test.r.

The "init.x0.2015" is so called steady state. Yasso15 is first run to
steady state before the test runs are made for monthly time step.

The output files is Yasso15ResultsMonth.txt that contains results for two cases.
1) the monthly temperature is from the sine function (see the beginning of y15month_test.r)
2) the monthly temperature is from the data 'weather_monthly.txt'.
In both cases the monthly litter fall is from awen_monthly.txt. 

E)
1. Yasso15ResultsMonth.xlsx
Excel file from  Yasso15ResultsMonth.txt. 
The comparison with one year time step is made.
