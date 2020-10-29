To repeat the yasso07 one month experiment
the following files are needed:

A)
1. awen_monthly.txt 
2. weather_monthly.txt

The monthly litter fall and climate data.

B)
3. y07_subroutinen.f90
4. y07_subroutine_month.f90

yasso07 implementation in fortran, the latter  one has change
on a single line (line 39) for one month. The loop
is replaced with a single call to temperature dependence.

C)
5. load_yasso07.r
6. load_yasso07_month.r

Load compiled fortran yasso implementations.

D)
7. y07test2015.r

To run the test install all the files in the same directory. Compile the fortran files.
Finally in R, run y07test2015.r. The initial state of yasso07 is in the line 48. The default
is so called steady state. Comment out the line 48 and the initial state will be zero
(AWEN=(0,0,0,0,0)).

The output files are OneMonthStep2015.txt that runs one year in 12 steps
and Yasso07ResultsMonth.txt that contains results for two cases.
1) the monthly temperature is from the sine function (see lines 3-16 in y07test2015.r and Yasso07.pdf)
2) the monthly temperature is from the data 'weather_monthly.txt'.


Yasso07Month2015.xlsx
Excel file from y07test2015.r that contains the two 
cases mentioned above for two initial conditions: 1) the initial state is zero (i.e. AWEN = (0,0,0,0,0) and 2) the initial state is so called steady state.
The comparison with one year time step is made.
