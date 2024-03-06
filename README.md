# Blume-Capel-3D
This code is able to simulate the Blume-Capel Model in 3D.

It can also, if chosen as e.g. D = -1000, simulate the Ising model, since for D -> -inf, the model converges to the Ising Model.

The parameters, which the user should input in the console are:

#INT > 0

size: size of the (3D) lattice

#INT > 0

#Recommendation for testing: 10000

n_steps: Amount of Monte Carlo Steps

#INT > 0 && n_steps % N_Bin == 0

#Recommendation for testing: 100

N_bin: Amount of Bins, used in the Jackknife routine

#Choose whether you want the Ising Model or Blume-Capel model

BC: ß_c = 0.387721735, D = 0.655

I: ß_c = 0.22165463

If you choose the BC model, you can also decide to change the 'D'.

There is also the possibility to change the temperature to a high temperature (for example ß = 0.1).

Then the program will calculate the high-temperature expansion of the susceptibility of the BC model at D = 0.641

and give you the value & comparison to the measured value from the MC simulation at the end.

###BOUNDARY CONDITIONS

This program can only calculate the correlations for the two following cases:

###ONLY PERIODIC BOUNDARY CONDITIONS

This is when you choose 'P'.
Then all directions will have periodic boundary conditions and the Correlation function
inside the bulk, (bulk-bulk) will be calculated.

--> Choose this option if you want to calculate the high-temp-susceptibility-expansion comparison and it makes sense to always choose this for the Ising model as well.
Generally if you are interested in calculating the Binder Cumulant, choose this option.


###OPEN BOUNDARY CONDITIONS

This is when you choose 'o'. 
Here you can decide which directions should have periodic boundary conditions, and which not.

Correlations will only be calculated for the case of surfaces in +/- z direction. 

That means to DISABLE periodic boundary conditions in +/-z direction you enter the following when asked:

+x = 1

-x = 1

+y = 1

-y = 1

+z = 0

-z = 0

1 == PERIODIC boundary conditions are turned ON

0 == PERIODIC boundary conditions are turned OFF --> OPEN BOUNDARY CONDITIONS


#External Field

Do you want an external field h != 0?

Then enter Yes, and chose a value.


#Other Observables(Y/N):

Choose 'Y' if you want to calculate the Magnetization, Susceptibility and the Binder Cumulant.

This is on by default, if you calculate the high-temp-expansion...

------------------------------------------------------------------------------------------------------------------------------------

This program will always calculate correlations if the above described boundary conditions are chosen.

These may not always generate meaningful results if for example the temperature is changed, the magnetic field is on...!


