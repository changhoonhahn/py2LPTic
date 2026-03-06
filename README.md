# modified version of 2LPTic
modified version of Roman's 2LPTic code, which can be found in
https://cosmo.nyu.edu/roman/2LPT/, used specifically for running the Quijote
simulations 

## Quijote modifications 
These are modifications made specifically for Quijote: 

```
This code can be used to generate the IC for N-body simulations using 2LPT.
We have downloaded the code from Roman Scoccimarro web page and performed
a few changes to use in the standard way we do.

The first modifications are these:

1) We have left
  
  Dplus = GrowthFactor(InitTime, 1.0);
  
in power.c instead of setting:

   Dplus = 1.0;

In this way we need to input the z=0 P(k) and set sigma8 to its value at z=0.
The code will compute the growth factor at the starting redshift and change
the amplitude by itself.



2) It seems that the code is designed to read a file containing 

log10(k), log10(delta2). Therefore, we have modified in power.c:

  PowerTable[NPowerTable].logk = k;
  PowerTable[NPowerTable].logD = p;

by

  PowerTable[NPowerTable].logk = log10(k);
  PowerTable[NPowerTable].logD = log10(p * 4 * M_PI * k * k * k);



3) From the way it computes the value of sigma_8 (see function sigma2_int
in power.c) there are some extra factors. To correct for this we have added
to power.c in line 235

   P = P / (8.0 * M_PI * M_PI * M_PI); //Normalized factor


4) We have modified the routines

-allvars.h
-allvars.c
-read_param.c
-main.c

to allow the user switching on/off the Rayleigh sampling of the modes





The code has the functionality to output the linear density field used to
generate the initial conditions. For this switch on the OUTPUT_DF flag in
the Makefile. The output will be three different files

-Coordinates_ptype_X.Y
-Amplitudes_ptype_X.Y
-Phases_ptype_X.Y

where X can be 0 (gas),1 (cdm) and 2 (neutrinos), and Y goes from 0 to the
effective number of cpus used when running the code

All the files have a header:
Nfiles(int32) Nmesh(int32) Local_nx(int32)

The header is followed by the data from Local_nx*Nmesh*(Nmesh/2 + 1) modes
In Coordinates is np.int64 and in Amplitudes and Phases is np.float32

Coordinates is a single long integer containing the position of the mode:
Coordinate = (i*Nmesh + j)*(Nmesh/2 + 1) + k, where i and j go from 0 to Nmesh
and k goes from 0 to Nmesh/2 + 1

Amplitudes contains a float with the value of the amplitude (A from A*exp(i*theta))
The units can be obtained by multiplying its number by BoxSize**(3.0/2.0)

Phases contains a float with the value of the phase (no units)

The Pk_library.pyx of the Pylians library contain a routine (Pk_NGenIC_IC_field)
to compute the Pk from these files
```
