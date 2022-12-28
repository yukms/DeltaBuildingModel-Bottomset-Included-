# DeltaBuildingModel-Bottomset-Included-
1D delta building model with bottomset

The effect of bottomset on fluviodeltaic land‐building process: Numerical modelling and physical experiment (2022)

**Authors**: Minsik Kim, Wonsuck Kim, and Wook-Hyun Nahm

**Corresponding Author**: Wonsuck Kim (delta@yonsei.ac.kr)
* This model is computing Delta building process. Especially bottomset is formed by advection settling model


inputdata_Loop.m
--------
This define input parameters and run the main calculation (fbt_delta_settling_Loop.m)

fbt_delta_settling_Loop.m
--------
This basic model is responsible for running repeated delta building under a given set of input conditions (water depth and sand/mud ratio).


inputdata_fieldapp.m
--------
This define input parameters and run the main calculation (fbt_delta_settling_fieldapp.m)


fbt_delta_settling_fieldapp.m
--------
This is the model script, responsible for field application to Wax Lake Delta, under a given set of retention rate change from 0 to 0.5
