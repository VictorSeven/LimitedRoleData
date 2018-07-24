# LimitedRoleData

Source code and data for our work, "Limited role of spatial self-structuring in emergent trade-offs during pathogen evolution". You can read a pre-print version [in ArXiV](https://arxiv.org/abs/1804.08463).

## Abstract

  Pathogen transmission and virulence are main evolutionary variables broadly assumed to be linked through trade-offs. In well-mixed populations, these trade-offs are often ascribed to physiological restrictions, while populations with spatial self-structuring might evolve emergent trade-offs. Here, we reexamine a spatially-explicit, SIR model of the latter kind proposed by Ballegooijen and Boerlijst with the aim of characterising the mechanisms causing the emergence of the trade-off and its structural robustness. Using invadability criteria, we establish the conditions under which an evolutionary feedback between transmission and virulence mediated by pattern formation can poise the system to a critical boundary separating a disordered state (without emergent trade-off) from a self-structured phase (where the trade-off emerges), and analytically calculate the functional shape of the boundary in a certain approximation. Beyond evolutionary parameters, the success of an invasion depends on the size and spatial structure of the invading and invaded populations. Spatial self-structuring is often destroyed when hosts are mobile, changing the evolutionary dynamics to those of a well-mixed population. In a metapopulation scenario, the systematic extinction of the pathogen in the disordered phase may counteract the disruptive effect of host mobility, favour pattern formation and therefore recover the emergent trade-off.

## How to use the information on this repository

### C++ Code

Here we present all the source code and data we used in our research. 

The file `virus.cpp` contains the source code of the simulations. It must be compiled for C++11 standard. When you execute the code, you must introduce flags as parameters:

```
./virus.cpp its beta ti random_seed D_counter file_name frames recover 
```

where:

- `its` is the number of iterations of the code. In order to get an evolutionary trajectory you should use values above 10^6.
- `beta` and `ti` are the values of the initial transmission rate and infection time, respectively.
- `random_seed` is the number passed to the RNG.
- `D_counter` controls the diffusion coefficient, which is D=1/D_counter. 
- `file_name` is the name of the output file(s).
- `frames` is the number of outputs you get when doing a video (note that each frame will have 4 files containing a LxL matrix. A big number of frames may produce a huge output! Use carefully!)
- `recover` is an index used to recover the simulation starting from a previous frame. In this case `file_name` should match the name of the existent simulation. If you don't want to recover a previous simulation, set this to -1. For example, setting recover to 50 and `file_name` to `test` would make the simulation start from file `video_test_50.txt` (which was generated automatically if you let the program to write the video).

Take in account that for performance reasons, videos are disabled by default. You need to uncomment the appropiate function in the main function in order to write them. The same applies to diffusion. The code is fully commented, and these functions are marked for ease of use.

### Graphs

Graphs were generated using a Jupyter notebook, with Python 2.7 and the Anaconda distribution. It may work under Python 3 also. In order to produce the outputs by yourself, it is necessary to downloaded and decompress the `data.zip` folder. Put the `graphs.ipynb` at the same level than the decompressed folder. Then you will be able to produce the graphs by yourself.

## Further information

For any other inquiry about the code, you can reach me at vbuendiar #at# onsager #dot# ugr #dot es. 






