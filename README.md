# wavelet-sbl
ECE 251C Project testing Wavelet Domain Processing in Sparse Bayesian Learning Framework. This of course is a very generic project statement, so below is a more in depth description of the project.

## Imports
A handful of packages have been added. A (potentially out of date) list is below:

```julia
Match.jl (Pkg.add('Match')) 
Parameters.jl (Pkg.add('Parameters'))
MAT.jl (Pkg.add('MAT'))
Plots.jl (Pkg.add('Plots'))
Revise.jl (Pkg.add('Revise'))
LinearAlgebra.jl (Pkg.add('LinearAlgebra'))
Peaks.jl (Pkg.add('Peaks'))


Pkg.add("Match"); Pkg.add("Parameters"); Pkg.add("MAT"); Pkg.add("Plots"); Pkg.add("LinearAlgebra"); Pkg.add("Peaks")
```
## Data
The data itself is from the SWellEx-96 Experiment. More information on the experiment can be found [here](https://swellex96.ucsd.edu/). The general gist is that we're working with a few array geometries. The one we're using here is a Vertical Line Array (VLA). Nine different tones were transmitted underwater to this array over 20 minutes, and the data we have is for 64 hydrophones (63 officially because one of the hydrophones was broken).

## Project
This project can be nicely broken into four different stages, and we'll go through each stage one by one.

### Stage 1 - Classical/Conventional Beamforming
This style of beamforming is the simplest for DoA (direction of arrival) source localization. We take the input data and get the output power correlation for each steering vector (a vector of expected delays based on the direction of the incoming signals). 

The direction of arrival should just be the direction with the highest power/correlation value. This is similar to an autocorrelation function calculation - for a sequence $r_{xx}(n)$, the largest value must be $r_{xx}(0)$ because that is the highest correlation between the two signals. 

The issue with conventional beamforming is that it will always struggle to measure multiple sources.