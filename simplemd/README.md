Modified SimpleMD version.

Compared with the original version hosted [here](https://github.com/GiovanniBussi/simplemd):
- Only c++ implementation
- Added stochastic velocity rescaling thermostat
- Added stochastic cell rescaling barostat

Thermostat can be activated with keyword `svr`. Notice that tau is computed as `1/(2*friction)`.

Barostat can be activated with keyword `barotype X` where `X=E,R,T` correspond to one of the integrators described in the manuscript.

The barostat can be applied with a multiple time step framework. The stride can be specified with the keyword `nstbaro N` where `N` is the stride (1 by default).

By default, the barostat does NOT scale velocities. The version scaling velocities can be triggered with the keyword `baroscalev T`.

A sample input is below:
```
inputfile liquid.xyz
outputfile output.xyz
temperature 1.5
tstep 0.005
friction 10.0
forcecutoff 2.5
listcutoff  3.0
nstep 10000000
nconfig 1000000 trajectory.xyz
nstat   10 new_energies_1.0_T_T_2.dat
pressure 1.0
taup 1.0
nstbaro 2
barotype T
baroscalev T
svr true
```
