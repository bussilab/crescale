Modified SimpleMD version. 

Compared with the original version hosted [here](https://github.com/GiovanniBussi/simplemd):
- Only c++ implementation
- Required PLUMED library 
- Added stochastic velocity rescaling thermostat
- Added stochastic cell rescaling barostat (isotropic and anisotropic implementations)

Thermostat can be activated with keyword `svr`. Notice that tau is computed as `1/(2*friction)`.

Barostat can be activated with keyword `barotype X`, where
- X=`I`: Isotropic barostat, Trotter-based integrator;
- X=`A`: Anisotropic barostat, Euler integrator (no energy drift available);
- X=`T`: Anisotropic barostat, time-reversible integrator;

The BAOAB - limit case implementation can be activated with the keyword `baoab`.

The barostat can be applied with a multiple time step framework. The stride can be specified with the keyword `nstbaro N` where `N` is the stride (1 by default).

By default, the barostat does NOT scale velocities. The version scaling velocities can be triggered with the keyword `baroscalev T`.

In the fully anisotropic case, box rotations can be activated with the keyword `rotations T` (by default they are eliminated). 

The keyword `pressure` is followed by the 9 components of the external stress tensor, in the following order: `xx,yy,zz,xy,xz,yx,yz,zx,zy`. In the isotropic case it is possible to use only one value, and in any case only the first element is read.

The input file is in gro format. The box matrix at the end of the file must be lower triangular if the anisotropic barostat is activated without rotations, namely: `xx,yy,zz,0,0,yx,0,zx,zy`.

The columns of the statistics file refer to:
1.  MD step;
2.  time; 
3.  instantaneous temperature; 
4.  potential energy;
5.  conserved energy; 
6.  volume; 
7.  virial; 
8.  box(xx); 
9.  box(yy); 
10. box(zz);
11. box(xy);
12. box(xz); 
13. box(yx);
14. box(yz);
15. box(zx); 
16. box(zy).

Entries 8-16 are only present in the anisotropic case.

A sample input is below:
```
inputfile crystal.gro
outputfile output.gro
temperature 0.1
tstep 0.005
friction 10.0
forcecutoff 2.5
listcutoff  3.0
nstep 10000000
nconfig 100000 trajectory.gro
nstat   10 energies.dat
pressure 1.0 1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0
taup 1.0
nstbaro 2
barotype A
baroscalev T
baoab T
rotations F
svr T

```
