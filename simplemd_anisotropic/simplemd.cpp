#include "Random.h"
#include <string>
#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <limits>
#include <plumed/tools/Vector.h>
#include <plumed/tools/Tensor.h>
#include <plumed/tools/Pbc.h>

using namespace std;
using namespace PLMD;


class SimpleMD
{

int iv[32];
int iy;
int iset;
double gset;
bool write_positions_first;
bool write_statistics_first;
int write_statistics_last_time_reopened;
FILE* write_statistics_fp;

public:
SimpleMD(){
  for(int i=0;i<32;i++) iv[i]=0.0;
  iy=0;
  iset=0;
  gset=0.0;
  write_positions_first=true;
  write_statistics_first=true;
  write_statistics_last_time_reopened=0;
  write_statistics_fp=NULL;
}

private:

void 
read_input(FILE*   fp,
           double& temperature,
           Tensor& pressure,
           double& tstep,
           double& friction,
           double& taup,
           double& forcecutoff,
           double& listcutoff,
           int&    nstep,
           int&    nconfig,
           int&    nstat,
           int&    nbox,
           bool&   svr,
           bool&   wrapatoms,
           int&    nstbaro,
           char&   barotype,
           bool&   baroscalev,
           bool&   baoab,
           bool&   rotations,
           string& inputfile,
           string& outputfile,
           string& trajfile,
           string& statfile,
           int&    maxneighbours,
           int&    idum)
{
  temperature=1.0;
  for(int i=0; i<3; i++) for(int j=0; j<3; j++) pressure[i][j]=0.0;
  tstep=0.005;
  friction=0.0;
  taup=0.0;
  forcecutoff=2.5;
  listcutoff=3.0;
  nstep=1;
  nconfig=10;
  nstat=1;
  nbox=1;
  svr=false;
  maxneighbours=1000;
  idum=0;
  wrapatoms=false;
  nstbaro=1;
  barotype='I';
  baroscalev=false;
  baoab=false;
  rotations=false;
  statfile="";
  trajfile="";
  outputfile="";
  inputfile="";

  string line;

  line.resize(256);
  char buffer[256];
  char buffer1[256];

  while(fgets(buffer,256,fp)){
    line=buffer;
    for(int i=0;i<line.length();++i) if(line[i]=='#' || line[i]=='\n') line.erase(i);
    for(int i=line.length()-1;i>=0;--i){
      if(line[i]!=' ')break;
      line.erase(i);
    }
    sscanf(line.c_str(),"%s",buffer);
    string keyword=buffer;
    if(keyword=="temperature")
      sscanf(line.c_str(),"%s %lf",buffer,&temperature);
    else if(keyword=="pressure")
      sscanf(line.c_str(),"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf",buffer,&pressure[0][0],
        &pressure[1][1],&pressure[2][2],&pressure[0][1],&pressure[0][2],
        &pressure[1][0],&pressure[1][2],&pressure[2][0],&pressure[2][1]);
    else if(keyword=="tstep")
      sscanf(line.c_str(),"%s %lf",buffer,&tstep);
    else if(keyword=="friction")
      sscanf(line.c_str(),"%s %lf",buffer,&friction);
    else if(keyword=="taup")
      sscanf(line.c_str(),"%s %lf",buffer,&taup);
    else if(keyword=="forcecutoff")
      sscanf(line.c_str(),"%s %lf",buffer,&forcecutoff);
    else if(keyword=="listcutoff")
      sscanf(line.c_str(),"%s %lf",buffer,&listcutoff);
    else if(keyword=="nstep")
      sscanf(line.c_str(),"%s %d",buffer,&nstep);
    else if(keyword=="nconfig")
    {
      sscanf(line.c_str(),"%s %d %s",buffer,&nconfig,buffer1);
      trajfile=buffer1;
    }
    else if(keyword=="nstat")
    {
      sscanf(line.c_str(),"%s %d %s",buffer,&nstat,buffer1);
      statfile=buffer1;
    }
    else if(keyword=="nstbaro")
    {
      sscanf(line.c_str(),"%s %d",buffer,&nstbaro);
    }
    else if(keyword=="barotype")
    {
      sscanf(line.c_str(),"%s %c",buffer,&barotype);
    }
    else if(keyword=="baroscalev")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      if(buffer1[0]=='T' || buffer1[0]=='t') baroscalev=true;
    }
    else if(keyword=="baoab")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      if(buffer1[0]=='T' || buffer1[0]=='t') baoab=true;
    }
    else if(keyword=="rotations")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      if(buffer1[0]=='T' || buffer1[0]=='t') rotations=true;
    }
    else if(keyword=="wrapatoms")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      if(buffer1[0]=='T' || buffer1[0]=='t') wrapatoms=true;
    }
    else if(keyword=="svr")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      if(buffer1[0]=='T' || buffer1[0]=='t') svr=true;
    }
    else if(keyword=="maxneighbours")
      sscanf(line.c_str(),"%s %d",buffer,&maxneighbours);
    else if(keyword=="inputfile")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      inputfile=buffer1;
    }
    else if(keyword=="outputfile")
    {
      sscanf(line.c_str(),"%s %s",buffer,buffer1);
      outputfile=buffer1;
    }
    else if(keyword=="idum")
      sscanf(line.c_str(),"%s %d",buffer,&idum);
    else{
      fprintf(stderr,"Unknown keywords :%s\n",keyword.c_str());
      exit(1);
    }
  }

  if(inputfile.length()==0){
      fprintf(stderr,"Specify input file\n");
      exit(1);
  }
  if(outputfile.length()==0){
      fprintf(stderr,"Specify output file\n");
      exit(1);
  }
  if(trajfile.length()==0){
      fprintf(stderr,"Specify traj file\n");
      exit(1);
  }
  if(statfile.length()==0){
      fprintf(stderr,"Specify stat file\n");
      exit(1);
  }
}

void read_natoms(const string& inputfile,int & natoms){
// read the number of atoms in file "input.gro"
  FILE* fp=fopen(inputfile.c_str(),"r");
  char buffer[256];
  fgets(buffer,256,fp);
  fscanf(fp,"%d",&natoms);
  fclose(fp);
}

void read_positions(const string& inputfile,int natoms,vector<Vector>& positions,Tensor& box,char& barotype,bool rotations){
// read positions and box from a file called inputfile
// natoms (input variable) and number of atoms in the file should be consistent
  FILE* fp=fopen(inputfile.c_str(),"r");
  char buffer1[256], buffer2[5], resname[5], atomname[5];
  int atomnumber;
  fgets(buffer1,256,fp);
  fgets(buffer2,5,fp);
  for(int i=0;i<natoms;i++){
    fscanf(fp,"%s %s %d %lf %lf %lf",resname,atomname,&atomnumber,&positions[i][0],&positions[i][1],&positions[i][2]);
// note: resname, atomname and atombumber are read but not used
  }
  fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&box[0][0],&box[1][1],&box[2][2],&box[0][1],
                                       &box[0][2],&box[1][0],&box[1][2],&box[2][0],&box[2][1]);
  if(barotype == 'A' && rotations == 'F' &&(box[0][1] != 0.0 || box[0][2] != 0.0 || box[1][2] != 0.0)){
    fprintf(stderr,"Error: use a lower-triangular box matrix at input if barotype == 'A' and rotations == 'F'.\n");
    exit(1);
  }
  fclose(fp);
}

void randomize_velocities(const int natoms,const double temperature,const vector<double>& masses,vector<Vector>& velocities,Random& random){
// randomize the velocities according to the temperature
  for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<3;i++)
      velocities[iatom][i]=sqrt(temperature/masses[iatom])*random.Gaussian();
}

void reset_com_velocity(const int natoms,const vector<double>&masses,vector<Vector>& velocities){
  Vector cm;
  double m;
  for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<3;i++){
    cm[i]+=velocities[iatom][i]*masses[iatom];
    m+=masses[i];
  }
  for(int i=0;i<3;i++) cm[i]/=m;
  for(int iatom=0;iatom<natoms;iatom++) for(int i=0;i<3;i++)
    velocities[iatom][i]-=cm[i];
}

void check_list(const int natoms,const vector<Vector>& positions,const vector<Vector>& positions0,const double listcutoff,
                const double forcecutoff,bool & recompute)
{
// check if the neighbour list have to be recomputed
  Vector displacement;  // displacement from positions0 to positions
  double delta2;        // square of the 'skin' thickness
  recompute=false;
  delta2=(0.5*(listcutoff-forcecutoff))*(0.5*(listcutoff-forcecutoff));
// if ANY atom moved more than half of the skin thickness, recompute is set to .true.
  for(int iatom=0;iatom<natoms;iatom++){
    for(int k=0;k<3;k++) displacement[k]=positions[iatom][k]-positions0[iatom][k];
    double s=0.0;
    for(int k=0;k<3;k++) s+=displacement[k]*displacement[k];
    if(s>delta2) recompute=true;
  }
}


void compute_list(const int natoms,const int listsize,const vector<Vector>& positions,const Tensor& box,const double listcutoff,
                  vector<int>& point,vector<int>& list){
// see Allen-Tildesey for a definition of point and list
  Vector distance;     // distance of the two atoms
  Vector distance_pbc; // minimum-image distance of the two atoms
  double listcutoff2;  // squared list cutoff
  listcutoff2=listcutoff*listcutoff;
  point[0]=0;
  Pbc pbc;
  pbc.setBox(box);
  for(int iatom=0;iatom<natoms-1;iatom++){
    point[iatom+1]=point[iatom];
    for(int jatom=iatom+1;jatom<natoms;jatom++){
      distance_pbc=pbc.distance(positions[jatom],positions[iatom]);
// if the interparticle distance is larger than the cutoff, skip
      double d2=0; for(int k=0;k<3;k++) d2+=distance_pbc[k]*distance_pbc[k];
      if(d2>listcutoff2)continue;
      if(point[iatom+1]>listsize){
// too many neighbours
        fprintf(stderr,"%s","Verlet list size exceeded\n");
        fprintf(stderr,"%s","Increase maxneighbours\n");
        exit(1);
      }
      list[point[iatom+1]]=jatom;
      point[iatom+1]++;
    }
  }
}

void compute_forces(const int natoms,const int listsize,const vector<Vector>& positions,const Tensor& box,
                    double forcecutoff,const vector<int>& point,const vector<int>& list,vector<Vector>& forces,
                    Tensor& virial, double& engconf)
{
  Vector distance;        // distance of the two atoms
  Vector distance_pbc;    // minimum-image distance of the two atoms
  double distance_pbc2;   // squared minimum-image distance
  double forcecutoff2;    // squared force cutoff
  Vector f;               // force
  double engcorrection;   // energy necessary shift the potential avoiding discontinuities

  forcecutoff2=forcecutoff*forcecutoff;
  engconf=0.0;
  virial*=0.0;
  Pbc pbc;
  pbc.setBox(box);
  for(int i=0;i<natoms;i++)for(int k=0;k<3;k++) forces[i][k]=0.0;
  engcorrection=4.0*(1.0/pow(forcecutoff2,6.0)-1.0/pow(forcecutoff2,3));
  for(int iatom=0;iatom<natoms-1;iatom++){
    for(int jlist=point[iatom];jlist<point[iatom+1];jlist++){
      int jatom=list[jlist];
      distance_pbc=pbc.distance(positions[jatom],positions[iatom]); 
      distance_pbc2=0.0; for(int k=0;k<3;k++) distance_pbc2+=distance_pbc[k]*distance_pbc[k];
// if the interparticle distance is larger than the cutoff, skip
      if(distance_pbc2>forcecutoff2) continue;
      double distance_pbc6=distance_pbc2*distance_pbc2*distance_pbc2;
      double distance_pbc8=distance_pbc6*distance_pbc2;
      double distance_pbc12=distance_pbc6*distance_pbc6;
      double distance_pbc14=distance_pbc12*distance_pbc2;
      engconf+=4.0*(1.0/distance_pbc12 - 1.0/distance_pbc6) - engcorrection;
      for(int k=0;k<3;k++) f[k]=2.0*distance_pbc[k]*4.0*(6.0/distance_pbc14-3.0/distance_pbc8);
// same force on the two atoms, with opposite sign:
      for(int k=0;k<3;k++) forces[iatom][k]+=f[k];
      for(int k=0;k<3;k++) forces[jatom][k]-=f[k];
      for(int ivir=0;ivir<3;ivir++)for(int jvir=0;jvir<3;jvir++) virial[ivir][jvir]+=f[ivir]*distance_pbc[jvir];
    }
  }
}

void compute_engkin(const int natoms,const vector<double>& masses,const vector<Vector>& velocities,Tensor& engkin)
{
// calculate the kinetic energy from the velocities
  engkin*=0.0;
  for(int iatom=0;iatom<natoms;iatom++){
    for(int j=0;j<3;j++)for(int k=0;k<3;k++) engkin[j][k]+=0.5*masses[iatom]*velocities[iatom][j]*velocities[iatom][k];
  }
}

Tensor compute_matrixexp(const Tensor& in)
// calculate the exponential matrix exp(in) using Padé approximation up to 6th order
{
  auto id=Tensor::identity();
  Tensor in2=matmul(in,in);
  Tensor in3=matmul(in2,in);
  return matmul(inverse(id - 0.5*in + 0.1*in2 - in3/120.0),id + 0.5*in + 0.1*in2 + in3/120.0);
}

void rotate_back_box(Tensor& box)
// rotate back the box to a lower triangular shape, eliminating rotations
{
  Tensor box_buffer=+box;
  box[0][0] = std::sqrt(box_buffer[0][0]*box_buffer[0][0] + box_buffer[0][1]*box_buffer[0][1] 
               + box_buffer[0][2]*box_buffer[0][2]);
  box[1][0] = (box_buffer[0][0]*box_buffer[1][0] + box_buffer[0][1]*box_buffer[1][1] 
               + box_buffer[0][2]*box_buffer[1][2]) / box[0][0];
  box[1][1] = std::sqrt(box_buffer[1][0]*box_buffer[1][0] + box_buffer[1][1]*box_buffer[1][1] 
               + box_buffer[1][2]*box_buffer[1][2] - box[1][0]*box[1][0]);
  box[2][0] = (box_buffer[0][0]*box_buffer[2][0] + box_buffer[0][1]*box_buffer[2][1] 
               + box_buffer[0][2]*box_buffer[2][2]) / box[0][0];
  box[2][1] = (box_buffer[1][0]*box_buffer[2][0] + box_buffer[1][1]*box_buffer[2][1] 
               + box_buffer[1][2]*box_buffer[2][2] - box[2][0]*box[1][0]) / box[1][1];
  box[2][2] = std::sqrt(box_buffer[2][0]*box_buffer[2][0] + box_buffer[2][1]*box_buffer[2][1] 
               + box_buffer[2][2]*box_buffer[2][2] - box[2][0]*box[2][0] - box[2][1]*box[2][1]);
  box[0][1] = 0.0;
  box[0][2] = 0.0;
  box[1][2] = 0.0;
}

double resamplekin_sumnoises(int nn,Random & random){
/*
  returns the sum of n independent gaussian noises squared
   (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
*/
  if(nn==0) {
    return 0.0;
  } else if(nn==1) {
    double rr=random.Gaussian();
    return rr*rr;
  } else if(nn%2==0) {
    return 2.0*random.Gamma(nn/2);
  } else {
    double rr=random.Gaussian();
    return 2.0*random.Gamma((nn-1)/2) + rr*rr;
  }
}

double resamplekin(double kk,double sigma,int ndeg,double taut,Random& random){
/*
  kk:    present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  sigma: target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  ndeg:  number of degrees of freedom of the atoms to be thermalized
  taut:  relaxation time of the thermostat, in units of 'how often this routine is called'
*/
  double factor,rr;
  if(taut>0.1){
    factor=exp(-1.0/taut);
  } else{
    factor=0.0;
  }
  rr = random.Gaussian();
  return kk + (1.0-factor)* (sigma*(resamplekin_sumnoises(ndeg-1,random)+rr*rr)/ndeg-kk)
            + 2.0*rr*sqrt(kk*sigma/ndeg*(1.0-factor)*factor);
}

void thermostat(const int natoms,const vector<double>& masses,const double dt,const double friction,
                const double temperature,vector<Vector>& velocities,double & engint,Random & random){
// Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
// it is a linear combination of old velocities and new, randomly chosen, velocity,
// with proper coefficients
  double c1,c2;
  c1=exp(-friction*dt);
  for(int iatom=0;iatom<natoms;iatom++){
    c2=sqrt((1.0-c1*c1)*temperature/masses[iatom]);
    for(int i=0;i<3;i++){
      engint+=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
      velocities[iatom][i]=c1*velocities[iatom][i]+c2*random.Gaussian();
      engint-=0.5*masses[iatom]*velocities[iatom][i]*velocities[iatom][i];
    }
  }
}

void gthermostat(const int natoms,const vector<double>& masses,const double dt,const double friction,
                const double temperature,vector<Vector>& velocities,double & engint,Random & random){
// Stochastic velocity rescaling
  Tensor engkin;
  compute_engkin(natoms,masses,velocities,engkin);
  double engkin_scalar=engkin[0][0]+engkin[1][1]+engkin[2][2];
  engint+=engkin_scalar;
  double r=sqrt(resamplekin(engkin_scalar,(3*natoms-3)*temperature*0.5,3*natoms-3,1.0/(2*dt*friction),random)/engkin_scalar);
  for(int iatom=0;iatom<natoms;iatom++){
    for(int i=0;i<3;i++){
      velocities[iatom][i]*=r;
    }
  }
  engint-=r*r*engkin_scalar;
}

void write_positions(const string& trajfile,int natoms,const int istep,const double tstep,
                     const vector<Vector>& positions,const Tensor& box,const bool wrapatoms)
{
// write positions on file trajfile (gro format)
// positions are appended at the end of the file
  vector<Vector> pos(1);
  FILE*fp;
  Pbc pbc;
  pbc.setBox(box);
  if(write_positions_first){
    fp=fopen(trajfile.c_str(),"w");
    write_positions_first=false;
  } else {
    fp=fopen(trajfile.c_str(),"a");
  }
  fprintf(fp,"t = %8.3f\n",istep*tstep);
  fprintf(fp,"%d\n",natoms);
  for(int iatom=0;iatom<natoms;iatom++){
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
    for(int k=0;k<3;k++) pos[0][k]=positions[iatom][k];
    if(wrapatoms) pbc.apply(pos);
    fprintf(fp,"%5dAr     Ar%5d%8.3f%8.3f%8.3f\n",iatom+1,iatom+1,pos[0][0],pos[0][1],pos[0][2]);
  }      
  fprintf(fp,"%10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n",box[0][0],
            box[1][1],box[2][2],box[0][1],box[0][2],box[1][0],box[1][2],box[2][0],box[2][1]);
  fclose(fp);
}

void write_final_positions(const string& outputfile,int natoms,const int nstep,const double tstep,
                           const vector<Vector>& positions,const Tensor& box,const bool wrapatoms)
{
// write positions on file outputfile (gro format)
  vector<Vector> pos(1);
  FILE* fp;
  Pbc pbc;
  pbc.setBox(box);
  fp=fopen(outputfile.c_str(),"w");
  fprintf(fp,"Final positions: t = %8.3f\n",nstep*tstep);
  fprintf(fp,"%d\n",natoms);
  for(int iatom=0;iatom<natoms;iatom++){
// usually, it is better not to apply pbc here, so that diffusion
// is more easily calculated from a trajectory file:
    for(int k=0;k<3;k++) pos[0][k]=positions[iatom][k];
    if(wrapatoms) pbc.apply(pos);
    fprintf(fp,"%5dAr     Ar%5d%8.3f%8.3f%8.3f\n",iatom+1,iatom+1,pos[0][0],pos[0][1],pos[0][2]);
  }      
  fprintf(fp,"%10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f\n",box[0][0],
            box[1][1],box[2][2],box[0][1],box[0][2],box[1][0],box[1][2],box[2][0],box[2][1]);
  fclose(fp);
}

void write_statistics(const string & statfile,const int istep,const double tstep, const int natoms,const Tensor engkin,
           const double engconf,const double engint,const double volume,const Tensor virial,const char barotype,const Tensor& box){
// write statistics on file statfile
  if(write_statistics_first){
// first time this routine is called, open the file
    write_statistics_fp=fopen(statfile.c_str(),"w");
    write_statistics_first=false;
  }
  if(istep-write_statistics_last_time_reopened>100){
// every 100 steps, reopen the file to flush the buffer
    fclose(write_statistics_fp);
    write_statistics_fp=fopen(statfile.c_str(),"a");
    write_statistics_last_time_reopened=istep;
  }
  double virialscalar=virial[0][0]+virial[1][1]+virial[2][2];
  double engkinscalar=engkin[0][0]+engkin[1][1]+engkin[2][2];
  fprintf(write_statistics_fp,"%d %f %f %f %f %f %f",
          istep,istep*tstep,2.0*engkinscalar/(3.0*natoms),engconf,engkinscalar+engconf+engint,volume,virialscalar);
  if(barotype=='A' || barotype=='T') fprintf(write_statistics_fp," %f %f %f %f %f %f %f %f %f",
                            box[0][0],box[1][1],box[2][2],box[0][1],box[0][2],box[1][0],box[1][2],box[2][0],box[2][1]);
  fprintf(write_statistics_fp,"\n");
}

public:
int main(FILE*in,FILE*out){
  int             natoms;        // number of atoms
  vector<Vector>  positions;     // atomic positions
  vector<Vector>  velocities;    // velocities
  vector<double>  masses;        // masses
  vector<Vector>  forces;        // forces
  Tensor          box;           // box
  Tensor          box_start;     // starting reference box
  Tensor          box_start_inv; // inverse of reference box

// neighbour list variables
// see Allen and Tildesey book for details
  int              listsize;     // size of the list array
  vector<int>      list;         // neighbour list
  vector<int>      point;        // pointer to neighbour list
  vector<Vector>   positions0;   // reference atomic positions, i.e. positions when the neighbour list

// input parameters
// all of them have a reasonable default value, set in read_input()
  double      tstep;             // simulation timestep
  double      temperature;       // temperature
  Tensor      pressure;          // external pressure (stress) tensor
  double      pressure_hydro;    // scalar external pressure
  double      friction;          // friction for Langevin dynamics (for NVE, use 0)
  double      taup;              // for barostat
  double      listcutoff;        // cutoff for neighbour list
  double      forcecutoff;       // cutoff for forces
  int         nstep;             // number of steps
  int         nconfig;           // stride for output of configurations
  int         nstat;             // stride for output of statistics
  int         nbox;              // stride for output of box coordinates
  int         maxneighbour;      // maximum average number of neighbours per atom
  int         idum;              // seed
  bool        wrapatoms;         // if true, atomic coordinates are written wrapped in minimal cell
  int         nstbaro;           // stride for barostat
  char        barotype;          // integrator for barostat: I (isotropic), A (anisotropic Euler), T (anisotropic time-reversible)
  bool        baroscalev;        // if true, scale velocities and use kinetic energy for internal pressure calculation
  bool        baoab;             // if true, barostat implemented in the BAOAB limit case
  bool        rotations;         // if true, rotations of the box are allowed
  bool        svr;               // use stochastic velocity rescaling instead of Langevin
  string      inputfile;         // name of file with starting configuration (gro)
  string      outputfile;        // name of file with final configuration (gro)
  string      trajfile;          // name of the trajectory file (gro)
  string      statfile;          // name of the file with statistics
  string      string;            // a string for parsing

  Tensor      engkin;            // kinetic energy
  double      engconf;           // configurational energy
  Tensor      virial;            // virial tensor 
  double      engint;            // integral for conserved energy in Langevin dynamics and SCR crescale barostat
  Tensor      pext;              // tensor diagonal external pressure

  bool recompute_list;           // control if the neighbour list have to be recomputed
  bool deviatoric;               // control if deviatoric correction is applied 
  double volume0;                // volume of reference initial box

  Random random;                              // random numbers stream
  vector<Tensor> random_pairs(2);             // random numbers for BAOAB limit case

  read_input(stdin,temperature,pressure,tstep,friction,taup,forcecutoff,
             listcutoff,nstep,nconfig,nstat,nbox,svr,wrapatoms,nstbaro,
             barotype,baroscalev,baoab,rotations,inputfile,outputfile,
             trajfile,statfile,maxneighbour,idum);

// number of atoms is read from file inputfile
  read_natoms(inputfile,natoms);

// write the parameters in output so they can be checked
  fprintf(stdout,"%s %s\n","Starting configuration           :",inputfile.c_str());
  fprintf(stdout,"%s %s\n","Final configuration              :",outputfile.c_str());
  fprintf(stdout,"%s %d\n","Number of atoms                  :",natoms);
  fprintf(stdout,"%s %f\n","Temperature                      :",temperature);
  fprintf(stdout,"%s %f %f %f %f %f %f %f %f %f\n",
                           "Pressure tensor                  :",pressure[0][0],pressure[1][1],pressure[2][2],
                   pressure[0][1],pressure[0][2],pressure[1][0],pressure[1][2],pressure[2][0],pressure[2][1]);
  fprintf(stdout,"%s %f\n","Time step                        :",tstep);
  fprintf(stdout,"%s %f\n","Friction                         :",friction);
  fprintf(stdout,"%s %f\n","Taup                             :",taup);
  fprintf(stdout,"%s %f\n","Cutoff for forces                :",forcecutoff);
  fprintf(stdout,"%s %f\n","Cutoff for neighbour list        :",listcutoff);
  fprintf(stdout,"%s %d\n","Number of steps                  :",nstep);
  fprintf(stdout,"%s %d\n","Stride for trajectory            :",nconfig);
  fprintf(stdout,"%s %s\n","Trajectory file                  :",trajfile.c_str());
  fprintf(stdout,"%s %d\n","Stride for statistics            :",nstat);
  fprintf(stdout,"%s %s\n","Statistics file                  :",statfile.c_str());
  fprintf(stdout,"%s %d\n","Max average number of neighbours :",maxneighbour);
  fprintf(stdout,"%s %d\n","Seed                             :",idum);
  fprintf(stdout,"%s %s\n","Are atoms wrapped on output?     :",(wrapatoms?"T":"F"));
  fprintf(stdout,"%s %s\n","Stochastic velocity rescaling    :",(svr?"T":"F"));

// Setting the seed
  random.setSeed(idum);

// Since each atom pair is counted once, the total number of pairs
// will be half of the number of neighbours times the number of atoms
  listsize=maxneighbour*natoms/2;

// allocation of dynamical arrays
  positions.resize(natoms);
  positions0.resize(natoms);
  velocities.resize(natoms);
  forces.resize(natoms);
  masses.resize(natoms);
  point.resize(natoms);
  list.resize(listsize);

// masses are hard-coded to 1
  for(unsigned i=0;i<natoms;++i) masses[i]=1.0;

// energy integral initialized to 0
  engint=0.0;

// positions and box coordinates are read from file inputfile
  read_positions(inputfile,natoms,positions,box,barotype,rotations);
  box_start=box;
  box_start_inv=inverse(box_start);

// velocities are randomized according to temperature
  randomize_velocities(natoms,temperature,masses,velocities,random);

// if using global thermostat, reset center of mass velocity
  if(svr) reset_com_velocity(natoms,masses,velocities);

// neighbour list are computed, and reference positions are saved
  compute_list(natoms,listsize,positions,box,listcutoff,point,list);

  fprintf(stdout,"List size: %d\n",point[natoms-1]);
  for(int iatom=0;iatom<natoms;++iatom) for(int k=0;k<3;++k) positions0[iatom][k]=positions[iatom][k];

// forces are computed before starting md
  compute_forces(natoms,listsize,positions,box,forcecutoff,point,list,forces,virial,engconf);

// scalar hydrostatic pressure 
  if(barotype=='I') pressure_hydro = pressure[0][0];
  else pressure_hydro = (pressure[0][0]+pressure[1][1]+pressure[2][2])/3.0; 
  pext=Tensor::identity()*pressure_hydro;

// deviatoric correction is activated if external pressure tensor is not proportional to the identity
  deviatoric=false;
  if(pressure[0][0]!=pressure[1][1] || pressure[1][1]!=pressure[2][2] || pressure[0][0]!=pressure[2][2] ||
     pressure[0][1]!=0.0 || pressure[0][2]!=0.0 || pressure[1][0]!=0.0 || pressure[1][2]!=0.0 || 
     pressure[2][0]!=0.0 || pressure[2][1]!=0.0) 
  {
    deviatoric=true;
    volume0=box_start.determinant();
  }

// initialize first number of random pairs for BAOAB-limit case integrator
  if(baoab) for(int i=0;i<3;i++)for(int j=0;j<3;j++) random_pairs[0][i][j]=random.Gaussian();

// here is the main md loop
// Langevin thermostat is applied before and after a velocity-Verlet integrator
// the overall structure is:
//   thermostat
//   update velocities
//   barostat
//   update positions
//   (eventually recompute neighbour list)
//   compute forces
//   update velocities
//   thermostat
//   (eventually dump output informations)
  for(int istep=0;istep<nstep;istep++){

// estimated isothermal compressibility
    double betaT=0.3;
    double prefac_det=-betaT/(3.0*taup);

    int ndeg=3*natoms;
    if(svr) ndeg-=3;

    double dlambda_save=0.0;      // save lambda increment to compute drift
    double b_save;                // save b scalar to compute drift
    Tensor deps_matrix_save;      // save argument of exp matrix to compute drift

    if(svr) gthermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);
    else thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);

// scaling matrices used in barostat implementation
    auto rscaling=Tensor::identity();
    auto vscaling=Tensor::identity();

    for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
      velocities[iatom][k]+=forces[iatom][k]*0.5*tstep/masses[iatom];


// Barostat: isotropic implementation
    if(taup>0.0 && istep%nstbaro==0 && barotype=='I')
    {
      double volume=box.determinant();
      double lambda=sqrt(volume);
      double lambda_D=0.25*temperature*betaT/taup;
      double pint_scalar=(virial[0][0]+virial[1][1]+virial[2][2])/(3*volume);
      if(baroscalev) {
        compute_engkin(natoms,masses,velocities,engkin);
        pint_scalar+=2.0/3.0*(engkin[0][0]+engkin[1][1]+engkin[2][2])/volume;
      } else {
        pint_scalar+=(ndeg/3.0)*temperature/volume;
      }
      double lambda_f=-2*lambda*(pressure_hydro-pint_scalar-0.5*temperature/volume);
      
      double lambda_random;
      if(baoab) {
        random_pairs[1][0][0]=random.Gaussian();
        lambda_random=0.5*(random_pairs[0][0][0] + random_pairs[1][0][0]);
        random_pairs[0][0][0]=random_pairs[1][0][0];
      } else lambda_random=random.Gaussian();

      double dlambda=lambda_D/temperature*lambda_f*nstbaro*tstep + sqrt(2.0*lambda_D*nstbaro*tstep)*lambda_random;
      dlambda_save=dlambda;

      engint+=0.5*dlambda*lambda_f - 0.25*lambda_D/temperature*nstbaro*tstep*lambda_f*lambda_f + 0.5*temperature*log(volume);
      if(!baroscalev) engint+=ndeg/3.0*temperature*log(volume);

      double depsilon=2*log1p(dlambda/sqrt(volume));
      double scaling=exp(depsilon/3.0);
      rscaling*=scaling;
      if(baroscalev) vscaling/=scaling; 
    }

// Barostat: anisotropic implementation, Euler integrator
    else if(taup>0.0 && istep%nstbaro==0 && barotype=='A')
    {
      double volume=determinant(box);
      double D_box=-prefac_det*temperature/volume;

      Tensor pint=virial/volume;
      if(baroscalev) {
        pint+=2*engkin/volume;
      } else {
        for(int i=0;i<3;i++) pint[i][i]+=(ndeg/3.0)*temperature/volume;
      }

      Tensor box_random;
      if(baoab) {
        for(int i=0;i<3;i++) { 
          for(int j=0;j<3;j++) {
            random_pairs[1][i][j]=random.Gaussian();
            box_random[i][j]=0.5*(random_pairs[0][i][j] + random_pairs[1][i][j]);
            random_pairs[0][i][j]=random_pairs[1][i][j];
          }
        }
      } else for(int i=0;i<3;i++) for(int j=0;j<3;j++) box_random[i][j]=random.Gaussian();

      rscaling+=prefac_det*(pext-pint-temperature/volume*Tensor::identity())*nstbaro*tstep + sqrt(2.0*D_box*nstbaro*tstep)*box_random;

      // deviatoric correction
      if(deviatoric) {     
        Tensor stress_dev=pressure-pressure_hydro*Tensor::identity();
        Tensor sigma=volume0*matmul(transpose(box_start_inv),matmul(stress_dev,box_start_inv));
        Tensor pressure_dev=matmul(transpose(box),matmul(0.5*(sigma+transpose(sigma)),box))/volume;
	rscaling+=prefac_det*pressure_dev*nstbaro*tstep;
      }

      // eliminate rotations
      if(!rotations) { 
        rotate_back_box(rscaling);
      }
      if(baroscalev) vscaling=inverse(rscaling);
    }
    
// Barostat: anisotropic implementation, time-reversible integrator 
    else if(taup>0.0 && istep%nstbaro==0 && barotype=='T')
    {
// 1) propagation of lambda for half timestep
      double volume=box.determinant(); 
      double lambda=sqrt(volume);
      double lambda_D=0.25*temperature*betaT/taup;
      double pint_scalar=(virial[0][0]+virial[1][1]+virial[2][2])/(3*volume);
      if(baroscalev) {
        compute_engkin(natoms,masses,velocities,engkin);
        pint_scalar+=2.0/3.0*(engkin[0][0]+engkin[1][1]+engkin[2][2])/volume;
      } else {
        pint_scalar+=(ndeg/3.0)*temperature/volume;
      }
      double lambda_f=-2*lambda*(pressure_hydro-pint_scalar-0.5*temperature/volume);
      // deviatoric correction
      Tensor pressure_dev;
      if(deviatoric) {     
        Tensor stress_dev=pressure-pressure_hydro*Tensor::identity();
        Tensor sigma=volume0*matmul(transpose(box_start_inv),matmul(stress_dev,box_start_inv));
        pressure_dev=+matmul(transpose(box),matmul(0.5*(sigma+transpose(sigma)),box))/volume;
        lambda_f-=2*lambda*(pressure_dev[0][0]+pressure_dev[1][1]+pressure_dev[2][2])/3.0;
      }
      double lambda_random=random.Gaussian();

      double dlambda=lambda_D/temperature*lambda_f*nstbaro*tstep/2.0 + sqrt(lambda_D*nstbaro*tstep)*lambda_random;
      dlambda_save=dlambda;

// 2) propagation of h for a full timestep, at fixed volume
      double volumenew=(lambda+dlambda)*(lambda+dlambda); 
      double b=sqrt(-2.0*prefac_det*temperature/volumenew);

      Tensor pint=virial/volume;
      if(baroscalev) {
        pint+=2*engkin/volume;
      } else {
        for(int d=0;d<3;d++) pint[d][d]+=(ndeg/3.0)*temperature/volume;
      }

      Tensor A=prefac_det*(pext-pint-temperature/volume*Tensor::identity());      
      // deviatoric correction
      if(deviatoric) {     
        A+=prefac_det*pressure_dev;
      }
      Tensor A2=(A-(A[0][0]+A[1][1]+A[2][2])/3.0*Tensor::identity());    
      Tensor box_random;
      for(int i=0;i<3;i++)for(int j=0;j<3;j++) box_random[i][j]=random.Gaussian();
      Tensor box_random2=box_random-(box_random[0][0]+box_random[1][1]+box_random[2][2])/3.0*Tensor::identity();

      Tensor deps_matrix=A2*nstbaro*tstep + b*box_random2*sqrt(nstbaro*tstep);
      rscaling=compute_matrixexp(deps_matrix);
      if(!rotations) { 
        rotate_back_box(rscaling);
      }
      for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
        engint+=-0.5*temperature*nstbaro*tstep/(b*b)*A2[i][j]*A2[i][j] + temperature/(b*b)*deps_matrix[i][j]*A2[i][j];
      }

// 3) propagation of lambda for half timestep
      lambda_random=random.Gaussian();
      dlambda+=lambda_D/temperature*lambda_f*nstbaro*tstep/2.0 + sqrt(lambda_D*nstbaro*tstep)*lambda_random;

      dlambda_save=dlambda;
      b_save=b;
      deps_matrix_save=+deps_matrix;

      engint+=0.5*dlambda*lambda_f - 0.25*lambda_D/temperature*nstbaro*tstep*lambda_f*lambda_f + 0.5*temperature*log(volume);
      if(!baroscalev) engint+=ndeg/3.0*temperature*log(volume);

      rscaling*=pow((lambda+dlambda)*(lambda+dlambda)/volume,1/3.0);
      if(baroscalev) vscaling=inverse(rscaling);
    }

// propagate positions
    for(int iatom=0;iatom<natoms;iatom++){
      positions[iatom]=matmul(positions[iatom],rscaling) + (matmul(vscaling,velocities[iatom]) + matmul(velocities[iatom],rscaling))*tstep/2.0;
    }
    box=matmul(box,rscaling);
    for(int iatom=0;iatom<natoms;iatom++) velocities[iatom]=matmul(vscaling,velocities[iatom]);

// a check is performed to decide whether to recalculate the neighbour list
    check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute_list);
    if(recompute_list){
      compute_list(natoms,listsize,positions,box,listcutoff,point,list);
      for(int iatom=0;iatom<natoms;++iatom) for(int k=0;k<3;++k) positions0[iatom][k]=positions[iatom][k];
      fprintf(stdout,"Neighbour list recomputed at step %d\n",istep);
      fprintf(stdout,"List size: %d\n",point[natoms-1]);
    }

    compute_forces(natoms,listsize,positions,box,forcecutoff,point,list,forces,virial,engconf);
    engconf+=pressure_hydro*box.determinant();

// Barostat: isotropic implementation, compute energy drift
    if(taup>0.0 && istep%nstbaro==0 && (barotype=='I' || barotype=='T'))
    {
      double volume=box.determinant();
      double lambda=sqrt(volume);
      double lambda_D=0.25*temperature*betaT/taup;
      double pint_scalar=(virial[0][0]+virial[1][1]+virial[2][2])/(3*volume);
      if(baroscalev) {
        compute_engkin(natoms,masses,velocities,engkin);
        pint_scalar+=2.0/3.0*(engkin[0][0]+engkin[1][1]+engkin[2][2])/volume;
      } else {
        pint_scalar+=(ndeg/3.0)*temperature/volume;
      }
      double lambda_f=-2*lambda*(pressure_hydro-pint_scalar-0.5*temperature/volume);
      // deviatoric correction
      Tensor pressure_dev;
      if(deviatoric) {     
        Tensor stress_dev=pressure-pressure_hydro*Tensor::identity();
        Tensor sigma=volume0*matmul(transpose(box_start_inv),matmul(stress_dev,box_start_inv));
        pressure_dev=+matmul(transpose(box),matmul(0.5*(sigma+transpose(sigma)),box))/volume;
        lambda_f-=2*lambda*(pressure_dev[0][0]+pressure_dev[1][1]+pressure_dev[2][2])/3.0;
      }

      engint+=0.5*dlambda_save*lambda_f + 0.25*lambda_D/temperature*nstbaro*tstep*lambda_f*lambda_f - 0.5*temperature*log(volume);
      if(!baroscalev) engint-=ndeg/3.0*temperature*log(volume);

// Barostat: anisotropic implementation (time-reversible), compute energy drift
      if (barotype=='T')
      {
        Tensor pint=virial/volume;
        if(baroscalev) {
          pint+=2*engkin/volume;
        } else {
          for(int d=0;d<3;d++) pint[d][d]+=(ndeg/3.0)*temperature/volume;
        }
        Tensor A=prefac_det*(pext-pint-temperature/volume*Tensor::identity());      
        // deviatoric correction
        if(deviatoric) {     
          A+=prefac_det*pressure_dev;
        }
        Tensor A2=(A-(A[0][0]+A[1][1]+A[2][2])/3.0*Tensor::identity());
        for(int i=0;i<3;i++)for(int j=0;j<3;j++) {
          engint+=0.5*temperature*nstbaro*tstep/(b_save*b_save)*A2[i][j]*A2[i][j] 
                  + temperature/(b_save*b_save)*deps_matrix_save[i][j]*A2[i][j];
        }
      }
    }

    for(int iatom=0;iatom<natoms;iatom++) for(int k=0;k<3;k++)
      velocities[iatom][k]+=forces[iatom][k]*0.5*tstep/masses[iatom];

    if(svr) gthermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);
    else thermostat(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,random);

    if(istep%1000==0 && svr) reset_com_velocity(natoms,masses,velocities);

// kinetic energy is calculated
    compute_engkin(natoms,masses,velocities,engkin);

// eventually, write positions and statistics
    if((istep+1)%nconfig==0) write_positions(trajfile,natoms,istep+1,tstep,positions,box,wrapatoms);
    if((istep+1)%nstat==0)   write_statistics(statfile,istep+1,tstep,natoms,engkin,engconf,engint,box.determinant(),virial,barotype,box);
  } 

  write_final_positions(outputfile,natoms,nstep,tstep,positions,box,wrapatoms);

// close the statistic file if it was open:
  if(write_statistics_fp) fclose(write_statistics_fp);

  return 0;
}

};

int main(){
  SimpleMD smd;
  return smd.main(stdin,stdout);
}


