/*
The goal of this program is to simulate P and Rayleigh waves muliple scattering due to perturbations on Lame parameter in a homogeneous medium.
The code is based on the Monte Carlo simulation. The idea is to simulate a wave particle (P/Rayleigh) propagating and scattering in a time period multiple times.

Variables and packages in the code:
gsl_rng.h is for generating random number.
Dz is the thickness for each layer in the medium. Although the medium is homogenous, I divide the medium into layers for computing Rayleigh-wave eigenfunctions at each depth.
Nlay is the number of layers.
dt is the time interval in one simulation.
nsteps is the total time steps in one simulation.
nwalks is the total simulation times.
dz is the depth range for detecting waves/particles.
dr is the radius range for detecting waves/particles.
MASTER master thread in MPI.
vp is the P-wave velocity of the medium.
vs is the S-wave velocity
vr is the Rayleigh-wave phase/group velocity of the medium. As the medium is homogenous, there is no dispersion.
n is the scatter density.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <mpi.h>
#define om  (2.0*M_PI)
#define Dz 0.1
#define dt 1E-1
#define nsteps 20000
#define nwalks 10000000
#define dz 5.
#define dr 5.
#define nr 80
#define MASTER 0
#define Nlay 1001
#define vp 5.0
#define vs 3.0
#define vr 2.7425
#define n 2


// sigR cross section of incident Rayleigh waves, equal sigRR+sighRP
// sigeP effective cross section of incidnet P wave
// intpisq integral of the square eigenfunctions in the whole depth.
double sigR,sigRR,sigRP,sigPP,sigeP,zsource,fre,intpisq;
// eigenfunctions and probability distribution of Rayleigh wave scattering to P wave in depth
double eigfun[3][Nlay],pRPz[Nlay-1];
// Energy/number of P (in a shallow slab), Rayleigh and total P waves at each time step.
double Eb[nsteps],Es[nsteps],Ebt[nsteps];
double Eblay[nsteps][Nlay],Eblaysum[nsteps][Nlay];
double  Envb[nsteps][nr],Envs[nsteps][nr],dv[nr],ds[nr];
double Ebsum[nsteps],Essum[nsteps],Ebtsum[nsteps];
double Envbsum[nsteps][nr],Envssum[nsteps][nr];
long long int ndiffusssum[nsteps][nr], ndiffuss[nsteps][nr];

// azi0 is the propagation azimuth
// ang0 is the propagation angle from the vertical downgoing direction.
// x,y,z are the partical location
// kx0, ky0, kz0 is the propagation slowness in x,y,z direction
// deltat is the traveltime before next time step
// freetime is the propagation time between two scatterers.
double azi0,ang0,x,y,z,t,kx0,ky0,kz0,deltat,freetime;
int it=0;
long long int idiffus=0;
// for random number
gsl_rng * r;
// wave/partical mode: P or Rayleigh
char mode;

//***************Rayleigh waves
// eigenfunctions for Rayleigh wave in different depths
void eigenfun(void)
{
 int ilay;
 double depth,kr,vrsq1,vrsq2;
 kr=om*fre/vr;
 vrsq1=vr*vr/vs/vs;
 vrsq2=vr*vr/vp/vp;
 for(ilay=0;ilay<Nlay;ilay++)
 {
    depth=ilay*Dz;
//  r1 horizontal displacement eigenfunction
    eigfun[0][ilay]=-kr*exp(-kr*sqrt(1-vrsq2)*depth)+kr*(2-vrsq1)/2*exp(-kr*sqrt(1-vrsq1)*depth);
//  r2 vertical displacement
    eigfun[1][ilay]=-kr*sqrt(1-vrsq2)*exp(-kr*sqrt(1-vrsq2)*depth)+kr*2*sqrt(1-vrsq2)/(2-vrsq1)*exp(-kr*sqrt(1-vrsq1)*depth);  
// 1/kr*dr2/dz
    eigfun[2][ilay]=kr*(1-vrsq2)*exp(-kr*sqrt(1-vrsq2)*depth)-kr*2*sqrt(1-vrsq2)*sqrt(1-vrsq1)/(2-vrsq1)*exp(-kr*sqrt(1-vrsq1)*depth);
 }
}

// intergral on eigenfunctions
void inteigfun(void)
{
  int ilay;
  intpisq=0;
// intergral on r1^2+r2^2 dz
  for(ilay=0;ilay<Nlay;ilay++)
    intpisq=intpisq+eigfun[0][ilay]*eigfun[0][ilay]+eigfun[1][ilay]*eigfun[1][ilay];
  intpisq=intpisq-eigfun[0][0]*eigfun[0][0]/2-eigfun[0][Nlay-1]*eigfun[0][Nlay-1]/2;
  intpisq=intpisq-eigfun[1][0]*eigfun[1][0]/2-eigfun[1][Nlay-1]*eigfun[1][Nlay-1]/2;
  intpisq=intpisq*Dz;  
}

// Probability for R->P depth distribution
void proRPdep(void)
{
 int ilay;
 double intg=0,intergrand;
// intergral on (r1+1/kr*dr2/dz)^2 dz
  for(ilay=0;ilay<Nlay-1;ilay++){
    intergrand = pow((eigfun[0][ilay]+eigfun[2][ilay]),2)+pow((eigfun[0][ilay+1]+eigfun[2][ilay+1]),2);
    intg=intg+intergrand/2;
    pRPz[ilay]=intg;
  }
  for(ilay=0;ilay<Nlay-1;ilay++)
    pRPz[ilay]=pRPz[ilay]/intg;
}

// find the depth for R->P depth
double RPdepz(double randn)
{
  int ilay;
  double depz=0;
// find the possibility bin
  for(ilay=0;ilay<Nlay-1;ilay++){
    if(ilay==0&&randn<pRPz[0]){
// interpolate the depth
      depz = randn/pRPz[0];
      break;}
    if(randn>=pRPz[ilay]&&randn<pRPz[ilay+1]){
      depz = (randn - pRPz[ilay])/(pRPz[ilay+1]-pRPz[ilay])+ilay+1;
      break;}
  }
  return (depz*Dz);
}
//***************Cross sections for lame parameters
// cross sections for P->Rayleigh
double csPR(double zz)
{   
  int iz;
  double cs,cs1,cs2,constant;
  iz = floor(zz/Dz);
  constant=pow(om*fre/vp,3)/4/pow(vr,4)*pow(vp*vp-2*vs*vs,2)/intpisq;
  if(iz>Nlay-1)
    iz = Nlay-1;
  cs1 = pow((eigfun[0][iz]+eigfun[2][iz]),2);
  if(iz<Nlay-1)
  {
     cs2=pow((eigfun[0][iz+1]+eigfun[2][iz+1]),2);
     cs = cs1 + (cs2-cs1)/Dz*(zz-(double)(iz*Dz));
  }
  else
    cs = cs1;
  cs = cs * constant;
  return cs;
}
// check intergral on sigPR
double intcsPR(void)
{
  double intg=0;
  int iz;
  for(iz=0;iz<Nlay;iz++)
    intg=intg+csPR((double)iz*Dz);
  intg=intg-csPR(0)/2-csPR(Nlay-1)/2;
  intg=intg*Dz;
  return intg;
}
// cross section for P->P
double csPP(void)
{
  return (pow(om*fre/vp,4)/4/M_PI*pow((vp*vp-2*vs*vs)/vp/vp,2));
}
// R->P
double csRP(void)
{
  int ilay;
  double intg=0;
// intergral on (r1+1/kr*dr2/dz)^2 dz
  for(ilay=0;ilay<Nlay;ilay++)
    intg=intg+(eigfun[0][ilay]+eigfun[2][ilay])*(eigfun[0][ilay]+eigfun[2][ilay]);
  intg=intg-(eigfun[0][0]+eigfun[2][0])*(eigfun[0][0]+eigfun[2][0])/2;
  intg=intg-(eigfun[0][Nlay-1]+eigfun[2][Nlay-1])*(eigfun[0][Nlay-1]+eigfun[2][Nlay-1])/2;
  intg=intg*Dz;    
  return  (pow(om*fre/vp,4)/4/M_PI/vp/pow(vr,3)*pow(vp*vp-2*vs*vs,2)/intpisq*intg);
}
// R->R
double csRR(void)
{
  int ilay;
  double intg=0;
// intergral on (r1+1/kr*dr2/dz)^4 dz
  for(ilay=0;ilay<Nlay;ilay++)
    intg=intg+pow(eigfun[0][ilay]+eigfun[2][ilay],4);
  intg=intg-pow(eigfun[0][0]+eigfun[2][0],4)/2;
  intg=intg-pow(eigfun[0][Nlay-1]+eigfun[2][Nlay-1],4)/2;
  intg=intg*Dz;
  return (pow(om*fre,3)/4/pow(vr,7)*pow(vp*vp-2*vs*vs,2)/intpisq/intpisq*intg);
}
// effective cross sections for incident P waves
double cseP(void)
{
   int ilay;
   double maxcs=0;
   for(ilay=0;ilay<Nlay;ilay++)
   {
     if(csPR((double)ilay*Dz)>maxcs)
       maxcs=csPR((double)ilay*Dz);
   }
   return (maxcs+sigPP);
}
// imaginary cross section / tricks in dealing with varied cross sections with depth for incident P waves
double csim(double z)
{
  return (sigeP - (sigPP + csPR(z)));
}

//*************Simulation
// initial surface wave energy
double e0s(double zz)
{
  return (om*fre*om*fre*vr*intpisq);
 }
// initial body wave energy
double e0b(double zz)
{
  return (om*fre*om*fre*vp);
}

// count body wave energy (partical number) 
void detectb(void)
{
  double r;
  int ir,iz;
  r = sqrt(x*x +y*y);
  ir = (int)(r/dr);
  iz = (int)(z/Dz);
  if ((ir < nr) && (it < nsteps) && (z < dz)){
    Envb[it][ir] +=1.0;
  }
// P wave energy in a shallow slab
  if ( (z < 20*Dz) && (it < nsteps))
    Eb[it] += 1.0/Dz/20;
// P wave enery in every layer
  if ( (iz < Nlay) && (it < nsteps))
    Eblay[it][iz] += 1.0/Dz;
// total P wave energy
  if (it < nsteps)
    Ebt[it] +=1.0;
}

// Rayleigh wave energy (partical number) 
void detects(void)
{
  double r;
  int ir;
  r = sqrt(x*x +y*y);
  ir = (int)(r/dr);
// Rayleigh wave energy in different offsets
  if ((ir < nr) && (it < nsteps)){
    Envs[it][ir] +=1.0;
    ndiffuss[it][ir] +=idiffus;
  }
// total Rayleigh wave energy
  if (it < nsteps)
    Es[it] +=1.0;
}

// initalize all settings, e.g. propagation direction (kx/y/z) and free time
void myinit(void)
{
  double pb, ran,rancos;
  t = 0.00000001;
  x=0.0; y=0.0; z = zsource;
  idiffus=0;
  pb = e0b(z)/(e0b(z) +e0s(z));
// generate a P/Rayleigh partical (wave)
  ran = gsl_rng_uniform_pos(r);
  if (ran<0){
    mode ='b';
// dipping angle 0-pi
    rancos = -1 + 2.0*gsl_rng_uniform_pos(r);
    ang0 = acos(rancos);
// azimuth 0-2pi 
    azi0 = 2.0*M_PI* gsl_rng_uniform_pos(r) ; 
    kx0 = sin(ang0)*cos(azi0);
    ky0 = sin(ang0)*sin(azi0); 
    kz0 = cos(ang0); 
// effective P-wave cross section
   freetime = -log(gsl_rng_uniform_pos(r))/sigeP/n/vp;
    detectb();
  }
  else{
    mode ='s';
    azi0 = 2.0*M_PI* gsl_rng_uniform_pos(r) ;
    kx0 = cos(azi0);
    ky0 = sin(azi0);
    kz0 = 0.0;
    freetime = -log(gsl_rng_uniform_pos(r))/sigR/n/vr;
    detects();
  }
 }

// P-wave scattering simulation
void scatterb(void)
{
  double pb,pbb,randn,rancos;
  int iz;
  iz = round(z/Dz);
// in case P waves propagte too deep
  if(iz>=Nlay)
    iz=Nlay-1;
// generate random number
  randn = gsl_rng_uniform(r);
// possibility for real scattering
  pb = (sigPP + csPR(z))/sigeP;
//****P->P/Rayleigh
  if (randn < pb){ 
    /*  Real scattering occurs */
    idiffus ++;
    pbb = sigPP/(sigPP+csPR(z));
    randn = gsl_rng_uniform(r);
// P-> P
    if (randn < pbb){
      mode = 'b';
      rancos = -1 + 2.0*gsl_rng_uniform_pos(r);
      ang0 = acos(rancos);      
    //  ang0 = M_PI*gsl_rng_uniform_pos(r) ; 
    // azimuth 0-2pi 
      azi0 = 2.0*M_PI* gsl_rng_uniform_pos(r) ; 
      kx0 = sin(ang0)*cos(azi0);
      ky0 = sin(ang0)*sin(azi0); 
      kz0 = cos(ang0); 
      freetime = -log(gsl_rng_uniform_pos(r))/sigeP/n/vp;
      return;
    }
//P->S
    else{
      mode = 's';
      azi0 = 2.0*M_PI* gsl_rng_uniform_pos(r) ;
      kx0 = cos(azi0);
      ky0 = sin(azi0);
      kz0 = 0.0;
      freetime = -log(gsl_rng_uniform_pos(r))/sigR/n/vr;
      return;
    }
  }
// // just a stop in the propagation path and will not change the propagation direction
  else { 
    /* Imaginary scattering occurs   */
    mode = 'b';
    freetime = -log(gsl_rng_uniform_pos(r))/sigeP/n/vp;
    return;
  }
}

// Surface/Rayleigh-wave scattering simulation
void scatters(void)
{
  double randn,rancos;
  randn = gsl_rng_uniform(r);
  idiffus++;
  if (randn < sigRR/sigR){ 
    mode='s';
    azi0=2.0*M_PI* gsl_rng_uniform_pos(r) ;
    kx0 = cos(azi0);
    ky0 = sin(azi0);
    kz0 = 0.0;
    freetime = -log(gsl_rng_uniform_pos(r))/sigR/n/vr;
    return;
  }
  else{
    mode = 'b';
// test how depth distribution changes affect 
//    z = Dz*(Nlay-1)*gsl_rng_uniform_pos(r);
    z = RPdepz(gsl_rng_uniform_pos(r));
//    ang0 = M_PI*gsl_rng_uniform_pos(r) ; 
    rancos = -1 + 2.0*gsl_rng_uniform_pos(r);
    ang0 = acos(rancos);
    azi0 = 2.0*M_PI* gsl_rng_uniform_pos(r) ; 
    kx0 = sin(ang0)*cos(azi0);
    ky0 = sin(ang0)*sin(azi0); 
    kz0 = cos(ang0); 
    freetime = -log(gsl_rng_uniform_pos(r))/sigeP/n/vp;
    return;
  }
}

// Body (P) wave propagation
void propagb(void)
{
  double dl,deltal;
  /******** Body wave case**********/
// no scattering
  if (freetime > deltat){
      dl = vp*deltat;
      z +=dl*kz0;
      x +=dl*kx0;
      y +=dl*ky0;
      t +=deltat;
      it++;
// spectral free surface boundary condition
      if (z<0.0){
	z=-z;
	kz0=-kz0;
      }
      detectb();
      freetime -=deltat;
      deltat = dt;
      return;
    }
// scattering
    else {
// before scattering, the former propagation direction and remaining distance
      deltal = vp*freetime;
      z +=deltal*kz0;
      x += deltal*kx0; 
      y += deltal*ky0; 
// encounter free surface
      if (z<0.0){
      z = -z;
      kz0 = - kz0;
      }
      t +=freetime;
      deltat -=freetime;  
// scattering 
      scatterb();
    }
}

// Surface-wave propagation
void propags(void)
{ 
  double dl,deltal;
  if (freetime > deltat){
    dl = vr*deltat;
    x +=dl*kx0;
    y +=dl*ky0;
    t +=deltat;
    it++;
    detects();
    freetime -=deltat;
    deltat = dt;
    return;
  }
// scattering
  else {
    deltal = vr*freetime;
    x +=deltal*kx0;
    y +=deltal*ky0;
    t +=freetime;
    deltat -=freetime;   
    scatters();
  }
}


int main(int argc,char *argv[])
{
  int istep,iz,ir,imft,npart;
  int taskid,numtasks,rc;
  FILE *fequip,*fenvb,*fenvs,*fparam,*fndiffuss,*fenvbz;
  gsl_rng_env_setup();
  r = gsl_rng_alloc(gsl_rng_gfsr4);

// MPI initilization
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  printf ("MPI task %d has started...\n", taskid);
  printf("taskid %d \n",taskid);

  gsl_rng_set(r,taskid);

// array initilization
  for(istep=0;istep<nsteps;istep++){
    Eb[istep]=0.0; Es[istep]=0.0;
    Ebt[istep] = 0.0;
    for(ir=0;ir<nr;ir++){
      Envb[istep][ir]=0.0;
      Envs[istep][ir]=0.0;
      ndiffuss[istep][ir] =0;
    }
    for(iz=0;iz<Nlay;iz++)
      Eblay[istep][iz] = 0;
  }
    
  for(ir=0;ir<nr;ir++){
    dv[ir] = dz*M_PI*dr*(2.0*ir+1.0)*dr;
    /* dv[jj] = dz*M_PI*dr*dr;*/
    ds[ir] = dv[ir]/dz; 
  }
  
// frequency
  fre=1;
// prepare eigenfunctions and cross sections
  eigenfun();
  inteigfun();
  sigPP=csPP();
  sigRR=csRR();
  sigRP=csRP();
  sigeP=cseP();
  sigR=sigRR+sigRP;
  proRPdep();
// source depth
  zsource = 1;
   /* exit(0);  */
  if (taskid == MASTER){
// output parameters
    printf("Checking part:");
    printf("zsource %e km I1 %e\n",zsource,intpisq);
    printf("Rayleigh-wave mean free time %le s\n",1/(sigRR+sigRP)/n/vr);
    printf("Theoretical equipartion %le\n",intcsPR()*vp/sigRP/vr);
    fparam=fopen("params3.dat","w");
    fprintf(fparam,"csPP %e csPR(0) %e csRP %e csRR %e cseP %e \n",sigPP,csPR(0),sigRP,sigRR,sigeP);
    fprintf(fparam,"Dz zsource nr nsteps %le %le %d %d \n",Dz,zsource,nr,nsteps);
    fprintf(fparam,"source %e \n",zsource);
    
    for(iz=0;iz<Nlay-1;iz++){
      fprintf(fparam,"%le %le %le %le \n",iz*Dz,csPR((double)iz*Dz),csim((double)iz*Dz),pRPz[iz]); 
    }
    fclose(fparam);
  }
 
// ************************************
// Simulation start
  for (npart=1;npart<=(nwalks/numtasks);npart++){
   if (taskid == MASTER){
// output every 1/10 of the whole simulations
     if ( (npart % (nwalks/numtasks/10)) == 0)
       printf("npart %d \n",npart);}
// Simulation initialization
   it =0;
   myinit();
   deltat= dt;
   do {
     if (mode == 'b')
       propagb();
     else
       propags();
   } while (it < nsteps);
  }
// Collect simulation results from all processors
  rc = MPI_Reduce(&Eb, &Ebsum,nsteps, MPI_DOUBLE, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);
  rc = MPI_Reduce(&Eblay, &Eblaysum,nsteps*Nlay, MPI_DOUBLE, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);
  rc = MPI_Reduce(&Es, &Essum,nsteps, MPI_DOUBLE, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);
  rc = MPI_Reduce(&Ebt, &Ebtsum,nsteps, MPI_DOUBLE, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);
  rc = MPI_Reduce(&Envb, &Envbsum,nsteps*nr, MPI_DOUBLE, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);
  rc = MPI_Reduce(&Envs, &Envssum,nsteps*nr, MPI_DOUBLE, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);

  rc = MPI_Reduce(&ndiffuss, &ndiffusssum,nsteps*nr, MPI_LONG_LONG_INT, MPI_SUM,
		  MASTER, MPI_COMM_WORLD);
  if (rc != MPI_SUCCESS)
    printf("%d: failure on mpc_reduce\n", taskid);

// ************************************
// Output  
  /* exit(0); */
  if (taskid == MASTER) {
    fequip = fopen("equip3.dat","w");
    fenvb = fopen("envb3.dat","w");
    fenvs = fopen("envs3.dat","w");
    fndiffuss = fopen("ndiffuss3.dat","w");
    fenvbz = fopen("envPdepz.dat","w");
    imft=0;
    for (istep=0;istep<nsteps;istep++){
      fprintf(fequip,"%e %e %e %e %e %e \n",istep*dt,Ebsum[istep]/nwalks,Essum[istep]/nwalks,Ebtsum[istep]/nwalks,Essum[istep]/Ebsum[istep],Essum[istep]/Ebtsum[istep]);
      for (ir=0;ir<nr;ir++){
	fprintf(fenvb,"%e ",Envbsum[istep][ir]/dv[ir]/nwalks);
	fprintf(fenvs,"%e ",Envssum[istep][ir]/ds[ir]/nwalks);
	if (Envssum[istep][ir] != 0.)
	  fprintf(fndiffuss,"%e ",1.*ndiffusssum[istep][ir]/Envssum[istep][ir]);
	else
	  fprintf(fndiffuss,"%e ",0.);
      }
      fprintf(fenvb,"\n");
      fprintf(fenvs,"\n");
      fprintf(fndiffuss,"\n");
      if((int)(istep*dt*(sigRR+sigRP)*n*vr)==imft){
      if(((int)(istep*dt*(sigRR+sigRP)*n*vr)) %10 == 0)
      {
        for (iz=0;iz<Nlay;iz++){
	  fprintf(fenvbz,"%e ",Eblaysum[istep][iz]/nwalks);
        }
        fprintf(fenvbz,"\n");
      }
       imft++;}
    }
    fclose(fequip);
    fclose(fenvb);
    fclose(fenvs);
    fclose(fenvbz);
  }
  gsl_rng_free(r);
  MPI_Finalize();
  return 0;
}

