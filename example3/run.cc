#include <cstdlib>
#include <vector>
#include <chrono>
#include <cmath>

#include <itensor/all.h>

#include "electronk.h"

using namespace itensor;

int
main()
  {
  // Parse command line arguments
  //double T=1;
  double U=2;
  int nx=4;
  int ny=4;
  int seed = 41;
  double t=1;

  // Input checks
  //assert(!close(T, 0));
  //double beta = 1. / T;
  itensor::seedRNG(seed);

  auto args = Args("Kmod", ny);
  args.add("ConserveQNs", true);
  int n_sites=nx*ny;

  auto sweeps = Sweeps(5);
  sweeps.maxdim() = 10000;
  sweeps.cutoff() = 1e-8;
  sweeps.noise() = 1e-8,0;


  /////////////////////////////////
  // K-space Hamiltonian

  SiteSet sitesk = ElectronK(n_sites, args);
  // SiteSet sitesk = Electron(n_sites);

  auto Ampo = AutoMPO(sitesk);

  // "hopping" in y-direction
  for (int x=0; x<nx; ++x)
  for (int ky=0; ky<ny; ++ky)
    {
    int s = x*ny + ky + 1; // itensor is 1-indexed
    double disp = -2*t*cos( (2 * M_PI / ny) * ky );
    if (std::abs(disp)>1e-12)
      {
      Ampo += disp,"Nup",s; 
      Ampo += disp,"Ndn",s;
      }
    }

  // hopping in x-direction
  for (int x=0; x<nx-1; ++x)
  for (int ky=0; ky<ny; ++ky)
    {
    int s1 = x*ny + ky + 1; // itensor is 1-indexed
    int s2 = (x+1)*ny + ky + 1;

    Ampo += -t,"Cdagup",s1,"Cup",s2; 
    Ampo += -t,"Cdagup",s2,"Cup",s1;
    Ampo += -t,"Cdagdn",s1,"Cdn",s2;
    Ampo += -t,"Cdagdn",s2,"Cdn",s1;
    }

  // Hubbard interaction
  for (int x=0; x<nx; ++x)
  for (int ky=0; ky<ny; ++ky)
    {
    for (int py=0; py<ny; ++py)
    for (int qy=0; qy<ny; ++qy)
      {
      int s1 = x*ny + (ky+qy+ny)%ny + 1;
      int s2 = x*ny + (py-qy+ny)%ny + 1;
      int s3 = x*ny + py + 1;
      int s4 = x*ny + ky + 1;
      Ampo += (U/ny),"Cdagdn",s1,"Cdagup",s2,"Cup",s3,"Cdn",s4;
      // Ampo += (U/ny),"Cdn",s4,"Cup",s3,"Cdagup",s2,"Cdagdn",s1;
      }
    }
  auto H = toMPO(Ampo);
  Print(maxLinkDim(H));

  // Create start state
  auto state = InitState(sitesk);
  for (auto i : range1(n_sites))
    {
    int x = (i-1)/ny;
    int y = (i-1)%ny;

    if(x%2==0)
      {
      if (y%2==0) state.set(i,"Up");
      else        state.set(i,"Dn");
      }
    else
      {
      if (y%2==0) state.set(i,"Dn");
      else        state.set(i,"Up");
      }

    }

  state.set(1,"Emp");
  state.set(ny+1,"Emp");
  auto psi0 = MPS(state);
  
  // auto psi0 = InitialState(sitesk, n_sites, true, seed);

  // auto psi0 = MPS(sitesk);
  psi0.position(1);
  Print(norm(psi0));
  Print(inner(psi0, H, psi0));

  auto [energyk, psi] = dmrg(H, psi0, sweeps, {"Quiet", true});
  Print(energyk);

  return 0;
  }
