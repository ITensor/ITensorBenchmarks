#include "itensor/all.h"

using namespace itensor;

int
main()
    {
    auto N = 100;
    auto Npart = N;

    auto t1 = 1.0;
    auto t2 = 0.2;
    auto U = 1.0;
    auto V1 = 0.5;

    auto quiet = true;

    auto sweeps = Sweeps(6);
    sweeps.maxdim() = 50,100,200,400,800;
    sweeps.mindim() = 10,20;
    sweeps.cutoff() = 1E-12;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,1E-10,0.0,1E-11,0.0;

    println(sweeps);

    //
    // Initialize the site degrees of freedom.
    //
    auto sites = Electron(N);

    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= N; ++i) 
        {
        ampo += U,"Nupdn",i;
        }
    for(int b = 1; b < N; ++b)
        {
        ampo += -t1,"Cdagup",b,"Cup",b+1;
        ampo += -t1,"Cdagup",b+1,"Cup",b;
        ampo += -t1,"Cdagdn",b,"Cdn",b+1;
        ampo += -t1,"Cdagdn",b+1,"Cdn",b;
        ampo += V1,"Ntot",b,"Ntot",b+1;
        }
    for(int b = 1; b < N-1; ++b)
        {
        ampo += -t2,"Cdagup",b,"Cup",b+2;
        ampo += -t2,"Cdagup",b+2,"Cup",b;
        ampo += -t2,"Cdagdn",b,"Cdn",b+2;
        ampo += -t2,"Cdagdn",b+2,"Cdn",b;
        }
    auto H = toMPO(ampo);

    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    int p = Npart;
    for(int i = N; i >= 1; --i) 
        {
        if(p > i)
            {
            println("Doubly occupying site ",i);
            state.set(i,"UpDn");
            p -= 2;
            }
        else
        if(p > 0)
            {
            println("Singly occupying site ",i);
            state.set(i,(i%2==1 ? "Up" : "Dn"));
            p -= 1;
            }
        else
            {
            state.set(i,"Emp");
            }
        }

    auto psi0 = MPS(state);

    Print(totalQN(psi0));

    //
    // Begin the DMRG calculation
    //
    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",quiet});

    //
    // Measure spin densities
    //
    Vector upd(N),dnd(N);
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        upd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Nup",j)*psi(j));
        dnd(j-1) = elt(dag(prime(psi(j),"Site"))*op(sites,"Ndn",j)*psi(j));
        }

    //println("Up Density:");
    //for(int j = 0; j < N; ++j)
    //    printfln("%d %.10f",1+j,upd(j));
    //println();

    //println("Dn Density:");
    //for(int j = 0; j < N; ++j)
    //    printfln("%d %.10f",1+j,dnd(j));
    //println();

    //println("Total Density:");
    //for(int j = 0; j < N; ++j)
    //    printfln("%d %.10f",1+j,(upd(j)+dnd(j)));
    //println();

    //
    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);

    return 0;
    }
