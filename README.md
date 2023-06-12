IFNPlugin 
=========

This code provides an implementation of the interleaved
flavour-neutralisation (IFN) algorithm from

> Flavoured jets with exact anti-kt kinematics and tests of infrared and collinear safety,
> by Fabrizio Caola, Radoslaw Grabarczyk, Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
> https://arxiv.org/abs/2023.NNNNN

The code is structured following the pattern of a FastJet contrib (but
not yet included in fjcontrib).

You will need `fastjet-config` to be in your path (otherwise edit the
Makefile), and then you can 

- build the library: `make -j`
- build the example: `make -j example`
- run the example: `./example < events/pythia8_Zq_vshort.dat`
- check the output is correct: `make -j check`

To learn how to use the library, the [`example.cc`](example.cc) code is
a good place to get started.

Main principles of IFN
----------------------

The IFN algorithm uses a base jet algorithm (anti-kt or C/A) and
generates a clustering sequence that is kinematically equivalent to that
base algorithm. At each clustering step, the algorithm checks whether it
needs to perform flavour "neutralisation" of the particles involved in
the clustering. The neutralisation can occur with other particles in the
event, not involved in the kinematic clustering. 

Code Structure
--------------

Various flavour-related utilities are to be found in [`FlavInfo.hh`](FlavInfo.hh):
```cpp
// given a particle assign a pdg_id to it, either via a FlavInfo 
// (stores a static flavour)
PseudoJet particle = ...;
particle.set_user_info(new FlavInfo(pdg_id));
// or via a FlavHistory (can track the evolution of the flavour)
// particle.set_user_info(new FlavHistory(pdg_id));

// Retrieve the amount of flavour of a given kind, e.g. amount of b-flavour
// (example given for a particle with a FlavHistory)
FlavInfo flav_info = FlavHistory::current_flavour_of(particle);
// net amount of b-flavour
int nb = flav_info[5];
cout << "nb = " << nb << ", full description: " << flav_info.description() << endl;
```

The main interface to the IFN algorithm is in [`IFNPlugin.hh`](IFNPlugin.hh):
```cpp
// create a base jet definition
JetDefinition base_jet_def(antikt_algorithm, R);
// create an IFNPlugin based on the jet definition, with 
// a given alpha (recommended values, either 1 or 2)
double alpha = 2.0;
JetDefinition jet_def(new IFNPlugin(base_jet_def, alpha));
jet_def.delete_plugin_when_unused();
```

The plugin can be used as standard with the FastJet package, e.g.:
```cpp
vector<PseudoJet> particles;
// ... fill the particles ...
auto jets = jet_def(particles);
for (const auto & jet : jets) {
  cout << "jet pt = " << jet.pt() 
       << ", flav = " << FlavHistory::current_flavour_of(jet) << endl;
}
```