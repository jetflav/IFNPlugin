// To run this example, use the following command:
//
//   ./example < events/pythia8_Zq_vshort.dat
//
//----------------------------------------------------------------------
// $Id$
//
// Copyright (c) 2023, Fabrizio Caola, Radoslaw Grabarczyk, 
// Maxwell Hutt, Gavin P. Salam, Ludovic Scyboz, and Jesse Thaler
//
//----------------------------------------------------------------------
// This file is part of FastJet contrib.
//
// It is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at
// your option) any later version.
//
// It is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
// License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this code. If not, see <http://www.gnu.org/licenses/>.
//----------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "fastjet/PseudoJet.hh"
#include <sstream>
#include "FlavNeutraliserPlugin.hh" // In external code, this may become fastjet/contrib/FlavNeutraliserPlugin.hh

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

// forward declaration to make things clearer
void read_event(vector<PseudoJet> &event);

//----------------------------------------------------------------------
int main(){

  // set up the IFNPlugin

  int n_events = 10;
  for (int iev = 0; iev < n_events; iev++) {
    //----------------------------------------------------------
    // read in input particles: see that routine for info 
    // on how to set up the PseudoJets with flavour information
    vector<PseudoJet> event;
    read_event(event);
    cout << "# read an event with " << event.size() << " particles" << endl;

    //----------------------------------------------------------
    // illustrate how this FlavNeutraliserPlugin contrib works

  }

  return 0;
}

// read in input particles and set up PseudoJets with flavour information
void read_event(vector<PseudoJet> &event){  
    // read in the input particles and their PDG IDs
    string line;
    double px, py, pz, E;
    int    pdg_id;
    event.resize(0);
    while(getline(cin,line)) {
      if(line[0] == '#') continue;

      istringstream iss(line);
      iss >> px >> py >> pz >> E >> pdg_id;
      // create a fastjet::PseudoJet with these components and put it onto
      // back of the input_particles vector
      PseudoJet p(px,py,pz,E);

      // assign information about flavour (will be deleted automatically)
      // COMMENT: we probably want to introduce a FlavHistory constructor 
      // that takes just a PDG ID. 
      p.set_user_info(new FlavHistory(FlavInfo(pdg_id)));
      event.push_back(p);

      if (cin.peek() == '\n' || cin.peek() == EOF) {
        getline(cin,line);
        break;
      }
    }
}
