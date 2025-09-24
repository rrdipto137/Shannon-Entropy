//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications, and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// -------------------------------------------------------------
//!
//! \file       CriticalitySteppingAction.hh
//! \author     Adapted for U-235 criticality simulation
//! \date       September 10, 2025
//!
//! \brief      Definition of the CriticalitySteppingAction class
//!
//! \details    Records production and loss events for criticality parameters
//!
// -------------------------------------------------------------

#ifndef CRITICALITYSTEPPINGACTION
#define CRITICALITYSTEPPINGACTION

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class CriticalitySteppingAction : public G4UserSteppingAction
{
public:
    // Constructor
    CriticalitySteppingAction();

    // Destructor
    virtual ~CriticalitySteppingAction();

    // Method to process each step
    virtual void UserSteppingAction(const G4Step* step);
};

#endif // CRITICALITYSTEPPINGACTION
