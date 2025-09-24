
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
//! \file       CriticalitySteppingAction.cc
//! \author     Adapted for U-235 criticality simulation
//! \date       September 10, 2025
//!
//! \brief      Implementation of the CriticalitySteppingAction class
//!
//! \details    Records production and loss events for k_eff calculation
//!
// -------------------------------------------------------------
#include "CriticalitySteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Neutron.hh"
#include "G4VProcess.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include <cmath> // For std::log10
#include <string> // For std::to_string
// Constructor
CriticalitySteppingAction::CriticalitySteppingAction()
: G4UserSteppingAction()
{
}
// Destructor
CriticalitySteppingAction::~CriticalitySteppingAction()
{
}
// Process each step
void CriticalitySteppingAction::UserSteppingAction(const G4Step* step)
{
G4Track* track = step->GetTrack();
if (track->GetDefinition() != G4Neutron::Definition()) return;
G4StepPoint* postPoint = step->GetPostStepPoint();
G4String process = postPoint->GetProcessDefinedStep()->GetProcessName();
G4double time = postPoint->GetGlobalTime();
G4AnalysisManager* man = G4AnalysisManager::Instance();
G4String preVol = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
G4String postVol = (postPoint->GetPhysicalVolume()) ? postPoint->GetPhysicalVolume()->GetName() : "";
// Flux tally (track length estimator)
if (preVol == "U235Sphere") {
G4double ekin = step->GetPreStepPoint()->GetKineticEnergy();
if (ekin > 0) {
G4double logE = std::log10(ekin / eV);
G4double length = step->GetStepLength() / cm;
man->FillH1(5, logE, length);
}
}
// Fission sites for Shannon entropy
if (process == "nFission") {
G4ThreeVector pos = postPoint->GetPosition() / cm;
G4int tbin = static_cast<G4int>(time / ns / 10.0); // 10 ns bins
if (tbin >= 0 && tbin < 30) {
G4String name = "fiss_t" + std::to_string(tbin);
G4int id = man->GetH3Id(name);
if (id >= 0) {
man->FillH3(id, pos.x(), pos.y(), pos.z(), 1.0);
} else {
G4cout << "Warning: Histogram " << name << " not found." << G4endl;
}
}
}
// Production: neutrons from fission or inelastic
if (process == "nFission" || process == "hadInelastic") {
const G4TrackVector* secondaries = step->GetSecondary();
G4int numNeut = 0;
for (const auto* sec : *secondaries) {
if (sec->GetDefinition() == G4Neutron::Definition()) numNeut++;
}
if (numNeut > 0) {
man->FillH1(2, time / ns, static_cast<G4double>(numNeut)); // prod_hist id=2
}
}
// Loss: absorption (if killed)
if ((process == "nCapture" || process == "nFission" || process == "hadInelastic") && track->GetTrackStatus() == fStopAndKill) {
man->FillH1(3, time / ns, 1.0); // loss_hist id=3
}
// Loss: escape from sphere to world
if (preVol == "U235Sphere" && postVol == "World") {
man->FillH1(3, time / ns, 1.0); // loss_hist
track->SetTrackStatus(fStopAndKill);
}
}
