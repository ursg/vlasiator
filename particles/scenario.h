#pragma once
#include <vector>
#include <map>
#include "particles.h"
#include "particleparameters.h"
#include "field.h"

// A "Scenario" describes one possible class of simulation setup, and
// Output quantities we are interested in.
// There are multiple different ways in which particles can be used on
// vlasiator data, this one abstracts these ways.

struct Scenario {

	// Fill in the initial particles, given electromagnetic- and velocity fields
  // Further parameters, depending on the scenario, are given by the parameter object
	virtual std::vector<Particle> initial_particles(Field& E, Field& B, Field& V) {return std::vector<Particle>();};

  // Do something when a new input timestep has been opened
  virtual void new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {};

  // Modify or analyze particle behaviour before they are moved by the particle pusher.
  virtual void before_push(std::vector<Particle>& particles, Field& E, Field& B, Field& V) {};

  // Modify or analzye particle behaviour after they are moved by the particle pusher
  virtual void after_push(int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V) {};

  // Analyze the final state and / or write output
  virtual void finalize(std::vector<Particle>& particles, Field& E, Field& B, Field& V) {};

  // Flags specifying which fields are required for this scenario
  bool needV; // Velocity field (it's always available for initialization)

  Scenario() : needV(false) {};
};





// Specific scenarios below here

// Trace a single particle, it's initial position and velocity given in the parameter file
struct singleParticleScenario : Scenario {
  std::vector<Particle> initial_particles(Field& E, Field& B, Field& V);
  void after_push(int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V);

  singleParticleScenario() {needV = false;};
};

// Sample a bunch of particles from a distribution, create them at a given point, then trace them
struct distributionScenario : Scenario {
  std::vector<Particle> initial_particles(Field& E, Field& B, Field& V);
  void new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V);
  void finalize(std::vector<Particle>& particles, Field& E, Field& B, Field& V);

  distributionScenario() {needV = false;};
};

// Inject particles in the tail continuously, see where they precipitate
struct precipitationScenario : Scenario {
  void new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V);
  void after_push(int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V);

  precipitationScenario() {needV = true;};
};

// For interactive usage from analysator: read positions and velocities from stdin, push those particles.
struct analysatorScenario : Scenario {
  std::vector<Particle> initial_particles(Field& E, Field& B, Field& V);
  void new_timestep(int input_file_counter, int step, double time, std::vector<Particle>& particles, Field& E, Field& B, Field& V);
  analysatorScenario() {needV = false;};
};



template<typename T> Scenario* create_scenario() {
  return new T;
}

Scenario* create_scenario(std::string name); 
