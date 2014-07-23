/*
 * Postprocessing particle trajectory analyzer
 * for Vlasiator
 *
 */
#include <mpi.h>
#include <iostream>
#include <random>
#include <string.h>
#include "particles.h"
#include "physconst.h"
#include "readfields.h"
#include "distribution.h"
#include "particleparameters.h"
#include "../readparameters.h"

// "mode" that the particle pusher operates in
enum pusher_mode {
	single = 0,    // Single particle trajectory tracking
	distribution  // Particle distribution created at a given point
};

void help_message() {
	std::cerr << "Usage: particle_post_pusher <filename> <mode> [args]" << std::endl;
	std::cerr << "Where <filename> is a .vlsv file providing field input and <mode> is one of:" << std::endl;
	std::cerr << std::endl;
	std::cerr  << "  single <x> <y> <z> <vx> <vy> <vz>" << std::endl;
	std::cerr  << "  distribution <x> <y> <z> <temperature> <num_particles>" << std::endl;
}

int main(int argc, char** argv) {

	MPI::Init(argc, argv);

	/* Parse commandline and config*/
	Readparameters parameters(argc, argv, MPI_COMM_WORLD);
	ParticleParameters::addParameters();
	parameters.parse();
	ParticleParameters::getParameters();

	/* Read starting fields from specified input file */
	std::string filename_pattern = ParticleParameters::input_filename_pattern;
	vlsvFieldReader<newVlsv::Reader>* input;
	if(ParticleParameters::dt < 0) {
		input = new vlsvBackwardFieldReader<newVlsv::Reader>;
	} else {
		input = new vlsvForwardFieldReader<newVlsv::Reader>;
	}
	input->open(filename_pattern, ParticleParameters::start_time, ParticleParameters::input_file_step);
	input->open_next_timestep(ParticleParameters::start_time);

	/* Init particles */
	std::vector<Particle> particles;
	pusher_mode mode;

	double dt=ParticleParameters::dt;
	double maxtime=ParticleParameters::end_time - ParticleParameters::start_time;
	int maxsteps = maxtime/dt;

	Field V = input->read_velocity_field();

	if(ParticleParameters::mode == "single") {

		Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);
		/* Look up builk velocity in the V-field */
		Vec3d bulk_vel = V(vpos);

		mode = single;
		particles.push_back(Particle(PhysicalConstantsSI::mp, PhysicalConstantsSI::e, vpos, bulk_vel));
	} else if(ParticleParameters::mode == "distribution") {

		mode = distribution;

		std::default_random_engine generator(ParticleParameters::random_seed);
		Distribution* velocity_distribution=ParticleParameters::distribution(generator);

		Vec3d vpos(ParticleParameters::init_x, ParticleParameters::init_y, ParticleParameters::init_z);

		/* Look up builk velocity in the V-field */
		Vec3d bulk_vel = V(vpos);

		for(unsigned int i=0; i< ParticleParameters::num_particles; i++) {
			/* Create a particle with velocity drawn from the given distribution ... */
			Particle p = velocity_distribution->next_particle();
			/* Shift it by the bulk velocity ... */
			p.v += bulk_vel;
			/* And put it in place. */
			p.x=vpos;
			particles.push_back(p);
		}

		delete velocity_distribution;
	} else {
		std::cerr << "Config option particles.mode not set correctly!" << std::endl;
		return 1;
	}

	/* Dump the initial state */
	if(mode == distribution) {
		write_particles(particles, "particles_initial.vlsv");
	}


	double output_step = ParticleParameters::output_file_step;
	double last_output_time = ParticleParameters::start_time;

	std::cerr << "Pushing " << particles.size() << " particles for " << maxsteps << " steps..." << std::endl;
        std::cerr << "[                                                                        ]\x0d[";

	/* Push them around */
	for(int step=0; step<maxsteps; step++) {

		bool newfile;
		/* Load newer fields, if neccessary */
		newfile = input->open_next_timestep(ParticleParameters::start_time + step*dt);

		Interpolated_Field cur_E = input->interpolate_E(ParticleParameters::start_time+step*dt);
		Interpolated_Field cur_B = input->interpolate_B(ParticleParameters::start_time+step*dt);

		#pragma omp parallel for
		for(unsigned int i=0; i< particles.size(); i++) {
			/* Get E- and B-Field at their position */
			Vec3d Eval,Bval;

			Eval = cur_E(particles[i].x);
			Bval = cur_B(particles[i].x);

			/* Push them around */
			particles[i].push(Bval,Eval,dt);
		}

		/* Write output */
		if(ParticleParameters::start_time + step*dt > last_output_time+output_step) {
			char filename_buffer[256];
			int filenum = floor(ParticleParameters::start_time+step*dt);
			snprintf(filename_buffer,256, ParticleParameters::output_filename_pattern.c_str(),filenum);
			write_particles(particles, filename_buffer);

			last_output_time += output_step;
		}

		/* Draw progress bar */
		if((step % (maxsteps/71))==0) {
			std::cerr << "=";
		}
	}

	/* Dump the final state */
	if(mode == distribution) {
		write_particles(particles, "particles_final.vlsv");
	}

	std::cerr << std::endl;

	delete input;
	MPI::Finalize();
	return 0;
}
