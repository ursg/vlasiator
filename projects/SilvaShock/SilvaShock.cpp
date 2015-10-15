/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "../../common.h"
#include "../../readparameters.h"
#include "../../backgroundfield/backgroundfield.h"

#include "SilvaShock.h"
#include "noise.h"

namespace projects {
   SilvaShock::SilvaShock(): TriAxisSearch() { }
   SilvaShock::~SilvaShock() { }



   bool SilvaShock::initialize(void) {return true;}

   void SilvaShock::addParameters() {
      typedef Readparameters RP;
      RP::add("SilvaShock.BX0u", "Upstream mag. field value (T)", 1.0e-9);
      RP::add("SilvaShock.BY0u", "Upstream mag. field value (T)", 2.0e-9);
      RP::add("SilvaShock.BZ0u", "Upstream mag. field value (T)", 3.0e-9);
      RP::add("SilvaShock.EX0u", "Upstream elec. field value (V/m)", 0.0);
      RP::add("SilvaShock.EY0u", "Upstream elec. field value (V/m)", 0.0);
      RP::add("SilvaShock.EZ0u", "Upstream elec. field value (V/m)", 0.0);
      RP::add("SilvaShock.VX0u", "Upstream Bulk velocity in x", 0.0);
      RP::add("SilvaShock.VY0u", "Upstream Bulk velocity in y", 0.0);
      RP::add("SilvaShock.VZ0u", "Upstream Bulk velocuty in z", 0.0);
      RP::add("SilvaShock.rhou", "Upstream Number density (m^-3)", 1.0e7);
      RP::add("SilvaShock.Temperatureu", "Upstream Temperature (K)", 2.0e6);

      RP::add("SilvaShock.BX0d", "Downstream mag. field value (T)", 1.0e-9);
      RP::add("SilvaShock.BY0d", "Downstream mag. field value (T)", 2.0e-9);
      RP::add("SilvaShock.BZ0d", "Downstream mag. field value (T)", 3.0e-9);
      RP::add("SilvaShock.EX0d", "Downstream elec. field value (V/m)", 0.0);
      RP::add("SilvaShock.EY0d", "Downstream elec. field value (V/m)", 0.0);
      RP::add("SilvaShock.EZ0d", "Downstream elec. field value (V/m)", 0.0);
      RP::add("SilvaShock.VX0d", "Downstream Bulk velocity in x", 0.0);
      RP::add("SilvaShock.VY0d", "Downstream Bulk velocity in y", 0.0);
      RP::add("SilvaShock.VZ0d", "Downstream Bulk velocuty in z", 0.0);
      RP::add("SilvaShock.rhod", "Downstream Number density (m^-3)", 1.0e7);
      RP::add("SilvaShock.Temperatured", "Downstream Temperature (K)", 2.0e6);

      RP::add("SilvaShock.nSpaceSamples", "Number of sampling points per spatial dimension", 2);
      RP::add("SilvaShock.nVelocitySamples", "Number of sampling points per velocity dimension", 5);
      RP::add("SilvaShock.maxwCutoff", "Cutoff for the maxwellian distribution", 1e-12);
      RP::add("SilvaShock.Scale_x", "Scale length in x (m)", 2.0e6);
      RP::add("SilvaShock.Scale_y", "Scale length in y (m)", 2.0e6);
      RP::add("SilvaShock.Width", "Shock Width (m)", 50000);

      RP::add("SilvaShock.NoiseStrength", "Strength of Perlin Noise on the magnetic field (T)", 1e-9);
      RP::add("SilvaShock.NoiseScale", "Largest scale Perlin Noise on the magnetic field (m)", 100000);
      RP::add("SilvaShock.NoiseOctaves", "Number of octaves Perlin Noise on the magnetic field (m)", 10);
   }

   void SilvaShock::getParameters() {
      typedef Readparameters RP;
      RP::get("SilvaShock.BX0u", BX0u);
      RP::get("SilvaShock.BY0u", BY0u);
      RP::get("SilvaShock.BZ0u", BZ0u);
      RP::get("SilvaShock.EX0u", EX0u);
      RP::get("SilvaShock.EY0u", EY0u);
      RP::get("SilvaShock.EZ0u", EZ0u);
      RP::get("SilvaShock.VX0u", VX0u);
      RP::get("SilvaShock.VY0u", VY0u);
      RP::get("SilvaShock.VZ0u", VZ0u);
      RP::get("SilvaShock.rhou", DENSITYu);
      RP::get("SilvaShock.Temperatureu", TEMPERATUREu);

      RP::get("SilvaShock.BX0d", BX0d);
      RP::get("SilvaShock.BY0d", BY0d);
      RP::get("SilvaShock.BZ0d", BZ0d);
      RP::get("SilvaShock.EX0d", EX0d);
      RP::get("SilvaShock.EY0d", EY0d);
      RP::get("SilvaShock.EZ0d", EZ0d);
      RP::get("SilvaShock.VX0d", VX0d);
      RP::get("SilvaShock.VY0d", VY0d);
      RP::get("SilvaShock.VZ0d", VZ0d);
      RP::get("SilvaShock.rhod", DENSITYd);
      RP::get("SilvaShock.Temperatured", TEMPERATUREd);

      RP::get("SilvaShock.nSpaceSamples", nSpaceSamples);
      RP::get("SilvaShock.nVelocitySamples", nVelocitySamples);
      RP::get("SilvaShock.maxwCutoff", maxwCutoff);
      RP::get("SilvaShock.Scale_x", SCA_X);
      RP::get("SilvaShock.Scale_y", SCA_Y);
      RP::get("SilvaShock.Width", Shockwidth);

      RP::get("SilvaShock.NoiseStrength", noise_strength);
      RP::get("SilvaShock.NoiseScale", noise_scale);
      RP::get("SilvaShock.NoiseOctaves", noise_octaves);

      std::cerr << "Noise Strength = " << noise_strength << std::endl;
   }

   vector<std::array<Real, 3>> SilvaShock::getV0(creal x, creal y, creal z) {
      Real VX0 = interpolate(VX0u,VX0d, x);
      Real VY0 = interpolate(VY0u,VY0d, x);
      Real VZ0 = interpolate(VZ0u,VZ0d, x);

      std::array<Real, 3> V0{{VX0, VY0, VZ0}};
      vector<std::array<Real, 3>> retval;
      retval.push_back(V0);

      return retval;
   }

   Real SilvaShock::getDistribValue(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz) {
      creal k = 1.3806505e-23; // Boltzmann
      creal mass = 1.67262171e-27; // m_p in kg

      Real VX0 = interpolate(VX0u,VX0d, x);
      Real VY0 = interpolate(VY0u,VY0d, x);
      Real VZ0 = interpolate(VZ0u,VZ0d, x);
      Real TEMPERATURE = interpolate(TEMPERATUREu,TEMPERATUREd, x);

      return exp(- mass * ((vx-VX0)*(vx-VX0) + (vy-VY0)*(vy-VY0)+ (vz-VZ0)*(vz-VZ0)) / (2.0 * k * TEMPERATURE));
      //*exp(-pow(x-Parameters::xmax/2.0, 2.0)/pow(this->SCA_X, 2.0))*exp(-pow(y-Parameters::ymax/4.0, 2.0)/pow(this->SCA_Y, 2.0));
   }


   Real SilvaShock::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx, creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz) {
      if(vx < Parameters::vxmin + 0.5 * dvx ||
         vy < Parameters::vymin + 0.5 * dvy ||
         vz < Parameters::vzmin + 0.5 * dvz ||
         vx > Parameters::vxmax - 1.5 * dvx ||
         vy > Parameters::vymax - 1.5 * dvy ||
         vz > Parameters::vzmax - 1.5 * dvz
      ) return 0.0;
      
      creal mass = physicalconstants::MASS_PROTON;
      creal q = physicalconstants::CHARGE;
      creal k = 1.3806505e-23; // Boltzmann
      creal mu0 = 1.25663706144e-6; // mu_0
      Real DENSITY = interpolate(DENSITYu,DENSITYd,x);

      creal d_x = dx / (nSpaceSamples-1);
      creal d_y = dy / (nSpaceSamples-1);
      creal d_z = dz / (nSpaceSamples-1);
      creal d_vx = dvx / (nVelocitySamples-1);
      creal d_vy = dvy / (nVelocitySamples-1);
      creal d_vz = dvz / (nVelocitySamples-1);
      Real avg = 0.0;
      
      for (uint i=0; i<nSpaceSamples; ++i)
      for (uint j=0; j<nSpaceSamples; ++j)
      for (uint k=0; k<nSpaceSamples; ++k)      
      for (uint vi=0; vi<nVelocitySamples; ++vi)
      for (uint vj=0; vj<nVelocitySamples; ++vj)
      for (uint vk=0; vk<nVelocitySamples; ++vk)
            {
         avg += getDistribValue(x+i*d_x, y+j*d_y, z+k*d_z, vx+vi*d_vx, vy+vj*d_vy, vz+vk*d_vz);
         }
      
      Real TEMPERATURE = interpolate(TEMPERATUREu,TEMPERATUREd, x);
      creal result = avg *DENSITY * pow(mass / (2.0 * M_PI * k * TEMPERATURE), 1.5) /
                     (nSpaceSamples*nSpaceSamples*nSpaceSamples) / 
   //   	            (Parameters::vzmax - Parameters::vzmin) / 
                     (nVelocitySamples*nVelocitySamples*nVelocitySamples);
               
               
      if(result < maxwCutoff) {
         return 0.0;
      } else {
         return result;
      }
   }

   void SilvaShock::calcCellParameters(Real* cellParams,creal& t) {
      creal x = cellParams[CellParams::XCRD];
      creal dx = cellParams[CellParams::DX];
      creal y = cellParams[CellParams::YCRD];
      creal dy = cellParams[CellParams::DY];
      creal z = cellParams[CellParams::ZCRD];
      creal dz = cellParams[CellParams::DZ];
      
      cellParams[CellParams::EX   ] = interpolate(EX0u,EX0d,x);
      cellParams[CellParams::EY   ] = interpolate(EY0u,EY0d,x);
      cellParams[CellParams::EZ   ] = interpolate(EZ0u,EZ0d,x);

      Vec3d pn = noise_strength * divergence_free_noise(Vec3d(x,y,z), Vec3d(noise_scale), noise_octaves);
      // Window noise
      pn *= .5*(1. - cos(2.*M_PI*(y-Parameters::ymin)/(Parameters::ymax - Parameters::ymin)));
      pn *= .5*(1. - cos(2.*M_PI*(x-Parameters::xmin)/(Parameters::xmax - Parameters::xmin)));

      cellParams[CellParams::PERBX   ] = interpolate(BX0u,BX0d,x) + pn[0];
      cellParams[CellParams::PERBY   ] = interpolate(BY0u,BY0d,x) + pn[1];
      cellParams[CellParams::PERBZ   ] = interpolate(BZ0u,BZ0d,x) + pn[2];
   }

   Real SilvaShock::interpolate(Real upstream, Real downstream, Real x) {
	x /= 2. * Shockwidth;
	Real a = .5 * (1. + tanh(x * 2. * M_PI));

	return upstream * a + downstream * (1. - a);
   }
}//namespace projects
