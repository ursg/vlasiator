/*
This file is part of Vlasiator.

Copyright 2011, 2012 Finnish Meteorological Institute

*/

#pragma once

#include "../../definitions.h"
#include "../projectTriAxisSearch.h"

namespace projects {
   class SilvaShock: public TriAxisSearch {
      public:
         SilvaShock();
         virtual ~SilvaShock();
      
         virtual bool initialize(void);
         static void addParameters(void);
         virtual void getParameters(void);
      
      protected:
         Real getDistribValue(
            creal& x,creal& y, creal& z,
            creal& vx, creal& vy, creal& vz
         ); 
         virtual vector<std::array<Real, 3>> getV0(creal x, creal y, creal z);
         virtual void calcCellParameters(Real* cellParams,creal& t);
         virtual Real calcPhaseSpaceDensity(
            creal& x, creal& y, creal& z,
            creal& dx, creal& dy, creal& dz,
            creal& vx, creal& vy, creal& vz,
            creal& dvx, creal& dvy, creal& dvz
         );

         // Interpolate between up- and downstream quantities
         // based on position
         Real interpolate(Real u, Real d, Real x);

         // Upstream bulk values
         Real BX0u;
         Real BY0u;
         Real BZ0u;
         Real EX0u;
         Real EY0u;
         Real EZ0u;
         Real VX0u;
         Real VY0u;
         Real VZ0u;
         Real DENSITYu;
         Real TEMPERATUREu;

         // Downstream bulk values
         Real BX0d;
         Real BY0d;
         Real BZ0d;
         Real EX0d;
         Real EY0d;
         Real EZ0d;
         Real VX0d;
         Real VY0d;
         Real VZ0d;
         Real DENSITYd;
         Real TEMPERATUREd;

         Real maxwCutoff;
         uint nSpaceSamples;
         uint nVelocitySamples;
         Real SCA_X;
         Real SCA_Y;
         Real Shockwidth;

         Real noise_strength;
         Real noise_scale;
         uint noise_octaves;
   } ; //class Shock
} // namespace projects
