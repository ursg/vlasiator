/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2014 Finnish Meteorological Institute
 */

#ifndef VELOCITY_MESH_CUDA_H
#define VELOCITY_MESH_CUDA_H

#include "common.h"
#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#warning "Integrate this Realf with the vec.h machinery"
#define Realf float
#else
#define HOST_DEVICE
#endif

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/functional.h>

//

namespace vmesh {
   template<typename GID>
   class VelocityMeshCuda {
   public:
      __host__ void init(Realf *h_data, GID *h_blockIDs, uint nBlocks);
      
   private:
      Realf dv;      
      uint nBlocks;
      Realf *data;
      GID *blockIDs;      
   };
   

   /*init on host side*/
   template<typename GID> __host__ void VelocityMeshCuda<GID>::init(Realf *h_data, GID *h_blockIDs, uint h_nBlocks) {
      cudaEvent_t evA, evB;
      cudaEventCreate(&evA);
      cudaEventCreate(&evB);
      
      nBlocks=h_nBlocks;
      cudaMalloc(&data, nBlocks * WID3 * sizeof(Realf));
      cudaMalloc(&blockIDs, nBlocks * sizeof(GID));
      
      
      cudaEventRecord(evA);      
      cudaMemcpy(data, h_data, nBlocks *  WID3 * sizeof(Realf), cudaMemcpyHostToDevice);
      cudaMemcpy(blockIDs, h_blockIDs, nBlocks * sizeof(GID), cudaMemcpyHostToDevice);
      cudaEventRecord(evB);
      cudaEventSynchronize(evB);
      uint bytes = nBlocks *  WID3 * sizeof(Realf) +nBlocks * sizeof(GID);
      float milliseconds=0;
      cudaEventElapsedTime(&milliseconds, evA, evB);  
      printf("Transferred %d blocks,  %d bytes in %g s to GPU: %g GB/s. \n", nBlocks, bytes, milliseconds * 1e-3, (bytes * 1e-9) / (milliseconds * 1e-3) );
      
      
   }

   template<typename GID> __host__ VelocityMeshCuda<GID>* initVelocityMeshCuda(Realf *h_data, GID *h_blockIDs, uint nBlocks) {
      VelocityMeshCuda<GID> h_vmesh;
      VelocityMeshCuda<GID> *d_vmesh;
      //allocate space on device for device resident class
      cudaMalloc((void **)&d_vmesh, sizeof(VelocityMeshCuda<GID>));      
      //init members on host
      h_vmesh.init(h_data, h_blockIDs, nBlocks);
      //copy all  members to device
      cudaMemcpy(d_vmesh, &h_vmesh, sizeof(VelocityMeshCuda<GID>), cudaMemcpyHostToDevice);
      
      
      return d_vmesh;
   }


}; // namespace vmesh

#endif
