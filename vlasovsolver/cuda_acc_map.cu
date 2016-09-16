#include "../velocity_mesh_cuda.h"
#include "cuda_acc_map.h"

#warning "Integrate this Realv with the vec.h machinery"
#define Realv float

__global__ void map_1d(vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_sourceVmesh,
                       vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_targetVmesh,
                       Realv intersection,
                       Realv intersection_di,
                       Realv intersection_dj,
                       Realv intersection_dk,
                       uint dimension);

__global__ void printOutput(vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_Vmesh,
                            int tag) {
   const vmesh::LocalID columnStart = d_Vmesh->columnStartLID[0];
   const vmesh::LocalID sortedBlockLID = d_Vmesh->sortedBlockLID[columnStart];
   printf("=== %i ===\n", tag);
   for (uint i=0; i<64; i+=8) {
      printf("%d_%e\t", i, d_Vmesh->data[sortedBlockLID * WID3 + i]);
   }
   printf("\n============\n");
}

bool map3DCuda(Realf **blockDatas,
               vmesh::GlobalID **blockIDs,
               vmesh::LocalID *nBlocks,
               Real *intersections,
               const uint nCells,
               const Realf blockSize[3],
               const vmesh::LocalID gridLength[3],
               const Realf gridMinLimits[3],
               const float sparsityThreshold
               ){
   fprintf(stderr," - - - Starting map3DCuda - - -\n");
   bool success = true;
   cudaStream_t streams[nCells];
   vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_sourceVmesh[nCells];
   vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* h_sourceVmesh[nCells];
   
   vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_targetVmesh[nCells];
   vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* h_targetVmesh[nCells];

   fprintf(stderr," - - - Allocating memory - - -\n");

   /*allocate memory for all cells, these operations are blocking*/
   //TODO: add checks/throttling to make sure there is enough memory on device...

   for (int i = 0; i < nCells; i++) {
      cudaStreamCreate(&streams[i]);
      vmesh::createVelocityMeshCuda(&(d_sourceVmesh[i]), &(h_sourceVmesh[i]), 
            nBlocks[i], gridLength, blockSize, gridMinLimits);
   }
   
   fprintf(stderr," - - - Meshes created. Now, the loop: - - -\n");
   for (int i = 0; i < nCells; i++) {
      fprintf(stderr," `-> transferDataHostToDevice\n");
      vmesh::transferDataHostToDevice(d_sourceVmesh[i], h_sourceVmesh[i], blockDatas[i], blockIDs[i], streams[i]);
     //Order Z X Y
      //DO Z (REMEMBER TO REAORDER INTERSECTIONS FOR OTHER DIMENSIONS)
      fprintf(stderr," `-> velocityBlocksInColumns\n");
      vmesh::sortVelocityBlocksInColumns(d_sourceVmesh[i], h_sourceVmesh[i], 2, streams[i]);
      
      fprintf(stderr," `-> createTargetMesh\n");
      vmesh::createTargetMesh(&(d_targetVmesh[i]), &(h_targetVmesh[i]),
                              d_sourceVmesh[i], h_sourceVmesh[i],
                              intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z],
                              intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z_DI],
                              intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z_DJ],
                              intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z_DK],
                              2, streams[i]);
      printOutput<<<1,1>>>(d_sourceVmesh[i], 0);
      map_1d<<<h_sourceVmesh[i]->nColumns, dim3 (4, 4, 1)>>>(d_sourceVmesh[i],
                                                             d_sourceVmesh[i],
                                                             intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z],
                                                             intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z_DI],
                                                             intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z_DJ],
                                                             intersections[i * AccelerationIntersections::N_INTERSECTIONS + AccelerationIntersections::Z_DK],
                                                             2);
      printOutput<<<1,1>>>(d_sourceVmesh[i], 1);
      vmesh::transferDataDeviceToHost(d_sourceVmesh[i], h_sourceVmesh[i], blockDatas[i], blockIDs[i], streams[i]);
      cudaDeviceSynchronize();
      
//      vmesh::sortVelocityBlocksInColumns(d_sourceVmesh[i], h_sourceVmesh[i], 0, streams[i]);
//      vmesh::sortVelocityBlocksInColumns(d_sourceVmesh[i], h_sourceVmesh[i], 1, streams[i]);
      cudaDeviceSynchronize();
      fprintf(stderr," `-> adjustVelocityBlocks\n");
   }

   vmesh::adjustVelocityBlocks(d_sourceVmesh, h_sourceVmesh, sparsityThreshold, streams, nCells);
   
   for (int i = 0; i < nCells; i++) {
      vmesh::destroyVelocityMeshCuda(d_sourceVmesh[i], h_sourceVmesh[i]);
      cudaStreamDestroy(streams[i]);
   }
   
   //cudaDeviceSynchronize();
   return success;
}

__global__ void map_1d(vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_sourceVmesh,
                       vmesh::VelocityMeshCuda<vmesh::GlobalID, vmesh::LocalID>* d_targetVmesh,
                       Realv intersection,
                       Realv intersection_di,
                       Realv intersection_dj,
                       Realv intersection_dk,
                       uint dimension) {
   if (   blockIdx.x > d_sourceVmesh->nColumns
       || threadIdx.x >= WID
       || threadIdx.y >= WID 
   ) {
      return;
   }
   const vmesh::LocalID columnStart = d_sourceVmesh->columnStartLID[blockIdx.x];
   
   for (uint k=0; k<d_sourceVmesh->columnSize(blockIdx.x); k++) {
      const vmesh::LocalID sortedBlockLID = d_sourceVmesh->sortedBlockLID[columnStart + k];
      const uint threadIdxBXY = sortedBlockLID * WID3 + threadIdx.x + threadIdx.y * WID;
      
      for (uint kz=0; kz<4; kz++) {
         d_targetVmesh->data[threadIdxBXY + kz * WID2] = 10.0 * d_sourceVmesh->data[threadIdxBXY + kz * WID2];
      }
   }
   
}
