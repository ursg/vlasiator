#include "cuda_build_acc_columns.cuh"

#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <thrust/execution_policy.h>

struct isZero{
   __host__ __device__ bool operator()(const uint x) {
      return x == 0;
   }
};

__host__ void buildAccColumns(OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>& vmesh, uint dimension, cudaStream_t stream) {

   size_t numBlocks = vmesh.size();

   thrust::device_vector<vmesh::GlobalID> sortedBlockMappedGID(numBlocks);
   thrust::device_vector<vmesh::LocalID> sortedBlockLID(numBlocks);

   // Sort velocity blocks by the relevant direction
   //thrust::sort_by_key(thrust::cuda::par.on(stream),
   //      sortedBlockMappedGID.begin(), sortedBlockMappedGID.end(),
   //      sortedBlockLID);

   // Tag which blocks are start-of-column
   thrust::device_vector<vmesh::LocalID> columnStartLID(numBlocks);

   // Stream-compact the column start indices
   auto newEnd = thrust::remove_if(thrust::cuda::par.on(stream), 
         columnStartLID.begin(), columnStartLID.end(),
         isZero());

   int numColumns = newEnd - columnStartLID.begin(); 

   // Run along the columns and identify their stop
}
