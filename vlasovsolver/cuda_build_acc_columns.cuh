#pragma once
#include "vec.h"
#include "../definitions.h"

#include "device_launch_parameters.h"
#include "cuda.h"
#include "cuda_runtime.h"

#include "../open_bucket_hashtable.h"

__host__ void buildAccColumns(OpenBucketHashtable<vmesh::GlobalID,vmesh::LocalID>& vmesh, uint dimension, cudaStream_t stream);
