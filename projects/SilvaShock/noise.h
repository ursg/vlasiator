#pragma once
#include "simplexnoise1234.h"
#include "vectorclass.h"
#include "vector3d.h"

Vec3d divergence_free_noise(Vec3d p, Vec3d scales, int num_octaves, double axscale=1.);
