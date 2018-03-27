//
// Created by ymlai on 5/3/2018.
//

#ifndef PARTICLE_FILTER_KDTREE_H
#define PARTICLE_FILTER_KDTREE_H

#include "helper_functions.h"
#include <stdlib.h>

#define DIM 2
class KDTree {
public:
  KDTree(std::vector<LandmarkObs> landmarks);

  ~KDTree();

  LandmarkObs find_nearest(LandmarkObs point);
};


#endif //PARTICLE_FILTER_KDTREE_H
