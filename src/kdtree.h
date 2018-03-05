//
// Created by ymlai on 5/3/2018.
//

#ifndef PARTICLE_FILTER_KDTREE_H
#define PARTICLE_FILTER_KDTREE_H

#include "helper_functions.h"

#define DIM 2
class KDTree {
public:
  KDTree(std::vector<LandmarkObs> landmarks);

  ~KDTree();

  LandmarkObs find_nearest(LandmarkObs point);

private:
  struct kd_node_t{
    double x[DIM];
    struct kd_node_t *left, *right;
  };

  struct kd_node_t* root;
  struct kd_node_t* landmarks_node;
};


#endif //PARTICLE_FILTER_KDTREE_H
