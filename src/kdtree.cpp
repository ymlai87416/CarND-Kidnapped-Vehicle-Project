//
// Created by ymlai on 5/3/2018.
//

#include "kdtree.h"
#include <stdio.h>
#include <string.h>

#define MAX_DIM 2

struct kd_node_t{
  double x[MAX_DIM];
  struct kd_node_t *left, *right;
};

struct kd_node_t* root;
struct kd_node_t* landmarks_node;

inline double
dist(struct kd_node_t *a, struct kd_node_t *b, int dim)
{
  double t, d = 0;
  while (dim--) {
    t = a->x[dim] - b->x[dim];
    d += t * t;
  }
  return d;
}
inline void swap(struct kd_node_t *x, struct kd_node_t *y) {
  double tmp[MAX_DIM];
  memcpy(tmp,  x->x, sizeof(tmp));
  memcpy(x->x, y->x, sizeof(tmp));
  memcpy(y->x, tmp,  sizeof(tmp));
}


/* see quickselect method */
struct kd_node_t*
find_median(struct kd_node_t *start, struct kd_node_t *end, int idx)
{
  if (end <= start) return NULL;
  if (end == start + 1)
    return start;

  struct kd_node_t *p, *store, *md = start + (end - start) / 2;
  double pivot;
  while (1) {
    pivot = md->x[idx];

    swap(md, end - 1);
    for (store = p = start; p < end; p++) {
      if (p->x[idx] < pivot) {
        if (p != store)
          swap(p, store);
        store++;
      }
    }
    swap(store, end - 1);

    /* median has duplicate values */
    if (store->x[idx] == md->x[idx])
      return md;

    if (store > md) end = store;
    else        start = store;
  }
}

struct kd_node_t*
make_tree(struct kd_node_t *t, int len, int i, int dim)
{
  struct kd_node_t *n;

  if (!len) return 0;

  if ((n = find_median(t, t + len, i))) {
    i = (i + 1) % dim;
    n->left  = make_tree(t, n - t, i, dim);
    n->right = make_tree(n + 1, t + len - (n + 1), i, dim);
  }
  return n;
}

/* global variable, so sue me */
int visited;

void nearest(struct kd_node_t *root, struct kd_node_t *nd, int i, int dim,
             struct kd_node_t **best, double *best_dist)
{
  double d, dx, dx2;

  if (!root) return;
  d = dist(root, nd, dim);
  dx = root->x[i] - nd->x[i];
  dx2 = dx * dx;

  visited ++;

  if (!*best || d < *best_dist) {
    *best_dist = d;
    *best = root;
  }

  /* if chance of exact match is high */
  if (!*best_dist) return;

  if (++i >= dim) i = 0;

  nearest(dx > 0 ? root->left : root->right, nd, i, dim, best, best_dist);
  if (dx2 >= *best_dist) return;
  nearest(dx > 0 ? root->right : root->left, nd, i, dim, best, best_dist);
}

KDTree::KDTree(std::vector<LandmarkObs> landmarks) {
  auto n = landmarks.size();
  landmarks_node = (struct kd_node_t*) calloc(n, sizeof(struct kd_node_t));

  for(int i=0; i<n; ++i) {
    landmarks_node[i].x[0] = landmarks[i].x;
    landmarks_node[i].x[1] = landmarks[i].y;
  }

  root = make_tree(landmarks_node, n, 0, 2);
}

LandmarkObs KDTree::find_nearest(LandmarkObs point) {
  struct kd_node_t point_t;
  struct kd_node_t *found;
  double best_dist = 0;

  point_t.x[0] = point.x;
  point_t.x[1] = point.y;

  nearest(root, &point_t, 0, 2, &found, &best_dist);
}

KDTree::~KDTree() {

  free(landmarks_node);
}