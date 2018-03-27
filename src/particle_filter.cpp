/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  num_particles = 100;


  double std_x, std_y, std_theta; // Standard deviations for x, y, and theta

  // TODO: Set standard deviations for x, y, and theta
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];

  // This line creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x, std_x);
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);


  for (int i = 0; i < num_particles; ++i) {
    double sample_x, sample_y, sample_theta;

    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);

    Particle p;
    p.id = i;
    p.x = sample_x;
    p.y = sample_y;
    p.theta = sample_theta;
    p.weight = 1;

    particles.push_back(p);
    // Print your samples to the terminal.
    //cout << "Sample " << i + 1 << " " << sample_x << " " << sample_y << " " << sample_theta << endl;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  normal_distribution<double> N_x(0, std_pos[0]);
  normal_distribution<double> N_y(0, std_pos[1]);
  normal_distribution<double> N_theta(0, std_pos[2]);

  for(int i=0; i<num_particles; ++i){
    double curr_x = particles[i].x, curr_y = particles[i].y, curr_theta = particles[i].theta;
    double next_x, next_y, next_theta;

    if(fabs(yaw_rate) < 1e-6){
      next_x = curr_x + velocity * delta_t * cos(curr_theta);
      next_y = curr_y + velocity * delta_t * sin(curr_theta);
      next_theta = curr_theta;
    }
    else{
      next_x = curr_x + velocity/yaw_rate * (sin(curr_theta + yaw_rate * delta_t) - sin(curr_theta));
      next_y = curr_y + velocity/yaw_rate * (cos(curr_theta) - cos(curr_theta + yaw_rate * delta_t));
      next_theta = curr_theta + yaw_rate * delta_t;
    }

    next_x += N_x(gen);
    next_y += N_y(gen);
    next_theta += N_theta(gen);

    particles[i].x = next_x;
    particles[i].y = next_y;
    particles[i].theta = next_theta;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  //can improve from kd tree
  for (unsigned int i = 0; i < observations.size(); i++) {

    // grab current observation
    LandmarkObs o = observations[i];

    // init minimum distance to maximum possible
    double min_dist = numeric_limits<double>::max();

    // init id of landmark from map placeholder to be associated with the observation
    int map_id = -1;

    for (unsigned int j = 0; j < predicted.size(); j++) {
      // grab current prediction
      LandmarkObs p = predicted[j];

      // get distance between current/predicted landmarks
      double cur_dist = dist(o.x, o.y, p.x, p.y);

      // find the predicted landmark nearest the current observed landmark
      if (cur_dist < min_dist) {
        min_dist = cur_dist;
        map_id = p.id;
      }
    }

    // set the observation's id to the nearest predicted landmark's id
    observations[i].id = map_id;
  }
}

vector<LandmarkObs> transformObservationFromCarToWorld(double particle_x, double particle_y, double particle_theta, const vector<LandmarkObs> &observations){
  vector<LandmarkObs> transformed_observations;

  for (int j = 0; j < observations.size(); j++) {
    LandmarkObs transformed_obs;
    transformed_obs.id = j;
    transformed_obs.x = particle_x + (cos(particle_theta) * observations[j].x) - (sin(particle_theta) * observations[j].y);
    transformed_obs.y = particle_y + (sin(particle_theta) * observations[j].x) + (cos(particle_theta) * observations[j].y);
    transformed_observations.push_back(transformed_obs);
  }
  return transformed_observations;
}

vector<LandmarkObs> filterLandmarkInSensorRange(double particle_x, double particle_y, double sensor_range, const Map &map_landmarks){
  vector<LandmarkObs> filtered_landmarks;

  //can improve from kd tree
  for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
    Map::single_landmark_s current_landmark = map_landmarks.landmark_list[j];

    if (dist(particle_x, particle_y, current_landmark.x_f, current_landmark.y_f) <= sensor_range){
      filtered_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
    }
  }

  return filtered_landmarks;
}

double updateParticleWeight(const vector<LandmarkObs> & transformed_observations, const vector<LandmarkObs> & observations, double std_landmark[2]){
  double weight = 1.0;

  double sigma_x = std_landmark[0];
  double sigma_y = std_landmark[1];
  double sigma_x_2 = pow(sigma_x, 2);
  double sigma_y_2 = pow(sigma_y, 2);
  double normalizer = (1.0/(2.0 * M_PI * sigma_x * sigma_y));
  int k, l;

  for (k = 0; k < transformed_observations.size(); k++) {
    double trans_obs_x = transformed_observations[k].x;
    double trans_obs_y = transformed_observations[k].y;
    double trans_obs_id = transformed_observations[k].id;
    double multi_prob;

    for (l = 0; l < observations.size(); l++) {
      double pred_landmark_x = observations[l].x;
      double pred_landmark_y = observations[l].y;
      double pred_landmark_id = observations[l].id;

      if (trans_obs_id == pred_landmark_id) {
        multi_prob = normalizer * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y_2))));
        weight *= multi_prob;
      }
    }
  }

  return weight;
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  double total_weight = 0.0;
  for(int i=0; i<num_particles; ++i){
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;

    vector<LandmarkObs> transformed_observation = transformObservationFromCarToWorld(particle_x, particle_y, particle_theta, observations);

    vector<LandmarkObs> filtered_observation = filterLandmarkInSensorRange(particle_x, particle_y, sensor_range, map_landmarks);

    dataAssociation(filtered_observation, transformed_observation);

    double new_weight = updateParticleWeight(transformed_observation, filtered_observation, std_landmark);

    particles[i].weight = new_weight;
    total_weight += particles[i].weight;
  }

  for(int i=0; i<particles.size(); ++i){
    particles[i].weight /= total_weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  vector<Particle> new_particles;

  // get all of the current weights
  weights.clear();
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  std::discrete_distribution<> d(weights.begin(), weights.end());
  for(int i=0; i<num_particles; ++i ) {
    new_particles.push_back(particles[d(gen)]);
  }

  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
