/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::cout << "Initialize" << std::endl;
  
  num_particles = 100;  // TODO: Set the number of particles

  std::default_random_engine gen;
  //sample center around x,y,theta from gps and their standard deviation
  std::normal_distribution<double> dist_x(x, std[0]); 
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  // for each particle
  for(int i = 0; i < num_particles; i++){
    
    // create a particle around that mean with given std
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;

    particles.push_back(p);
    weights.push_back(p.weight);

  }

  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  std::cout << "Prediction" << std::endl;
  std::default_random_engine gen;

  // for each particles
  for(int i = 0; i < num_particles; i++){

    double prev_x = particles[i].x;
    double prev_y = particles[i].y;
    double prev_theta = particles[i].theta;
    double curr_x, curr_y, curr_theta;

    // set new position after movement
    // if yaw_rate is really small, just ignore
    if(fabs(yaw_rate) < 0.0001){
      curr_x = prev_x + velocity * delta_t * cos(prev_theta);
      curr_y = prev_y + velocity * delta_t * sin(prev_theta);
      curr_theta = prev_theta;
    }
    else{
      curr_x = prev_x + velocity * (sin(prev_theta + yaw_rate * delta_t) - sin(prev_theta)) / yaw_rate;
      curr_y = prev_y + velocity * (cos(prev_theta) - cos(prev_theta + yaw_rate * delta_t)) / yaw_rate;
      curr_theta = prev_theta + yaw_rate * delta_t;
    }
    
    // add some noise
    std::normal_distribution<double> dist_x(curr_x, std_pos[0]); 
    std::normal_distribution<double> dist_y(curr_y, std_pos[1]);
    std::normal_distribution<double> dist_theta(curr_theta, std_pos[2]);

    // update the state of the particle
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  std::cout << "Update Weights" << std::endl;

  // for each particle
  for(int i = 0; i < num_particles; i++){
    double px = particles[i].x;
    double py = particles[i].y;
    double ptheta = particles[i].theta;

    // find only those landmark within sensing range, store them in predictions
    vector<LandmarkObs> prediction;
    for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
      int lm_id = map_landmarks.landmark_list[j].id_i;
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      //if landmarks is beyond the particle sensor distance, ignore it
      if(dist(lm_x, lm_y, px, py) <= sensor_range){
        prediction.push_back(LandmarkObs{lm_id, lm_x, lm_y});
      }
    }
    
    double final_weight = 1.0;

    // update weight

    // for each observation from real car
    for(int j = 0; j < observations.size(); j++){

      // convert observation into map coor
      double obs_x, obs_y, obs_theta;
      obs_x = px + cos(ptheta)*observations[j].x - sin(ptheta)*observations[j].y;
      obs_y = py + cos(ptheta)*observations[j].y + sin(ptheta)*observations[j].x;

      double nearest_x, nearest_y, pred_dist_min;
      pred_dist_min = 10000000;

      // find the closest predition to that observation map coor
      for(int k = 0; k < prediction.size(); k++){
        double calc_dist = dist(obs_x, obs_y, prediction[k].x, prediction[k].y);
        if(calc_dist < pred_dist_min){
          pred_dist_min = calc_dist;
          nearest_x = prediction[k].x;
          nearest_y = prediction[k].y;
        }
      }

      // using multivariate gaussian dist to find prob from the error, then update weight
      double prob = 1/(2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-1 * (pow(obs_x - nearest_x,2)/(2*pow(std_landmark[0],2)) + pow(obs_y - nearest_y,2)/(2*pow(std_landmark[1],2))));
      final_weight *= prob;
    }
    // update particle weight
    particles[i].weight = final_weight;
    weights[i] = final_weight;
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::cout << "Resampling" << std::endl;
  vector<Particle> new_particles;

  std::default_random_engine generator;
  std::discrete_distribution<int> distribution(weights.begin(), weights.end());
  for(int i = 0; i < num_particles; i++){
    int num = distribution(generator);
    Particle p = particles[num];
    p.id = i;
    new_particles.push_back(p);
  }

  particles = new_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}