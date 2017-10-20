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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  std::cout << "ParticleFilter::init()" << std::endl;
	
	num_particles = 10;
	weights.resize(num_particles);
	
	std::default_random_engine generator;	 
	std::normal_distribution<double> nd_x(x, std[0]);
	std::normal_distribution<double> nd_y(y, std[1]);
	std::normal_distribution<double> nd_theta(theta, std[2]);
	
	for (int i=0; i<num_particles; i++) {
		// Create new particle and push to class vector
		Particle new_particle;
		
		new_particle.id = i;
		new_particle.x = nd_x(generator);
		new_particle.y = nd_y(generator);
		new_particle.theta = nd_theta(generator);
		new_particle.weight = 1.0;

    particles.push_back(new_particle);
	}

	std::cout << "  Number of particles: " << std::to_string(num_particles) << std::endl;

	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  std::cout << "ParticleFilter::prediction()" << std::endl;
	
	std::default_random_engine generator;
	std::normal_distribution<double> nd_x(0.0, std_pos[0]);
	std::normal_distribution<double> nd_y(0.0, std_pos[1]);
	std::normal_distribution<double> nd_theta(0.0, std_pos[2]);

	const double yaw_rate_min = 1e-4;
	if (fabs(yaw_rate) < yaw_rate_min) {
		yaw_rate = ((yaw_rate<0) ? (-1*yaw_rate_min) : (yaw_rate_min));
	}

	for (auto&& it_particle : particles) {
				// Predict and convert coordinate system from vehicle to map
				double map_coord_x = sin(it_particle.theta + yaw_rate * delta_t) - sin(it_particle.theta);
				double map_coord_y = cos(it_particle.theta) - cos(it_particle.theta + yaw_rate * delta_t);
				
				// Calculate the predicted coordinates
        it_particle.x += (velocity / yaw_rate) * map_coord_x + nd_x(generator);
        it_particle.y += (velocity / yaw_rate) * map_coord_y + nd_y(generator);
        it_particle.theta += yaw_rate * delta_t + nd_theta(generator);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  std::cout << "ParticleFilter::dataAssociation()" << std::endl;
	
	std::vector<LandmarkObs> associations;
	
	// Find closest association
	for (auto& it_observations : observations) {
		double min_dist = std::numeric_limits<double>::max();
		LandmarkObs min_association;
		
		// Iterate through predicted and calculate distance
		for (auto& it_predicted : predicted) {
			double dist = pow(it_predicted.x - it_observations.x, 2) + pow(it_predicted.y - it_observations.y, 2);
			if (dist < min_dist) {
				min_dist = dist;
				min_association = it_predicted;
			}
		}
		
		// Push closest association to vector
		associations.push_back(min_association);
	}
	
	// Print association
	for (auto& it_associations : associations) {
		std::cout << it_associations.x << " " << it_associations.y << std::endl;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
  std::cout << "ParticleFilter::updateWeights()" << std::endl;

	const double var_x = pow(std_landmark[0], 2);
	const double var_y = pow(std_landmark[1], 2);
	const double covar_xy = std_landmark[0] * std_landmark[1];
	double weights_sum = 0;	
	
	// Iterate through particles and update weight
	//for (auto& it_particles : particles) {
	// Use old iterator style to be able to use std::distance() later on
	for (std::vector<Particle>::iterator it_particles = particles.begin(); 
	     it_particles != particles.end(); 
			 it_particles++) {
		double weight = 1.0;
		
		// Iterate through observations
		for (auto& it_observations : observations) {
			// Calcluate the predicted coordinates
			double predicted_x = it_observations.x * cos(it_particles->theta) - 
												   it_observations.y * sin(it_particles->theta) + 
													 it_particles->x;
			double predicted_y = it_observations.x * sin(it_particles->theta) + 
													 it_observations.y * cos(it_particles->theta) + 
													 it_particles->y;			
	
			// Find nearest landmark
			Map::single_landmark_s nearest_landmark;
			double min_distance = std::numeric_limits<double>::max();
			double distance;
			
			for (auto& it_landmark_list : map_landmarks.landmark_list) {
				distance = fabs(predicted_x - it_landmark_list.x_f) + fabs(predicted_y - it_landmark_list.y_f);
			
				if (distance < min_distance) {
					min_distance = distance;
					nearest_landmark = it_landmark_list;
				}
			}
			
			// Update weight with double variate gaussian
			double x_diff = predicted_x - nearest_landmark.x_f;
			double y_diff = predicted_y - nearest_landmark.y_f;
			double num = exp(-0.5*((x_diff * x_diff)/var_x + (y_diff * y_diff)/var_y));
			double denom = 2*M_PI*covar_xy;

			weight = weight * (num / denom);
		} 
		
		it_particles->weight = weight;
		int idx_it_particles = std::distance(particles.begin(), it_particles);
		weights[idx_it_particles] = weight;
	}
}
		
void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  std::cout << "ParticleFilter::resample()" << std::endl;
	
	std::default_random_engine generator;
	std::discrete_distribution<> weights_pmf(weights.begin(), weights.end());	
	std::vector<Particle> newParticles;

	for (int i = 0; i < num_particles; ++i) {
			newParticles.push_back(particles[weights_pmf(generator)]);
	}	
	
	particles = newParticles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates
	
	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

std::string ParticleFilter::getAssociations(Particle best)
{
	std::vector<int> v = best.associations;
	std::stringstream ss;
    copy( v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
    std::string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
std::string ParticleFilter::getSenseX(Particle best)
{
	std::vector<double> v = best.sense_x;
	std::stringstream ss;
    copy( v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
    std::string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
std::string ParticleFilter::getSenseY(Particle best)
{
	std::vector<double> v = best.sense_y;
	std::stringstream ss;
    copy( v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
    std::string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
