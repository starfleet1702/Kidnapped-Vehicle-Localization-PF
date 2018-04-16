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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// DONE1: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	
	std::cout<<"Init Start"<<std::endl;
	
	//Number of Particles to draw
	num_particles = 50;
	std::default_random_engine gen;
	
	//Normal Distributions to represent uncertatiniy in GPS Measurements
	normal_distribution<double> dist_noise_x (x,std[0]);
	normal_distribution<double> dist_noise_y (y,std[1]);
	normal_distribution<double> dist_noise_theta (theta,std[2]);
	
	//Initializing Particles
	for(unsigned int i=0; i<num_particles; i++){
		Particle p;
		p.id = i;
		p.x = x + dist_noise_x(gen);
		p.y = y + dist_noise_y(gen);
		p.theta = theta + dist_noise_theta(gen);
		p.weight = 1.0;
		particles.push_back(p);
	}
	
	//done with Initialization
	is_initialized = true;
	
	std::cout<<"Init End"<<endl;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// DONE: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise std::normal_distribution and std::default_random_engine is used
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	std::cout<<"prediction Start"<<std::endl;
	
	double min_threshold = 0.0001;
	double v_by_yaw_rate = (velocity/yaw_rate);
	double yaw_rate_mul_dt =  yaw_rate*delta_t;
	double v_mul_dt = velocity*delta_t;
	
	std::default_random_engine gen;
	
	normal_distribution<double> dist_noise_x (0,std_pos[0]);
	normal_distribution<double> dist_noise_y (0,std_pos[1]);
	normal_distribution<double> dist_noise_theta (0,std_pos[2]);
	
	for(unsigned int i=0; i<num_particles; i++){
		Particle p = particles[i];
		//Considering bicycle CTRV motion model
		if(fabs(yaw_rate) > min_threshold){
			p.x += v_by_yaw_rate*(sin(p.theta+yaw_rate_mul_dt)-sin(p.theta))+dist_noise_x(gen);
			p.y += v_by_yaw_rate*(cos(p.theta)-cos(p.theta+yaw_rate_mul_dt))+dist_noise_y(gen);
			p.theta += yaw_rate_mul_dt+dist_noise_theta(gen);
		}else{
			//if yaw_rate too small, motion will be CV Linear motion
			p.x += v_mul_dt*cos(p.theta)+dist_noise_x(gen);			
			p.y += v_mul_dt*sin(p.theta)+dist_noise_y(gen);
		    //FIXME : No need to update theta as TR is small?
			p.theta += yaw_rate_mul_dt + dist_noise_theta(gen);
		}
		particles[i] = p;
	}
	
	std::cout<<"prediction End"<<std::endl;
}

void ParticleFilter::dataAssociation(std::vector<Map::single_landmark_s> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	//Nearest Neighbourhood Algorithm
	
	std::cout<<"dataAssociation Start"<<std::endl;
	
	for(unsigned int i=0; i<observations.size(); i++){
		double min_dist = 0.0;
		for(unsigned int j=0; j<predicted.size(); j++){
			double cdist = dist(observations[i].x,observations[i].y,predicted[j].x_f,predicted[j].y_f);
			if(cdist < min_dist){
				min_dist = cdist;
				observations[i].id = j;//predicted[j].id_i;
			}
		}
	}
	
	std::cout<<"dataAssociation End"<<std::endl;
	
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
	
	std::cout<<"updateWeights Start"<<std::endl;
	
	for(unsigned int i=0; i<num_particles; i++){
		
		
		
		// Particle p = particles[i];
		std::vector<Map::single_landmark_s> landmarks_in_range;
		
		for(unsigned int j=0; j<map_landmarks.landmark_list.size(); j++){
			
			
			if(dist(particles[i].x,particles[i].y,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f)<=sensor_range){
				Map::single_landmark_s landmark;
				landmark.id_i = map_landmarks.landmark_list[j].id_i;
				landmark.x_f = map_landmarks.landmark_list[j].x_f;
				landmark.y_f = map_landmarks.landmark_list[j].y_f;
				
				landmarks_in_range.push_back(landmark);
			}
		}
		std::vector<LandmarkObs> obs_in_map_cord;
		
		//transforming observations to map cordinates 
		for(unsigned int j=0; j<observations.size(); j++){
			LandmarkObs obs;
			obs.x = particles[i].x + cos(particles[i].theta)*(observations[j].x) - sin(particles[i].theta)*(observations[j].y);
			obs.y = particles[i].y + sin(particles[i].theta)*(observations[j].x) + cos(particles[i].theta)*(observations[j].y);
			obs_in_map_cord.push_back(obs);
		}
		
		dataAssociation(landmarks_in_range,obs_in_map_cord);
		
		
		for(unsigned int j=0; j<obs_in_map_cord.size(); j++){
			double x = obs_in_map_cord[j].x;
			double y = obs_in_map_cord[j].y;
			double mu_x = landmarks_in_range[obs_in_map_cord[j].id].x_f;
			double mu_y = landmarks_in_range[obs_in_map_cord[j].id].y_f;
			double gauss_norm = (1.0/(2.0*M_PI*std_landmark[0]*std_landmark[1]));	
			std::cout<<"gauss_norm : "<<gauss_norm<<std::endl;
		particles[i].weight *=gauss_norm*exp(-(((0.5*(mu_x-x)*(mu_x-x))/std_landmark[0]*std_landmark[0])+((0.5*(mu_y-y)*(mu_y-y))/std_landmark[1]*std_landmark[1])));
		}
		
		std::vector<int> associations;
        std::vector<double> sense_x;
		std::vector<double> sense_y;
		
		for(unsigned int j=0; j<obs_in_map_cord.size(); j++){
			associations.push_back(obs_in_map_cord[j].id);
			sense_x.push_back(obs_in_map_cord[j].x);
			sense_y.push_back(obs_in_map_cord[j].y);
		}
		
		SetAssociations(particles[i],associations,sense_x,sense_y);
		
		std::cout<<"updateWeights End"<<std::endl;
	}
	
}

void ParticleFilter::resample() {
	
	std::cout<<"resample Start"<<std::endl;
	// DONE: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	//=============Roulette Wheel Resampling Algorithm============
	
	//Defining a RANDOM disributions
	std::default_random_engine gen;
	
	std::uniform_int_distribution<int> uni_dist_num_particle(0,num_particles-1);
	double max_weight = getBestParticle().weight;
	std::uniform_int_distribution<int> uni_dist_2wm(0,(2*max_weight)-1);
	
	//Choosing a RANDOM Index from uniform disribution
	unsigned int index = uni_dist_num_particle(gen);
	
	double beta = 0.0;
	std::vector<Particle> resampled_particles;

		
	for(unsigned int i=0; i<num_particles; i++){
		beta += uni_dist_2wm(gen);
		std::cout<<"new beta : "<<beta<<std::endl;
		while(beta>particles[index].weight){
	//		std::cout<<"Weight : "<<particles[index].weight<<std::endl;
			beta-=particles[index].weight;
	//		std::cout<<"reduced beta : "<<beta<<std::endl;
			index = (index+1)%num_particles;
	//		std::cout<<"New Index : "<<index<<std::endl;
		}
		resampled_particles.push_back(particles[index]);
	}
	particles=resampled_particles;
	std::cout<<"resample End"<<std::endl;

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
	
	//FIXED : Added Return Statement
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

/*
* returns a best particle i.e the particle having the highest/maximum weight
*/
Particle ParticleFilter::getBestParticle(){
	double max_weight = -1.0;
	Particle best_particle;
	
	for(unsigned int i=0; i<num_particles; i++){
		if(particles[i].weight>max_weight){
			max_weight = particles[i].weight;
			best_particle = particles[i];
		}
	}

	return best_particle;	
}
