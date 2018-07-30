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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	no_of_particles = 200;

	//Sensor Noise Normal Distribution...
	n_distribution<double> n_x_init(0, std[0]);
	n_distribution<double> n_y_init(0, std[1]);
	n_distribution<double> n_theta_init(0, std[2]);

	//initialize particles...
	for(int i = 0; i < no_of_particles; i++)
	{
		Particle p;
		p.id = i;
		p.x = x;
		p.y = y;
		p.theta = theta;
		p.weight = 1.0;

	//Noise...
		p.x += n_x_init(gen);
		p.y += n_y_init(gen);
		p.theta += n_theta_init(gen);

		particles.push_back(p);
	}

	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	// Sensor Noise Normal Distributions...
	normal_distribution<double> n_x(0, std_pos[0]);
	normal_distribution<double> n_y(0, std_pos[1]);
	normal_distribution<double> n_theta(0, std_pos[2]);

	for (int i = 0; i < no_of_particles; i++)
	{
		// New State Calculations...
		if(fabs(yaw_rate) < 0.00001){
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		}
		else{
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			particles[i].theta += yaw_rate * delta_t;
		}

		//Noise...
		particles[i].x += n_x(gen);
		particles[i].y += n_y(gen);
		particles[i].theta += n_theta(gen);
	}	

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


	for (unsigned int  i = 0; i <observations.size(); i++){

		// Current Observations...
		LandmarkObs Ob = observations[i];


		//Initialising Minimum to Maximum Distance Possible...
		double min_dist = numeric_limits<double>::max();

		//Initialising Landmark id...
		int map_id = -1;

		for (unsigned int j = 0; j < predicted.size(); j++)
		{
			//Current Prediction...
			LandmarkObs Op = predicted[j];

			//Distance b/w current and predicted landmark...
			double current_dist < min_dist = dist(Ob.x, Ob.y, Op.x, Op.y); 

			//Predicting the nearest landmark to current observed landmark...
			if(current_dist < min_dist){
				min_dist = current_dist;
				map_id = p.id;
			}
		}

		//Setting Observation's Id as per nearest predicted landmark's id
		observations[i].id = map_id;
	}
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


	//For every particle...
	for (int i = 0; i < no_of_particles; i++)
	{
		// Particles x and y coordinates...
		double Op_x = particles[i].x;
		double Op_y = particles[i].y;
		double Op_theta = particles[i].theta;

		//Vector to store Landmarks Location Predictions within sensor range...
		vector<LandmarkObs> predictions;

		//Map Landmarks...
		for(unsigned int j = 0; j< map_landmarks.landmark_list.size();j++){

			// To fetch x, y coordinates and id...
			float lml_x = map_landmarks.landmark_list[j].x_f;
			float lml_y = map_landmarks.landmark_list[j].y_f;
			int lml_id = map_landmarks.landmark_list[j].id_i;

			//To Conider the Landmarks which are only within sensor range of particle...
			if(fabs(lml_x - Op_x) <= sensor_range && fabs(lml_y - Op_y) <= sensor_range){

				//Adding prediction to vector...
				predictions.push_back(LandmarkObs{lml_id, lml_x, lml_y});
			}

	}

	//Creating the copy of list of observations transformed from vehicle coordinates to map coordinates...
	vector<LandmarkObs> transformed_os;
	for (unsigned int j = 0; j observations.size(); j++){
		double Ot_x = cos(Op_theta)*observations[j].x - sin(Op_theta)*observations[j].y + p_x;
		double Ot_y = sin(Op_theta)*observations[j].x + cos(Op_theta)*observations[j].y + p_y;
		transformed_os.push_back(LandmarkObs{ observations[j].id, Ot_x, Ot_y });
		}

	//To Perform dataAssociation for the predicted and transformed observations of the current particles...
	dataAssociation(predictions, transformed_os);

	//Reinitialize weight...
	particles[i].weight = 1.0;

	for (unsigned int j = 0; j < count; j++)
	{
		//Placeholders for Observations and Associated Predicted Coordinates...
		double Oo_x, Oo_y, Opr_x, Opr_y;
		Oo_x = transformed_os[j].x;
		Oo_y = transformed_os[j].y;

		int associated_prediction = transformed_os[j].id;

		//To fetch x and y coordinates of prediction of current observation
		for (unsigned int k = 0; k < predictions.size(); k++){
			if (predictions[k].id == associated_prediction){
				Opr_x = predictions[k].x;
				Opr_y = predictions[k].y;
			}
		}

	    //Calculation of weight of current observation using Multivariate Guassian
	    double Os_x = std_landmark[0];
	    double Os_y = std_landmark[1];
	    double obs_w = ( 1/(2*M_PI*Os_x*Os_y)) * exp( -( pow(Opr_x - Oo_x,2)/(2*pow(Os_x, 2)) + (pow(Opr_y - Oo_y , 2)/(2*pow(Os_y, 2)))));

	    //Current Observation Weight * Total Observation Weight...
	    particles[i].weight *=obs_w;
	}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
