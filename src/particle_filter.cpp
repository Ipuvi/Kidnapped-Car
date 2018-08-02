#include <iostream>
#include <sstream>
#include <algorithm>
#include <random>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static int Num_Particles = 100;

// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of x, y, theta and their uncertainties from GPS) and all weights to 1. 
void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// Add random Gaussian noise to each particle.

	//Normal Distribution for theta, x and y...
	std::normal_distribution<double> n_x_(x, std[0]);
	std::normal_distribution<double> n_y_(y, std[1]);
	std::normal_distribution<double> n_theta_(theta, std[2]);
	std::default_random_engine gen;
	
	//Resize weight and vector of particles
	num_particles = Num_Particles;
	particles.resize(num_particles);

	//initialize particles...
	for(auto& Op: particles)
	{
	   Op.x = n_x_(gen);
	   Op.y = n_y_(gen);
	   Op.theta = n_theta_(gen);
	   Op.weight = 1;
	}

	is_initialized = true;

}


// TODO: Add measurements to each particle and add random Gaussian noise.

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	std::default_random_engine gen;

	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//Create Normal Guassian Noise Distributions...
	
	std::normal_distribution<double> N_x(0, std_pos[0]);
	std::normal_distribution<double> N_y(0, std_pos[1]);
	std::normal_distribution<double> N_theta(0, std_pos[2]);
	


	for (auto& Op: particles)
	{
	  // New State Calculations...
	  if( fabs(yaw_rate) < 0.0001)
 	   {
	    Op.x += velocity * delta_t * cos(Op.theta);
	    Op.y += velocity * delta_t * sin(Op.theta);
		}
	  else
    	   {
	    Op.x += velocity / yaw_rate * (sin( Op.theta + yaw_rate*delta_t) - sin(Op.theta) );
	    Op.y += velocity / yaw_rate * (cos( Op.theta) - cos(Op.theta + yaw_rate*delta_t ) );
	    Op.theta += yaw_rate * delta_t;
	        }
		
	  //Predicted Particles...
	  //Sensor Noise...
	    Op.x += N_x(gen);
	    Op.y += N_y(gen);
	    Op.theta += N_theta(gen);
	 }	

}

// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	// Current Observations...
	for (auto &Obs_o: observations){

	//Initialising Minimum to Maximum Distance Possible...
	double min_dist = std::numeric_limits<float>::max();

	//Current Prediction...
	for (const auto &Pred_o: predicted)
	{
	//Distance b/w current and predicted landmark...
	double current_dist = dist(Obs_o.x, Obs_o.y, Pred_o.x, Pred_o.y); 

	//Predicting the nearest landmark to current observed landmark...
	if(current_dist < min_dist){
	  min_dist = current_dist;
	//Setting Observation's Id as per nearest predicted landmark's id
	  Obs_o.id = Pred_o.id;
	}
      }

   }
}


// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], std::vector<LandmarkObs> observations, Map map_landmarks) 
{
// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located according to the MAP'S coordinate system. You will need to transform between the two systems.

	//For every particle...
	for (auto& Op: particles)
	{	
	//Update Weight...
	   Op.weight = 1.0;
		
	//Vector to store Landmarks Location Predictions within sensor range...
	vector<LandmarkObs> predictions;

	//Map Landmarks...
	for(const auto &lml: map_landmarks.landmark_list)
	{
	  double distance = dist(Op.x, Op.y, lml.x_f, lml.y_f);
		
	  //To Consider the Landmarks which are only within sensor range of particle...
	  if(distance < sensor_range){
	     //Adding prediction to vector...
	     predictions.push_back(LandmarkObs{lml.id_i, lml.x_f, lml.y_f});
	}
      }


//Creating the copy of list of observations transformed from vehicle coordinates to map coordinates...
	vector<LandmarkObs> transformed_obs;
	double cos_value = cos(Op.theta);
	double sin_value = sin(Op.theta);

	for (const auto &Obs_o: observations){
	 LandmarkObs temp;
	 temp.x = cos_value * Obs_o.x - sin(Op.theta) * Obs_o.y + Op.x;
	 temp.y = sin_value * Obs_o.x + cos(Op.theta) * Obs_o.y + Op.y;
	 transformed_obs.push_back(temp);
		}

	//To Perform dataAssociation for the predicted and transformed observations of the current particles...
	dataAssociation(predictions, transformed_obs);

	//Compute  Weight...
	for (const auto& trans_obs: transformed_obs)
	{
	 //Placeholders for Observations and Associated Predicted Coordinates...
	 Map::single_landmark_s landmark = map_landmarks.landmark_list.at(trans_obs.id-1);
	 double trans_obs_x = pow(trans_obs.x - landmark.x_f, 2) / (2 * pow(std_landmark[0], 2));
	 double trans_obs_y = pow(trans_obs.y - landmark.y_f, 2) / (2 * pow(std_landmark[1], 2));
		
  	 //Calculation of weight of current observation using Multivariate Guassian
	 double w = exp(-(trans_obs_x + trans_obs_y)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);

	 //Current Observation Weight * Total Observation Weight...
	 Op.weight *= w;

	}
  	weights.push_back(Op.weight);
       }
}


// TODO: Resample particles with replacement with probability proportional to their weight. 
void ParticleFilter::resample() {
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> dist(weights.begin(), weights.end());

	//
	vector<Particle> new_particles;
	new_particles.resize(num_particles);

	//Fetching Current Weights...

	//
	for (int i = 0; i < num_particles; i++){
	 int ids = dist(gen);
	 new_particles[i] = particles[ids];
	}
	
	//
	particles = new_particles;

	//
	weights.clear();
	
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
     //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
     // sense_x: the associations x mapping already converted to world coordinates
     // associations: The landmark id that goes along with each listed association
     // sense_y: the associations y mapping already converted to world coordinates


    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

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
