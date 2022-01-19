#pragma once 

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>  // for the file

#include <eigen3/Eigen/Dense>

using namespace std;


struct SceneOptions {
        double max_iter = 100;  // # Number of iterations 
        
        int n_fov = 1;   // number of elements in FoV
        vector<double> FoV = {100};           // # FoV

        int n_parallax = 1; 
        vector<double> max_parallax = {2};  // # max parallax
        
        double min_depth = 1.0; 
        
        double max_depth = 8.; 
        
        int n_outliers = 0; 
        vector<double> outlier_fraction;
        
        int n_focal = 1; 
        vector<double> focal_length = {800}; 
        
        int n_noise = 1; 
        vector<double> noise = {0.5}; 
        
        int n_corr = 1; 
        vector<int> number_correspondences = {100}; 
        
        bool use_custom_t = false; 
        
        Eigen::Matrix<double, 3, 1> custom_t; 
        
        double max_rotation = 0.5; 
        
        int method_init = 1;
        
        bool use_all_corr = false;
        
        double max_rotation_pert = 0.0; 
        
};


template <typename Scalar=double>
void readSetParams(ifstream & in_file, int & N, vector<Scalar> & arr);




bool readOptionsFromFile(string name_file, 
                         SceneOptions & options);



/*double computeCost(Eigen::Matrix<double, 9, 9> & C , Eigen::Matrix<double, 3, 3> & sol)
{        
        Eigen::Matrix<double, 9, 1> e = Eigen::Map<Eigen::Matrix<double, 9, 1>>(sol.data(), 9, 1);
        // e.normalize();
        return (e.transpose() * C * e);
};
*/


