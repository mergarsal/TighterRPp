#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <chrono>  // timer
#include <sstream>
#include <fstream>  // for the file

#include <math.h>  // rnd, acos


#include "PrimalTypes.h"
#include "PrimalUtils.h"
#include "PrimalLeft.h"
#include "PrimalRight.h"
#include "PrimalBoth.h"
#include "PrimalAdjugate.h"

#include "./../test/generateCorrespondences.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

#define PI 3.1415916


double computeCost(Eigen::Matrix<double, 9, 9> & C , Eigen::Matrix<double, 3, 3> & sol)
{        
        Eigen::Matrix<double, 9, 1> e = Eigen::Map<Eigen::Matrix<double, 9, 1>>(sol.data(), 9, 1);
        return (e.transpose() * C * e);
};


double distRot(Eigen::Matrix<double, 3, 3> & Rh, Eigen::Matrix<double, 3, 3> & R)
{
        double param = ((Rh.transpose() * R).trace() - 1) * 0.5;
        return (acos (param) * 180.0 / PI);
};


double distTrans(Eigen::Matrix<double, 3, 1> & th, Eigen::Matrix<double, 3, 1> & t)
{
        Eigen::Vector3d th_norm = th.normalized();
        Eigen::Vector3d t_norm = t.normalized(); 
        
        double e_pos = (acos (th.dot(t)) * 180.0 / PI); 
        double e_neg = (acos (-th.dot(t)) * 180.0 / PI);
        return (std::min(e_pos, e_neg)); 
};




int main(int argc, char** argv)
{

        std::srand(std::time(nullptr));

        /* Create points */

        UtilsRelPose::Vector3 translation;
        UtilsRelPose::Matrix3 rotation;
        UtilsRelPose::bearing_vectors_t points_correspondences;


        int number_corr = 100; 
        
        // define struct with params
        UtilsRelPose::StrInputGenProblem str_in = UtilsRelPose::StrInputGenProblem(); 
        str_in.FoV = 100; 
        str_in.min_depth = 2; 
        str_in.max_depth = 8; 
        str_in.focal_length = 800; 
        str_in.n_points = number_corr; 
        str_in.noise = 0.5; 
        str_in.parallax = 2.0; 
        str_in.outlier_fraction = 0;
        str_in.max_rotation = 0.5;


        // generate problem
        UtilsRelPose::StrOutputGenProblem str_out = UtilsRelPose::createSyntheticExperiment(str_in);  
        // extract data

        translation = str_out.translation; 
        rotation = str_out.rotation; 




        for (int i=0; i < number_corr; i++)
        {
                UtilsRelPose::CorrespondingFeatures correspondence = str_out.points_correspondences[i];
                points_correspondences.push_back(correspondence);                                                                                             
        }      




        UtilsRelPose::Matrix9 C = UtilsRelPose::constructDataMatrix(points_correspondences);
       
        // n_point
        UtilsRelPose::Matrix3 E = UtilsRelPose::computeEfromRt(rotation, translation);
        UtilsRelPose::Vector3 q; 
        translation.normalize(); 


        /*  ------------------------------------------------  */
        /* Run SDP estimation */                   


        /*  LEFT  */
        // 1. create problem

        auto start_left = high_resolution_clock::now();                                                                                  
        RelPose::PrimalLeft sdp_left =  RelPose::PrimalLeft(C, false);
                                                
        // 2. Run solver 
        RelPose::PrimalLeftResult res_left = sdp_left.getResult(false);
                                              
        // 3. Obtain relative pose from SDP solution 
        sdp_left.getSolutionFromResult(res_left, false); 
        auto duration_left = duration_cast<microseconds>(high_resolution_clock::now() - start_left);                                             



        // sdp_left -> printResult(res_left);                                              
        UtilsRelPose::Matrix3 Eleft = res_left.E_opt;                                              
        double fleft = computeCost(C, Eleft);
        UtilsRelPose::Matrix3 Rleft; 
        UtilsRelPose::Vector3 tleft; 
        UtilsRelPose::computeRtfromE(points_correspondences, Eleft, Rleft, tleft);
        double rot_left = distRot(Rleft, rotation); 
        double trans_left = distTrans(tleft, translation); 



        /*  RIGHT  */
        // 1. create problem 

        auto start_right = high_resolution_clock::now();                                                                                                                            
        RelPose::PrimalRight sdp_right = RelPose::PrimalRight(C, false);
                                                
        // 2. Run solver 
        RelPose::PrimalRightResult res_right = sdp_right.getResult(false);
                                              
        // 3. Obtain relative pose from SDP solution 
        sdp_right.getSolutionFromResult(res_right, false);                                              
        auto duration_right = duration_cast<microseconds>(high_resolution_clock::now() - start_right);                                 


        // sdp_left -> printResult(res_left);                                              
        UtilsRelPose::Matrix3 Eright = res_right.E_opt; 
        double fright = computeCost(C, Eright);
        UtilsRelPose::Matrix3 Rright; 
        UtilsRelPose::Vector3 tright; 
        UtilsRelPose::computeRtfromE(points_correspondences, Eright, Rright, tright);
        double rot_right = distRot(Rright, rotation); 
        double trans_right = distTrans(tright, translation); 


        /*  BOTH  */
        // 1. create problem


        auto start_both = high_resolution_clock::now();                                                                                                                             
        RelPose::PrimalBoth sdp_both = RelPose::PrimalBoth(C, false);
                                                
        // 2. Run solver 
        RelPose::PrimalBothResult res_both = sdp_both.getResult(false);
                                             
        // 3. Obtain relative pose from SDP solution 
        sdp_both.getSolutionFromResult(res_both, false);                                              
        auto duration_both = duration_cast<microseconds>(high_resolution_clock::now() - start_both);                                  



        // sdp_left -> printResult(res_left);                                              
        UtilsRelPose::Matrix3 Eboth = res_both.E_opt; 
        double fboth = computeCost(C, Eboth);
        UtilsRelPose::Matrix3 Rboth; 
        UtilsRelPose::Vector3 tboth; 
        UtilsRelPose::computeRtfromE(points_correspondences, Eboth, Rboth, tboth);
        double rot_both = distRot(Rboth, rotation); 
        double trans_both = distTrans(tboth, translation); 


        /*  ADJUGATE  */
        // 1. create problem

        auto start_adj = high_resolution_clock::now();                                                                                                                             
        RelPose::PrimalAdjugate sdp_adj = RelPose::PrimalAdjugate(C, false);
                                                
        // 2. Run solver 
        RelPose::PrimalAdjugateResult res_adj = sdp_adj.getResult(false);
                                              
        // 3. Obtain relative pose from SDP solution 
        sdp_adj.getSolutionFromResult(res_adj, false);                                              
        auto duration_adj = duration_cast<microseconds>(high_resolution_clock::now() - start_adj);                                 


        // sdp_left -> printResult(res_left);                                              
        UtilsRelPose::Matrix3 Eadj = res_adj.E_opt; 
        double fadj = computeCost(C, Eadj);
        UtilsRelPose::Matrix3 Radj; 
        UtilsRelPose::Vector3 tadj; 
        UtilsRelPose::computeRtfromE(points_correspondences, Eadj, Radj, tadj);
        double rot_adj = distRot(Radj, rotation); 
        double trans_adj = distTrans(tadj, translation); 


                                       
        std::cout << "Error in rotation (degrees):\n"; 
        std::cout << "Left: " << rot_left << std::endl; 
        std::cout << "Right: " << rot_right << std::endl; 
        std::cout << "Both: " << rot_both << std::endl; 
        std::cout << "Adjugate: " << rot_adj << std::endl; 


        std::cout << "Error in translation (degrees):\n"; 
        std::cout << "Left: " << trans_left << std::endl; 
        std::cout << "Right: " << trans_right << std::endl; 
        std::cout << "Both: " << trans_both << std::endl; 
        std::cout << "Adjugate: " << trans_adj << std::endl; 

        std::cout << "Cost for solution:\n"; 
        std::cout << "Left: " << fleft << std::endl; 
        std::cout << "Right: " << fright << std::endl; 
        std::cout << "Both: " << fboth << std::endl; 
        std::cout << "Adjugate: " << fadj << std::endl; 
        
        std::cout << "Time for solvers (microseconds):\n"; 
        std::cout << "Left: " << (double) duration_left.count() << std::endl; 
        std::cout << "Right: " << (double) duration_right.count() << std::endl; 
        std::cout << "Both: " << (double) duration_both.count() << std::endl; 
        std::cout << "Adjugate: " << (double) duration_adj.count() << std::endl; 


  return 0;

}




