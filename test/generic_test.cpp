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

// include helper 
#include "experimentsHelper.h"
#include "generateCorrespondences.h"

using namespace std;
using namespace Eigen;
using namespace std::chrono;

#define PI 3.1415916

/* 
SAVE: 
1. file with cost
2. file with error in rotation
3. file with error in translation
4. file with times for the certifiers 
*/


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
   /* Read params from input */ 
   
   string name_in_params = "basic_params.txt"; 
   
   /* Read the name of the file */
   if (argc > 1)
        name_in_params = argv[1]; 
        
  
  SceneOptions options; 
   
  std::cout << "Test for prior\n"; 
  std::cout << "Input for the test: " << name_in_params << std::endl; 
  
  
  // read params from file
  bool valid_options = readOptionsFromFile(name_in_params, options);   
  

  std::srand(std::time(nullptr));
  
  
    for (size_t noise_id=0; noise_id < options.n_noise; noise_id++)
    {
        double noise_i = options.noise[noise_id];

        for (size_t fov_id = 0; fov_id < options.n_fov; fov_id++)
        {
                double fov_i = options.FoV[fov_id];
        
        
                for (size_t par_id = 0; par_id < options.n_parallax; par_id++)
                {
                        double par_i = options.max_parallax[par_id];    
        
        
                        for (size_t focal_id = 0; focal_id < options.n_focal; focal_id++)
                        {
                                double focal_i = options.focal_length[focal_id];
                                
                               for (size_t nr_id = 0; nr_id < options.n_corr; nr_id++)
                               {
                                        double number_corr = options.number_correspondences[nr_id]; 
                                         
                                         // This file saves all our resutls
                                         // rank of the returned solution for SDP
                                         auto name_file_cost = "cost_SDP_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream fcost(name_file_cost);
                                         
                                         auto name_file_primal = "primal_SDP_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream fprimal(name_file_primal);
                                         
                                         
                                         auto name_file_dual = "dual_SDP_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream fdual(name_file_dual);
                                         
                                         auto name_file_gap = "gap_SDP_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream fgap(name_file_gap);
                                         
                                         // is opt for each method
                                         auto name_file_rot = "rot_SDP_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream frot(name_file_rot);

                                         // min. eigenvalue of Hessian for each certifier
                                         auto name_file_trans = "trans_SDP_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream ftrans(name_file_trans);
                                                                                  
                                         // time 
                                         auto name_file_time = "time_" + std::to_string(noise_i) 
                                                                     + "_fov_" + std::to_string((int)fov_i) 
                                                                     + "_par_" + std::to_string((int)par_i) 
                                                                     + "_focal_" + std::to_string((int)focal_i) 
                                                                     + "_N_" + std::to_string((int)number_corr) 
                                                                     + ".txt";
                                         std::ofstream ftimes(name_file_time);
                                         
                              
         
        
                                       for (size_t n_iter = 0; n_iter < options.max_iter; n_iter++)
                                       {

                                               /* Create points */
                                  
                                               std::cout << "Iter: " << n_iter << std::endl; 
                                               
                                               UtilsRelPose::Vector3 translation;
                                               UtilsRelPose::Matrix3 rotation;
                                               UtilsRelPose::bearing_vectors_t points_correspondences;
                                                
                                               
                                               
                                               // define struct with params
                                               UtilsRelPose::StrInputGenProblem str_in = UtilsRelPose::StrInputGenProblem(); 
                                               str_in.FoV = fov_i; 
                                               str_in.min_depth = options.min_depth; 
                                               str_in.max_depth = options.max_depth; 
                                               str_in.focal_length = focal_i; 
                                               str_in.n_points = number_corr; 
                                               str_in.noise = noise_i; 
                                               str_in.parallax = par_i; 
                                               str_in.outlier_fraction = options.outlier_fraction[0];
                                               str_in.max_rotation = options.max_rotation;
                                               
                                               
                                               str_in.custom_t = options.custom_t; 
                                               
                                               
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
                                               // normalize C 
                                               // C /= points_correspondences.size(); 
                                              
                                                double s = 0;
                                                double s2 = 0;
                                                for (int i = 0; i < 9; i++)
                                                {
                                                    for (int j = 0; j < 9; j++)
                                                    {
                                                            double tmp = C(i, j);
                                                            s += tmp;
                                                            s2 += tmp*tmp;
                                                    }
                                                }
                                                double EX = s / 81;
                                                double EX2 = s2 / 81;
                                                double c_std = sqrt(EX2 - EX*EX);
                                                for (int i = 0; i < 9; i++)
                                                {
                                                    for (int j = 0; j < 9; j++)
                                                    {
                                                        C(i, j) /= c_std;
                                                    }    
                                                }
                                                
       
       
       
                  
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
                                             
                                         
                                                                               
                                             double c = C.trace(); 
                              
      
                                              /*  ------------------------------------------------  */
                                              // Save results
                                            
                                              /* Costs */ 
                                              
                                              fcost << fleft / c<< ","; 
                                              fcost << fright / c<< ","; 
                                              fcost << fboth / c<< ","; 
                                              fcost << fadj/ c;
                                              fcost << std::endl; 
                                              
                                              /* Primal */
                                              
                                              fprimal << res_left.f_opt / c << ","; 
                                              fprimal << res_right.f_opt / c << ","; 
                                              fprimal << res_both.f_opt / c << ","; 
                                              fprimal << res_adj.f_opt/ c;
                                              fprimal << std::endl; 
                                              
                                              /* Dual */
                                              
                                              fdual << res_left.d_opt / c<< ","; 
                                              fdual << res_right.d_opt / c<< ","; 
                                              fdual << res_both.d_opt / c<< ","; 
                                              fdual << res_adj.d_opt/ c;
                                              fdual << std::endl; 
                                              
                                              /* Dual gap */
                                              
                                              fgap << fleft / c  - res_left.f_opt / c << ","; 
                                              fgap << fright / c - res_right.f_opt / c << ","; 
                                              fgap << fboth / c  - res_both.f_opt / c << ","; 
                                              fgap << fadj / c   - res_adj.f_opt/ c;
                                              fgap << std::endl; 
                                              
                                              
                                              /* Error rot */
                                              
                                              frot << rot_left << ","; 
                                              frot << rot_right << ","; 
                                              frot << rot_both << ","; 
                                              frot << rot_adj << ","; 
                                              frot << std::endl;
                                              
                                              
                                              /* Error trans */
                                              
                                              ftrans << trans_left << ","; 
                                              ftrans << trans_right << ","; 
                                              ftrans << trans_both << ","; 
                                              ftrans << trans_adj << ","; 
                                              ftrans << std::endl;
                                              
                                              
                                              /* Time required by each certifier */
                                              
                                              ftimes << (double)duration_left.count()  << ","; 
                                              ftimes << (double)duration_right.count() << ",";  
                                              ftimes << (double)duration_both.count() << ",";  
                                              ftimes << (double)duration_adj.count() << ",";  
                                              ftimes << std::endl; 
                                              
                                              
                                               
                                        }  // end for each iter
                                      
                                      
                                      
                                      // close file with cost
                                      fcost.close();
                                      // close file with error in rotation
                                      frot.close();
                                      // close file with error in translation
                                      ftrans.close();
                                      // close file with index of relaxation 
                                      ftimes.close(); 
                                      
                                }  // end for number correspondences  
      
                         }  // enf for focal
       
                 }  // end for parallax
      
        }  // end for fov 
      
      }  // end for noise

  return 0;

}




