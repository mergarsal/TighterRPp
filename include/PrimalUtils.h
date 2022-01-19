#pragma once

#include <eigen3/Eigen/Dense>


#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

#include "PrimalTypes.h"


namespace RelPose
{
        
        template <typename Vector, typename MatrixT, typename VectorT>
        struct PrimalResult
        {
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                // fill by getResult
                double stop_iter; 
                double f_opt; 
                double d_opt; 
                double primal_error;
                double dual_error; 
        
        
                Vector dual_point; 
                
                Matrix9 primal_Xe;
                MatrixT primal_Xt;   
                
                // modified by getSolutionFromResult
                Matrix3 E_opt; 
                VectorT t_opt;       
                
                int rank_Xe; 
                std::vector<int> rank_Xt;            
                        
                double stable_rank_Xe; 
                std::vector<double> stable_rank_Xt; 
                
                Vector9 eigenvalues_Xe; 
                VectorT eigenvalues_Xt;  
                
        
        };  // end of PrimalResult

        
        template <typename Vector, typename MatrixT, typename VectorT>
        class PrimalClassBase
        {
                public:
                        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                        
                        PrimalClassBase(void) {}; 
                        
                        PrimalClassBase(const Matrix9 & C, const bool debug = false){}; 
                        
                        ~PrimalClassBase(void) {  };
                             
                        virtual PrimalResult<Vector, MatrixT, VectorT> getResult(bool debug = false){};
                                                              
                        virtual void getSolutionFromResult(PrimalResult<Vector, MatrixT, VectorT> & res, bool debug = false){}; 
                                                              
                        void printResult(const PrimalResult<Vector, MatrixT, VectorT> & result); 
        
                        
        
        
        };  // end of Class PrimalClassBase

        int computeRankX(const Eigen::Matrix<double, Eigen::Dynamic, 1> & e_vector,double max_eig); 
         
        double computeStableRankX(const Eigen::MatrixXd & X_matrix, 
                                                                double max_eig);
                                                                
                                                                
}  // end of namespace PrimalRotPrior
