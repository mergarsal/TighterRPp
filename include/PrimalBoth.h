#pragma once

#include <eigen3/Eigen/Dense>

#include "PrimalUtils.h"

#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

#include "PrimalTypes.h"



namespace RelPose
{        

        typedef PrimalResult<Vector13, Matrix6, Vector6> PrimalBothResult;
        
           
        class PrimalBoth : public PrimalClassBase<Vector13, Matrix6, Vector6>
        {
        
        public: 
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                
                PrimalBoth(void) {};
                
                PrimalBoth(const Matrix9 & C, const bool debug = false);  
                                
                ~PrimalBoth(void) { Problem1_.terminate(); }
                             
                PrimalBothResult getResult(bool debug = false);
                                                      
                void getSolutionFromResult(PrimalBothResult & res, bool debug = false); 
                                                      
                                                      
                                                      
        private: 
                       int dimQ_; 
                       SDPA Problem1_; 
                       int dimE_; 
                       int dimT_; 
                       int nConstraints_;
                       bool debug_; 
                       Matrix9 C_; 
        };  // end of class left
               

        
}  // end of namespace PrimalRotPrior
