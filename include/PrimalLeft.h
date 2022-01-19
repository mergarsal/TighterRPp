#pragma once

#include <eigen3/Eigen/Dense>

#include "PrimalUtils.h"

#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

#include "PrimalTypes.h"



namespace RelPose
{        

        typedef PrimalResult<Vector7, Matrix3, Vector3> PrimalLeftResult;
        
           
        class PrimalLeft : public PrimalClassBase<Vector7, Matrix3, Vector3>
        {
        
        public: 
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                
                PrimalLeft(void) {};
                
                PrimalLeft(const Matrix9 & C, const bool debug = false);  
                                                
                ~PrimalLeft(void) { Problem1_.terminate();  }
                             
                PrimalLeftResult getResult(bool debug = false);
                                                      
                void getSolutionFromResult(PrimalLeftResult & res, bool debug = false); 
                
                
        private:
        
                       SDPA Problem1_; 
                       int dimE_; 
                       int dimT_; 
                       int nConstraints_;
                       bool debug_; 
                       Matrix9 C_; 
                                                      
        };  // end of class left
               

        
}  // end of namespace PrimalRotPrior
