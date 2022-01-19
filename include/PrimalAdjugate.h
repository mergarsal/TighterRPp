#pragma once

#include <eigen3/Eigen/Dense>

#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

#include "PrimalTypes.h"
#include "PrimalUtils.h"


namespace RelPose
{        

        typedef PrimalResult<Vector22, Matrix6, Vector6> PrimalAdjugateResult;
        
           
        class PrimalAdjugate : public PrimalClassBase<Vector22, Matrix6, Vector6>
        {
        
        public: 
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW
                
                PrimalAdjugate(void) {};
                
                PrimalAdjugate(const Matrix9 & C, const bool debug = false);  
                                
                ~PrimalAdjugate(void) { Problem1_.terminate();  }
                             
                PrimalAdjugateResult getResult(bool debug = false);
                                                      
                void getSolutionFromResult(PrimalAdjugateResult & res, bool debug = false); 
                                                      
        private:
        
                       SDPA Problem1_; 
                       int dimE_; 
                       int dimT_; 
                       int nConstraints_;
                       bool debug_; 
                       Matrix9 C_;                                               
        };  // end of class left
               

        
}  // end of namespace PrimalRotPrior
