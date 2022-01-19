
#include <cstdio>
#include <cstdlib>
#include <sdpa_call.h>

// for eigendecomposition
#include <eigen3/Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include "PrimalTypes.h"
#include "PrimalUtils.h"
#include "PrimalAdjugate.h"


namespace RelPose
{


PrimalAdjugate::PrimalAdjugate(const Matrix9& C, const bool debug)
{ 
 
     
  Problem1_ = SDPA(); 
  
  
  dimE_ = 9; 
  dimT_ = 6;
  nConstraints_ = 6 + 6 + 9 + 1;  // NOTE: we remove the last one. See note below
  // 6 for E E^T and 6 for E^T E, 9 - 1 adjugate, 1 norm E 
  int nBlock = 2;    // number blocks: E; tq
  
  
  // display info from solver
  if (debug) Problem1_.setDisplay(stdout);
  
  
  // All parameteres are renewed
  Problem1_.setParameterType(SDPA::PARAMETER_DEFAULT);
  // Problem1_.printParameters(stdout);

  if (debug) std::cout << "[SDP-ADJ] Create type solution\n";
  
  Problem1_.inputConstraintNumber(nConstraints_);
  Problem1_.inputBlockNumber(nBlock);
  // bLOCKsTRUCT :: 9(SDP)  6(SDP) 
  Problem1_.inputBlockSize(1,dimE_);
  Problem1_.inputBlockSize(2,dimT_);
  Problem1_.inputBlockType(1,SDPA::SDP);
  Problem1_.inputBlockType(2,SDPA::SDP);

  Problem1_.initializeUpperTriangleSpace();
  
  if (debug) std::cout << "[SDP-ADJ] Define vector b\n";
  //cVECT = {1,1, ..., 0}
  Problem1_.inputCVec(1,1);
  Problem1_.inputCVec(2,1);
  
  // NOTE: needed??
  for (int i=3; i <= nConstraints_; i++) Problem1_.inputCVec(i,0); 


  // ---------  Input F_0  (cost) --------------------

  if (debug) std::cout << "[SDP-ADJ] Define data matrix\n";
  // 1st block: C
  for (int i = 1; i <= 9; i++)
  {
        for (int j = i; j <= 9; j++)
        {
                 // NOTE: SDPA consider the problem as MAX
                 Problem1_.inputElement(0, 1, i, j, - C(i - 1,j - 1)); 
        }
        
  }
  
  

  /* Define constraints */
  
  // auxiliary variables
  size_t e1=1, e4=2, e7=3, e2=4, e5=5, e8=6, e3=7, e6=8, e9=9;
  // Recall that we have THREE blocks
  int t1=1,t2=2,t3=3,q1=4, q2=5, q3=6;
  
  // Norm t
  if (debug) std::cout << "[SDP-ADJ] Define constraint norm t\n";
  // 2nd block
  // # constraint, block, i1, j1, value
  Problem1_.inputElement(1, 2, t1, t1, 1);
  Problem1_.inputElement(1, 2, t2, t2, 1);
  Problem1_.inputElement(1, 2, t3, t3, 1); 
 
  // 3rd block
  Problem1_.inputElement(2, 2, q1, q1, 1);
  Problem1_.inputElement(2, 2, q2, q2, 1);
  Problem1_.inputElement(2, 2, q3, q3, 1); 
  
  
  if (debug) std::cout << "[SDP-ADJ] Define constraint E*E^T\n";
  // ---------  Input F_2  --------------------

  
  /* Constraint for essential set */
  size_t id_const = 2; 
  // # constraint, block, i1, j1, value 
  /** LEFT **/
  /* First set */ 
  /*
  id_const++;
  Problem1_.inputElement(id_const, 1, e2, e2, 1);
  Problem1_.inputElement(id_const, 1, e1, e1, 1);
  Problem1_.inputElement(id_const, 1, e6, e6, 1);
  Problem1_.inputElement(id_const, 2, t1, t1, -1);
  Problem1_.inputElement(id_const, 2, t3, t3, -1);
  */
  id_const++;      
  Problem1_.inputElement(id_const, 1, e4, e4, 1);
  Problem1_.inputElement(id_const, 1, e5, e5, 1);
  Problem1_.inputElement(id_const, 1, e6, e6, 1);
  Problem1_.inputElement(id_const, 2, t1, t1, -1);
  Problem1_.inputElement(id_const, 2, t3, t3, -1);

  id_const++;
  Problem1_.inputElement(id_const, 1, e7, e7, 1);
  Problem1_.inputElement(id_const, 1, e8, e8, 1);
  Problem1_.inputElement(id_const, 1, e9, e9, 1);
  Problem1_.inputElement(id_const, 2, t1, t1, -1);
  Problem1_.inputElement(id_const, 2, t2, t2, -1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e4, 1);
  Problem1_.inputElement(id_const, 1, e2, e5, 1);
  Problem1_.inputElement(id_const, 1, e3, e6, 1);
  Problem1_.inputElement(id_const, 2, t1, t2, 1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e7, 1);
  Problem1_.inputElement(id_const, 1, e2, e8, 1);
  Problem1_.inputElement(id_const, 1, e3, e9, 1);
  Problem1_.inputElement(id_const, 2, t1, t3, 1);

  id_const++;
  Problem1_.inputElement(id_const, 1, e4, e7, 1);
  Problem1_.inputElement(id_const, 1, e5, e8, 1);
  Problem1_.inputElement(id_const, 1, e6, e9, 1);
  Problem1_.inputElement(id_const, 2, t2, t3, 1);
  

  /** RIGHT **/
  /* Go for the other formulation */
  /* First set */ 
  /*
  id_const++;
  Problem1_.inputElement(id_const, 1, e3, e3, 1);
  Problem1_.inputElement(id_const, 1, e6, e6, 1);
  Problem1_.inputElement(id_const, 1, e9, e9, 1);
  Problem1_.inputElement(id_const, 2, q1, q1, -1);
  Problem1_.inputElement(id_const, 2, q2, q2, -1);
  */
        
        
  id_const++;
  Problem1_.inputElement(id_const, 1, e2, e2, 1);
  Problem1_.inputElement(id_const, 1, e5, e5, 1);
  Problem1_.inputElement(id_const, 1, e8, e8, 1);
  Problem1_.inputElement(id_const, 2, q1, q1, -1);
  Problem1_.inputElement(id_const, 2, q3, q3, -1);
  
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e1, 1);
  Problem1_.inputElement(id_const, 1, e4, e4, 1);
  Problem1_.inputElement(id_const, 1, e7, e7, 1);
  Problem1_.inputElement(id_const, 2, q2, q2, -1);
  Problem1_.inputElement(id_const, 2, q3, q3, -1);
  
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e2, 1);
  Problem1_.inputElement(id_const, 1, e4, e5, 1);
  Problem1_.inputElement(id_const, 1, e7, e8, 1);
  Problem1_.inputElement(id_const, 2, q1, q2, 1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e3, 1);
  Problem1_.inputElement(id_const, 1, e4, e6, 1);
  Problem1_.inputElement(id_const, 1, e7, e9, 1);
  Problem1_.inputElement(id_const, 2, q1, q3, 1);

  id_const++;
  Problem1_.inputElement(id_const, 1, e2, e3, 1);
  Problem1_.inputElement(id_const, 1, e5, e6, 1);
  Problem1_.inputElement(id_const, 1, e8, e9, 1);
  Problem1_.inputElement(id_const, 2, q2, q3, 1);
  
  
  /** ADJUGATE **/
  /* Norm E constraint */
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e1, 1);
  Problem1_.inputElement(id_const, 1, e2, e2, 1);
  Problem1_.inputElement(id_const, 1, e3, e3, 1);
  Problem1_.inputElement(id_const, 1, e4, e4, 1);
  Problem1_.inputElement(id_const, 1, e5, e5, 1);
  Problem1_.inputElement(id_const, 1, e6, e6, 1);
  Problem1_.inputElement(id_const, 1, e7, e7, 1);
  Problem1_.inputElement(id_const, 1, e8, e8, 1);
  Problem1_.inputElement(id_const, 1, e9, e9, 1);
  // Update entry in vector b
  Problem1_.inputCVec(id_const,2); 
  
  /* Adjugate constraints */
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e5, e9, 1);
  Problem1_.inputElement(id_const, 1, e6, e8, -1);
  Problem1_.inputElement(id_const, 2, t1, q1, -1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e3, e8, 1);
  Problem1_.inputElement(id_const, 1, e2, e9, -1);
  Problem1_.inputElement(id_const, 2, t2, q1, -1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e2, e6, 1);
  Problem1_.inputElement(id_const, 1, e3, e5, -1);
  Problem1_.inputElement(id_const, 2, t3, q1, -1); 


  id_const++;
  Problem1_.inputElement(id_const, 1, e6, e7, 1);
  Problem1_.inputElement(id_const, 1, e4, e9, -1);
  Problem1_.inputElement(id_const, 2, t1, q2, -1); 

  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e9, 1);
  Problem1_.inputElement(id_const, 1, e3, e7, -1);
  Problem1_.inputElement(id_const, 2, t2, q2, -1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e3, e4, 1);
  Problem1_.inputElement(id_const, 1, e1, e6, -1);
  Problem1_.inputElement(id_const, 2, t3, q2, -1); 


  id_const++;
  Problem1_.inputElement(id_const, 1, e4, e8, 1);
  Problem1_.inputElement(id_const, 1, e5, e7, -1);
  Problem1_.inputElement(id_const, 2, t1, q3, -1); 

  id_const++;
  Problem1_.inputElement(id_const, 1, e2, e7, 1);
  Problem1_.inputElement(id_const, 1, e1, e8, -1);
  Problem1_.inputElement(id_const, 2, t2, q3, -1);
  
  id_const++;
  Problem1_.inputElement(id_const, 1, e1, e5, 1);
  Problem1_.inputElement(id_const, 1, e2, e4, -1);
  Problem1_.inputElement(id_const, 2, t3, q3, -1); 
 

  if (debug) std::cout << "Number of constraints: " << id_const << std::endl;
  
}



PrimalAdjugateResult PrimalAdjugate::getResult(bool debug){


  /* Solve problem now */ 
  if (debug) std::cout << "[SDP-ADJ] Getting result\n";
  Problem1_.initializeUpperTriangle();
  if (debug) std::cout << "[SDP-ADJ] Init solver\n";
  Problem1_.initializeSolve();
  
  // if necessary, dump input data and initial point
  if (debug) std::cout << "[SDP-ADJ] Solve!\n";
  Problem1_.solve();
  
  //
  if (debug) std::cout << "[SDP-ADJ] Saving results\n";
  PrimalAdjugateResult res = PrimalAdjugateResult();
  res.stop_iter = Problem1_.getIteration(); 
  res.f_opt = - Problem1_.getPrimalObj(); 
  res.d_opt = - Problem1_.getDualObj();
  res.primal_error = Problem1_.getPrimalError();
  res.dual_error = Problem1_.getDualError();
  
  if (debug) std::cout << "[SDP-ADJ] Saving lagrange multipliers\n";
  auto dual_point = Problem1_.getResultXVec(); 
  for (int i=0; i < nConstraints_; i++)
        res.dual_point(i) = dual_point[i]; 
  
  

  if (debug) std::cout << "[SDP-ADJ] Saving blocks for solutions\n";
  auto Xe = Problem1_.getResultYMat(1);  // NOTE: we need to reshape this matrix (6 x 6)
  auto Xt = Problem1_.getResultYMat(2); 

  // std::cout << "Printing Xe\n"; 
    
  for (int i=0; i<dimE_; ++i) {
    for (int j=0; j<dimE_; ++j) {
      res.primal_Xe(i, j) = Xe[i * dimE_ + j];  // 6 x 6
      }
      }
      
      std::cout << std::endl;
   
  for (int i=0; i<dimT_; ++i) {
    for (int j=0; j<dimT_; ++j) {
      res.primal_Xt(i,j) = Xt[i * dimT_ + j];  // 3 x 3
      }
      }
     
   
  return res; 
}


/** 
Get results for this problem
**/  

void PrimalAdjugate::getSolutionFromResult(PrimalAdjugateResult & res, bool debug){



    auto max_index = [] (const Eigen::Matrix<double, Eigen::Dynamic, 1> & e_vector, 
                                                        double val_rot) -> int
                {
                int i = 0;  
                while (i < e_vector.rows())
                {
                        if (e_vector[i] == val_rot) 
                                break;   
                        else i++; 
                }
                return i;
                };
                
                
                
                
                
        // Extract solution from blocks
        // This function modifies some of the fields in res
        // E_opt, t_opt, q_opt
        Matrix9 Xe = res.primal_Xe; 
        Matrix6 Xt = res.primal_Xt; 
        
       if (debug)  std::cout << "Extracting data from SDP solution\n"; 
         
         
        // Extract essential matrix
        
        Eigen::SelfAdjointEigenSolver<Matrix9> eigen_solver_M(Xe);
                
        Vector9 eigenvalues_Xe = eigen_solver_M.eigenvalues().real();
        
        double val_rot = eigenvalues_Xe.maxCoeff();  
                
                // std::cout << "Matrix M:\n" << M << std::endl;
                // std::cout << "Eigenvalues:\n" << eigenvalues_M1_rot << std::endl;
                int mm = max_index(eigenvalues_Xe, val_rot); 
                // std::cout << "Index from lambda: " << mm << std::endl;
                Vector9 ee = eigen_solver_M.eigenvectors().col(mm);
                
                if (debug)
                {
                        // std::cout << "U for Xe\n" << eigen_solver_M.eigenvectors() << std::endl;
                        std::cout << "Eigenvalues of Xe:\n" << eigenvalues_Xe.transpose(); 
                        std::cout << "\nMaximum eigenvalue = " << val_rot; 
                        std::cout << " with index " << mm << std::endl;
                        std::cout << "Eigenvector:\n" << ee.transpose() << std::endl;
                
                }
                
                // reshape ee 
                Matrix3 E_opt = Eigen::Map<Matrix3> (ee.data(), 3, 3);
            
            
                res.eigenvalues_Xe = eigenvalues_Xe; 
                res.rank_Xe = computeRankX(eigenvalues_Xe, val_rot); 
                res.stable_rank_Xe = computeStableRankX(Xe, val_rot); 
                
                
                
               // std::cout << "E from SVD before projection:\n" << E_opt << std::endl;
               const Eigen::JacobiSVD<Matrix3> svd(E_opt, Eigen::ComputeFullU | Eigen::ComputeFullV);
               // std::cout << "Singular values of E before projection:\n" << svd.singularValues() << std::endl;
            
                Matrix3 D = Matrix3::Zero(); 
                D(0,0) = 1; 
                D(1,1) = 1; 
                // project onto essential matrix space 
                res.E_opt = svd.matrixU() * D * svd.matrixV().transpose();
                
                
                
        // do the same for Xt
        
        Eigen::SelfAdjointEigenSolver<Matrix6> eigen_solver_Mt(Xt);
                
        Vector6 eigenvalues_Xt = eigen_solver_Mt.eigenvalues().real();
        
        double val_Xt = eigenvalues_Xt.maxCoeff();  
                
                int mm_tt = max_index(eigenvalues_Xt, val_Xt); 
                Vector6 t_opt = eigen_solver_Mt.eigenvectors().col(mm_tt);
                
                if (debug)
                {
                        // std::cout << "U for Xt\n" << eigen_solver_Mt.eigenvectors() << std::endl;
                        std::cout << "Eigenvalues of Xt:\n" << eigenvalues_Xt.transpose(); 
                        std::cout << "\nMaximum eigenvalue = " << val_Xt; 
                        std::cout << " with index " << mm_tt << std::endl;
                        std::cout << "Eigenvector:\n" << t_opt.transpose() << std::endl;
                
                }
                
                
                // normalize
                res.t_opt = t_opt.normalized(); 
                
                 
                res.eigenvalues_Xt = eigenvalues_Xt; 
                
                
                res.rank_Xt.push_back(computeRankX(eigenvalues_Xt, val_Xt)); 
                
                res.stable_rank_Xt.push_back(computeStableRankX(Xt, val_Xt)); 
                

        return;
}







}  // end namesapce SDPRelPose

