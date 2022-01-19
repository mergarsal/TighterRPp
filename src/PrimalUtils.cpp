# include "PrimalUtils.h"


namespace RelPose
{
        template <typename Vector, typename MatrixT, typename VectorT>
        void PrimalClassBase<Vector, MatrixT, VectorT>::printResult(const PrimalResult<Vector, MatrixT, VectorT> & res)
        {
        
         
  /* Show results */
	
  fprintf(stdout, "\nStop iteration = %d\n",
	  res.stop_iter);
  fprintf(stdout, "objValPrimal   = %+10.6e\n",
	  res.f_opt);
  fprintf(stdout, "objValDual     = %+10.6e\n",
	  res.d_opt);
  fprintf(stdout, "p. feas. error = %+10.6e\n",
	  res.primal_error);
  fprintf(stdout, "d. feas. error = %+10.6e\n\n",
	  res.dual_error);

  
  

  // Problem1_.printComputationTime(stdout);


  std::cout << "Spectral properties of the solution\n"; 
  
  std::cout << "Rank Xe: " << res.rank_Xe; 
  std::cout << "\nStable rank Xe: " << res.stable_rank_Xe;
  std::cout << "\nRank Xt: " << res.rank_Xt; 
  std::cout << "\nStable rank Xt: " << res.stable_rank_Xt;
  
        
                
        }
        
        
        // Estimate rank matrix 
        int computeRankX(const Eigen::Matrix<double, Eigen::Dynamic, 1> & e_vector,double max_eig)
        {
        
        // consider for rank all the eigenvalues with 
        // value greater than max_eig
        int rank = 0; 
        
        for (int i = 0; i < e_vector.rows(); i++) 
                if (e_vector(i) >= 0.001 * max_eig) rank++;       
        
        return rank; 
        }
        
        
        // Estimate stable rank matrix 
        double computeStableRankX(const Eigen::MatrixXd & X_matrix, 
                                                                double max_eig)
        {
        
        // consider for rank all the eigenvalues with 
        // value greater than max_eig
        double norm_frob = 0; 
        
        for (int i = 0; i < X_matrix.rows(); i++) 
        {
                for (int j = 0; j < X_matrix.cols(); j++)
                        norm_frob += X_matrix(i, j) * X_matrix(i, j);        
        }
                
                     
        
        return norm_frob / ( max_eig * max_eig); 
        }
        
        


}  // end of namespace PrimalRotPrior
