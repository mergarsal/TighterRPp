#pragma once 


#include <Eigen/Core>




namespace RelPose{

        typedef Eigen::Matrix<double, 3, 1> Vector3;
        
        typedef Eigen::Matrix<double, 4, 1> Vector4;

        typedef Eigen::Matrix<double, 5, 1> Vector5;
        
        typedef Eigen::Matrix<double, 6, 1> Vector6;
        
        typedef Eigen::Matrix<double, 7, 1> Vector7;
        
        typedef Eigen::Matrix<double, 8, 1> Vector8;
        
        typedef Eigen::Matrix<double, 9, 1> Vector9;
        
        typedef Eigen::Matrix<double, 10, 1> Vector10;
        
        typedef Eigen::Matrix<double, 12, 1> Vector12;
        
        typedef Eigen::Matrix<double, 13, 1> Vector13;
        
        typedef Eigen::Matrix<double, 22, 1> Vector22;
        
        
        typedef Eigen::Matrix<double, 2, 2> Matrix2;
        
        typedef Eigen::Matrix<double, 3, 3> Matrix3;
        
        typedef Eigen::Matrix<double, 6, 6> Matrix6;
        
        typedef Eigen::Matrix<double, 9, 9> Matrix9;
        
        typedef Eigen::Matrix<double, 12, 12> Matrix12;

        typedef Eigen::Matrix<double, 3, 4> Matrix34;



        /** A simple struct that contains the elements of a corresponding fatures as bearing (unit) vectors */
       /*
        struct CorrespondingFeatures {
                EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
          Vector3 bearing_vector_0, bearing_vector_1;

          double weight_;
          CorrespondingFeatures() {}

          CorrespondingFeatures(const Vector3 & bearing_vector_0,
                                const Vector3 & bearing_vector_1,
                                double weight_match = 1.0)
                                : weight_(weight_match), bearing_vector_0(bearing_vector_0),
                                bearing_vector_1(bearing_vector_1) {}

        }; // end of CorrespondingFeatures struct

       // Define types
       typedef std::vector<CorrespondingFeatures> bearing_vectors_t;
       */
       
       
} // end of essential namespace
