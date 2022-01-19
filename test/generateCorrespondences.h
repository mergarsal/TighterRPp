/*
Original from Opengv.
Adaptation for the (central case) relative pose problem
*/

#pragma once

#include <stdlib.h>

#include <functional>


#include <eigen3/Eigen/Dense>

namespace UtilsRelPose
{

    typedef Eigen::Matrix<double, 4, 1> Vector4;
    typedef Eigen::Matrix<double, 6, 1> Vector6;
    typedef Eigen::Matrix<double, 3, 1> Vector3;
    typedef Eigen::Matrix<double, 9, 1> Vector9;
    typedef Eigen::Matrix<double, 12, 1> Vector12;
    typedef Eigen::Matrix<double, 3, 3> Matrix3;
    typedef Eigen::Matrix<double, 4, 4> Matrix4;
    typedef Eigen::Matrix<double, 9, 9> Matrix9;
    typedef Eigen::Matrix<double, 3, 9> Matrix39;
    typedef Eigen::Matrix<double, 9, 3> Matrix93;
    typedef Eigen::Matrix<double, 3, 4> Matrix34;
    typedef Eigen::Matrix<double, 12, 12> Matrix12;
  
  
    /** A simple struct that contains the elements of a corresponding fatures as bearing (unit) vectors */
    struct CorrespondingFeatures {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Vector3 bearing_vector_0, bearing_vector_1;

        double weight_;
        /** Simple default constructor; does nothing */
        CorrespondingFeatures() {}

        CorrespondingFeatures(const Vector3 & bearing_vector_0,
                              const Vector3 & bearing_vector_1,
                              double weight_match = 1.0)
                                : weight_(weight_match), bearing_vector_0(bearing_vector_0),
                                bearing_vector_1(bearing_vector_1) {}

    }; // end of CorrespondingFeatures struct

    // Define types
    typedef std::vector<CorrespondingFeatures> bearing_vectors_t;
    
    typedef Eigen::VectorXd weights_t;

    /* Generate random translation */
    using GenerateTranslation = std::function<Vector3(const double max_parallax, const Vector3& direction_parallax)>;
    /* Generate random rotation */
    using GenerateRotation = std::function<Matrix3(const double max_angle, const Vector3& rot_angle)>;

                  
    /** A simple struct that contains the elements of a corresponding fatures as bearing (unit) vectors */
    struct StrInputGenProblem {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        double FoV = 100;                    // in degrees
        double min_depth = 1.0;              // in meters
        double max_depth = 8.0;              // in meters 
        bool allow_coplanar = false;         // allow to generate synthetic scenes with coplanar points
        double max_X = 20.0;                 // in meters: this value is the absolute value. Max value for X-axis allowed 
        double focal_length = 800;           // in pixels
        size_t n_points = 100;
        double noise = 0.1;                  // in pixels
        double outlier_fraction = 0.0;
        double min_epipolar_error_sq = 0.0;  // min epipolar error for the outliers

        // params for relpose
        double parallax = 2.0;               // in meters
        double max_rotation = 0.5;           // in degrees

        /* Not used */
        bool use_custom_t = false;           // true if you want to provide t
        bool use_custom_R = false;           // true if you want to provide R
        Vector3 custom_t; 
        Matrix3 custom_R;  
        /* Until here */

        Vector3 dir_trans; 
        Vector3 dir_rotation;

        // constructor
        StrInputGenProblem(){}; 

    }; // end of StrInputGenProblem struct
                  

    struct StrOutputGenProblem {
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        Vector3  translation;
        Matrix3  rotation;

        bearing_vectors_t  points_correspondences;
        Eigen::MatrixXd  points_3D;
        std::vector<int>  indices_outliers;

        // constructor
        StrOutputGenProblem(size_t n_points) {
        translation.setZero(); 
        rotation.setZero(); 
        points_correspondences; 
        points_3D = Eigen::MatrixXd(3, n_points); 
        indices_outliers; 
        }; 

    }; // end of StrOutputGenProblem struct
                
    // Generate 3D point inside a frustum
    Vector3 generateRandomPointTruncated(double FoV, double min_depth, double max_depth, bool allow_coplanar, double max_X );
    // add noise in plane            
    Vector3 addNoiseGaussian(Vector3 clean_point, double noise_level);
    // generate random translation with the specified norm
    Vector3 generateRandomTranslationDefault( double max_parallax, const Vector3 & dir_parallax);
    // generate random rotation with maxAngle
    Matrix3 generateRandomRotationDefault( double maxAngle, const Vector3 & dir_rot);

    Matrix9 constructDataMatrix(const bearing_vectors_t & bearing_vectors); 
    
    Matrix3 computeEfromRt(const Matrix3 & R, const Vector3 & t); 
    
    void computeRtfromE(const bearing_vectors_t & points,
                        const Matrix3& E, Matrix3 & R,
                        Vector3 & t ); 
                        
                        
    StrOutputGenProblem createSyntheticExperiment(const StrInputGenProblem & in_param, 
                                  const GenerateTranslation& generateTranslation = generateRandomTranslationDefault,
                                  const GenerateRotation& generateRotation = generateRandomRotationDefault); 

}  // end of namespace UtilsRelPose
