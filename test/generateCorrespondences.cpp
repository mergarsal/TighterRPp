#include "generateCorrespondences.h"

#include <vector>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <memory>

#include <Eigen/Core>

#include <Eigen/SVD> // required for SVD decomposition in unlift method

#include <functional>


#define PI 3.14159
#define SQRT2 1.41421356237



namespace UtilsRelPose
{


    void computeRtfromE(const bearing_vectors_t & points,
                        const Matrix3& E, Matrix3 & R,
                        Vector3 & t )
    {

        Eigen::JacobiSVD<Matrix3> svd(E, Eigen::ComputeFullU | Eigen::ComputeFullV);

        Matrix3 W;
        W.setZero();
        W(0, 1) = -1; W(1, 0) = 1; W(2, 2) = 1;

        Matrix3 R1, R2, V, U;
        U = svd.matrixU();
        V = svd.matrixV();

        // Rotation matrix
        R1 = U * W * V.transpose();

        R2 = U * W.transpose() * V.transpose();


        // for R1, R2 in SO(3)
        if (R1.determinant() < 0)       R1 = -R1;
        if (R2.determinant() < 0)       R2 = -R2;

        Vector3 t1;
        t1 = U.col(2);

        // Cheirality

        unsigned int N = points.size();
        Eigen::MatrixXd X(9, N);
        Eigen::MatrixXd Y(6, N);
        X.setZero();
        Y.setZero();


        for (int i = 0; i < N; i++)
        {
            Vector9 temp;
            Vector3 v1 = points[i].bearing_vector_1;
            Vector3 v0 = points[i].bearing_vector_0;
            double weight_i = points[i].weight_;

            temp.setZero();
            for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 0) = weight_i * v1[j] * v0;

            X.col(i) = temp;

            Y.block<3, 1>(0, i) = v0 * v1.norm(); Y.block<3, 1>(3, i) = v1 * v0.norm();
        }

        Matrix3 ER1 = E * E.transpose() * R1;
        Matrix3 ER2 = E * E.transpose() * R2;

        Eigen::MatrixXd M11(N, 1);
        Eigen::MatrixXd M12(N, 1);
        M11 = X.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1> > (ER1.data(), 9, 1);
        M12 = X.transpose() * Eigen::Map<Eigen::Matrix<double, 9, 1> > (ER2.data(), 9, 1);

        double n_non_zeros_M11 = (M11.array() > 0).count();
        double n_non_zeros_M12 = (M12.array() > 0).count();

        if (n_non_zeros_M11 >= n_non_zeros_M12)      R = R1;
        else                                         R = R2;

        // Compute the translation

        Eigen::MatrixXd M21(N, 1), M22(N, 1);
        Eigen::MatrixXd R_eye(6, 3); R_eye.block<3,3>(0,0) = - R.transpose(); R_eye.block<3,3>(3, 0) = Matrix3::Identity();


        M21 = Y.transpose() * R_eye * t1;

        double n_non_zeros_M21 = (M21.array() > 0).count();
        double n_non_zeros_M22 = (-M21.array() > 0).count();



        if (n_non_zeros_M21 >= n_non_zeros_M22)    t = t1;
        else                                       t = -t1;

    }





    Matrix9 constructDataMatrix(const bearing_vectors_t & bearing_vectors)
    {
        Matrix9 C;
        Matrix9 temp;
        // clean output matrix
        C.setZero();

        for (int i = 0; i < bearing_vectors.size(); i++)
        {
            // clean data
            temp.setZero();
            // C_i.setZero();
            const Vector3 v1 = bearing_vectors[i].bearing_vector_1;
            const Vector3 v0 = bearing_vectors[i].bearing_vector_0;
            const double weight = bearing_vectors[i].weight_;
            for (int j = 0; j < 3; j++)    temp.block<3, 1>(j*3, 1) = v1[j] * v0;
            C += weight * temp * temp.transpose();
        }
        return 0.5 * (C + C.transpose());
    }
    
    
    Matrix3 computeEfromRt(const Matrix3 & R, const Vector3 & t)
        {
            Matrix3 t_skew = Matrix3::Zero();          
            t_skew(0, 1) = -t(2);   t_skew(0, 2) = t(1);
            t_skew(1, 0) = t(2);    t_skew(1, 2) = -t(0);
            t_skew(2, 0) = -t(1);   t_skew(2, 1) = t(0);
         
                   
            Matrix3 E = t_skew * R;
            return (E / (E.norm()) * SQRT2);
        }
    
    
    Vector3 generateRandomTranslationDefault(const double parallax, 
                                                      const Vector3 & dir_parallax)
    {               
        Vector3 translation;
        translation[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
        translation[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0;
        translation[2] = -(((double) rand())/ ((double) RAND_MAX));

        // std::cout << "Translation created!\n";
        return (parallax * translation.normalized());
    }
        

        
    // From OpenGV
    Matrix3 generateRandomRotationDefault(const double maxAngle, 
                                  const Vector3 & dir_rot)
    {

        Vector3 rpy;
        rpy[0] = ((double) std::rand())/ ((double) RAND_MAX);
        rpy[1] = ((double) std::rand())/ ((double) RAND_MAX);
        rpy[2] = ((double) std::rand())/ ((double) RAND_MAX);

        rpy[0] = maxAngle*2.0*(rpy[0]-0.5);
        rpy[1] = maxAngle*2.0*(rpy[1]-0.5);
        rpy[2] = maxAngle*2.0*(rpy[2]-0.5);

        Matrix3 R1;
        R1(0,0) = 1.0;
        R1(0,1) = 0.0;
        R1(0,2) = 0.0;
        R1(1,0) = 0.0;
        R1(1,1) = cos(rpy[0]);
        R1(1,2) = -sin(rpy[0]);
        R1(2,0) = 0.0;
        R1(2,1) = -R1(1,2);
        R1(2,2) = R1(1,1);

        Matrix3 R2;
        R2(0,0) = cos(rpy[1]);
        R2(0,1) = 0.0;
        R2(0,2) = sin(rpy[1]);
        R2(1,0) = 0.0;
        R2(1,1) = 1.0;
        R2(1,2) = 0.0;
        R2(2,0) = -R2(0,2);
        R2(2,1) = 0.0;
        R2(2,2) = R2(0,0);

        Matrix3 R3;
        R3(0,0) = cos(rpy[2]);
        R3(0,1) = -sin(rpy[2]);
        R3(0,2) = 0.0;
        R3(1,0) =-R3(0,1);
        R3(1,1) = R3(0,0);
        R3(1,2) = 0.0;
        R3(2,0) = 0.0;
        R3(2,1) = 0.0;
        R3(2,2) = 1.0;

        Matrix3 rotation = R3 * R2 * R1;
        rotation.col(0) = rotation.col(0) / rotation.col(0).norm();
        rotation.col(2) = rotation.col(0).cross(rotation.col(1));
        rotation.col(2) = rotation.col(2) / rotation.col(2).norm();
        rotation.col(1) = rotation.col(2).cross(rotation.col(0));
        rotation.col(1) = rotation.col(1) / rotation.col(1).norm();

        return rotation;
    }
        

    Vector3 generateRandomPointTruncated(double FoV, double min_depth, 
                                         double max_depth, bool allow_coplanar, double max_X)
    {

        Vector3 point3D;
        // Note: tan(X) applies mod(X, pi) before computing the actual tan
        // depth
        point3D[2] = min_depth + (max_depth - min_depth) * ( ((double) std::rand() / (double) RAND_MAX));

        double xmax = (tan(FoV * 0.5 * PI / 180)) * point3D[2];
        if (xmax <= 0) xmax *= -1;  


        if (!allow_coplanar)
        // update xmax
        xmax = std::min(xmax, max_X);

        point3D[0] = xmax * (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;
        point3D[1] = xmax * (((double) rand()) / ((double) RAND_MAX)-0.5)*2.0;

        return point3D;
    }


    Vector3 addNoiseGaussian(Vector3 clean_point, double noise_level)
    {
        Vector3 noisy_point = clean_point; 
        double last_entry = noisy_point(2); 
        noisy_point(0) /= last_entry; 
        noisy_point(1) /= last_entry; 
        noisy_point(2) = 1;


        double x_noise = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 * noise_level;
        double y_noise = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 * noise_level;

        noisy_point(0) += x_noise; 
        noisy_point(1) += y_noise;

        return noisy_point; 
    }


    StrOutputGenProblem createSyntheticExperiment(const StrInputGenProblem & in_param, 
                                                  const GenerateTranslation& generateTranslation, 
                                                  const GenerateRotation& generateRotation)
    {
      // struct for the ouput 
      StrOutputGenProblem str_out(in_param.n_points);
      
      // generate a random pointcloud
      str_out.points_3D.setZero();

      
      for( size_t i = 0; i < in_param.n_points; i++ )
            str_out.points_3D.col(i) = generateRandomPointTruncated(in_param.FoV, in_param.min_depth,
                                    in_param.max_depth, in_param.allow_coplanar, in_param.max_X );
                                                    

      // 2. Create the relative rotation and translation
      // 2.1 Fix first frame
      Matrix3 R1 = Matrix3::Identity();
      Vector3 T1 = Vector3::Zero();

      // 2.2 Fix second frame
      Matrix3 rotation = Matrix3::Zero();
      Vector3 T2 = Vector3::Zero();

      bool valid_position = false;
      Eigen::MatrixXd points3D_frame2(3, in_param.n_points); 
      points3D_frame2.setZero();
      
      int max_n_iters = 200, n_iters = 0;

      
      
    do
    {
        // create random translation
        // note: we always create a backward movement

        T2 = generateTranslation(in_param.parallax, in_param.dir_trans);


        // do the sae with rotation 
        rotation = generateRotation(in_param.max_rotation, in_param.dir_rotation);


        // compute relative coordinates
        points3D_frame2.setZero();
        for( size_t i = 0; i < in_param.n_points; i++ )
        {

          Vector3 temp = rotation.transpose() * (str_out.points_3D.col(i)- T2);
          points3D_frame2.col(i) = Eigen::Map<Vector3>(temp.data(), 3, 1);
         
        }

    
        // check condition for FoV
        Eigen::VectorXd ratio_X = (points3D_frame2.row(0)).cwiseQuotient((points3D_frame2.row(2)));
        
        double max_ratio_x = (ratio_X.array().abs()).maxCoeff();


        // check if any point was outside the FoV

        if (max_ratio_x > tan(in_param.FoV * 0.5 * PI / 180))
        {
            std::cout << "At least one point did not fulfill the condition\n"; 
            n_iters++;
            continue;
        }

        Eigen::VectorXd ratio_Y = (points3D_frame2.row(1)).cwiseQuotient((points3D_frame2.row(2)));
        
        double max_ratio_y = (ratio_Y.array().abs()).maxCoeff();
        // check if any point was outside the FoV

        if (max_ratio_y > tan(in_param.FoV * 0.5 * PI / 180))
        {
            std::cout << "At least one point did not fulfill the condition\n"; 
            n_iters++;
            continue;
        }

        // if we have arrive here, the position is valid
        valid_position = true;

    }while(!valid_position && (n_iters <= max_n_iters));

    // check invalid rotation 
    if ((!valid_position) && (n_iters > max_n_iters))
        std::cout << "[ERROR] Pose is not valid.\nSome points do not lie in the FOV\n";


    // save pose
    str_out.translation = T2.normalized(); 
    str_out.rotation = rotation;


    // Generate the correspondences
    Matrix3 K; 
    K.setZero();
    K(0, 0) = in_param.focal_length; 
    K(1, 1) = in_param.focal_length; 
    K(0, 2) = tan(in_param.FoV * 0.5 * PI / 180) * in_param.focal_length; 
    K(1, 2) = tan(in_param.FoV * 0.5 * PI / 180) * in_param.focal_length; 
    K(2, 2) = 1; 
  
    Matrix3 invK = K.inverse(); 


    // std::cout << " Creating the correspondences\n";
    for( size_t i = 0; i < in_param.n_points; i++ )
        {
        Vector3 obs1, obs2, o1, o2;
        obs1 = str_out.points_3D.col(i);
        obs2 = points3D_frame2.col(i);
            
        double last_o1 = obs1(2), last_o2 = obs2(2); 


        // observations in hom. coord.
        o1 << obs1(0) / last_o1, obs1(1) / last_o1, 1;  
        o2 << obs2(0) / last_o2, obs2(1) / last_o2, 1;  
          
        //add noise
        obs1 = K * o1;
        obs2 = K * o2;

        if(in_param.noise > 0.0 )
        {
          obs1 = addNoiseGaussian(obs1, in_param.noise);
          obs2 = addNoiseGaussian(obs2, in_param.noise);
        }

        CorrespondingFeatures correspondence;
        correspondence.bearing_vector_0 = (invK * obs1).normalized();
        correspondence.bearing_vector_1 = (invK * obs2).normalized(); 
        correspondence.weight_ = 1.0;

        // add it to the std::vector
        str_out.points_correspondences.push_back(correspondence);
        }

    /*   OUTLIERS   */

    // add outliers
    size_t number_outliers = (size_t) floor(in_param.outlier_fraction * in_param.n_points);
    size_t max_number_iters = 50;
    for(size_t i = 0; i < number_outliers; i++)
    {
        size_t i_iter = 0;
        bool valid_outlier = false;
        do
            {
                
            Vector3 outlier_correspondence = generateRandomPointTruncated(170, in_param.min_depth,
                                        in_param.max_depth, true, in_param.max_X);


            //normalize the bearing vector
            outlier_correspondence = rotation.transpose() * (outlier_correspondence - T2);
            outlier_correspondence.normalize();


            str_out.points_correspondences[i].bearing_vector_1 = outlier_correspondence;
            str_out.indices_outliers.push_back(i);

            valid_outlier = true;

            if ((!valid_outlier) && (i_iter > max_number_iters))
                {
                // break the loop
                outlier_correspondence[0] = 1.0;
                outlier_correspondence[1] = 1.0;
                outlier_correspondence[2] = 1.0;

                outlier_correspondence.normalize();

                str_out.points_correspondences[i].bearing_vector_1 = outlier_correspondence;
                 str_out.indices_outliers.push_back(i);

                valid_outlier = true;
                }

            // increase the counter
            i_iter += 1;


            }while(!valid_outlier);

    }  // end of: for each outlier

  return str_out; 

}


}  // end of namespace UtilsRelPose
