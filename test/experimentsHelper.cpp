#include "experimentsHelper.h"



template <typename Scalar=double>
void readSetParams(ifstream & in_file, int & N, vector<Scalar> & arr)
{
        
        string line_i, line2; 
        
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        // get nine components
        std::stringstream ss(line_i);
        
        getline(ss, line2, ',');
        N = (int)atof(line2.c_str());

        // clear arr if it is not empty 
        if (!arr.empty()) arr.clear();
        
        for (size_t j_col = 0; j_col < N; j_col++)
            {
              getline(ss, line2, ',');
              arr.push_back((double)atof(line2.c_str()));
            }
            
        // show
        std::cout << "Set of params read:\n";
        for (size_t j_col = 0; j_col < N; j_col++)
            {
              std::cout << arr[j_col] << " "; 
            }
        std::cout << std::endl;
        
}



bool readOptionsFromFile(string name_file, 
                         SceneOptions & options)
{
        /* SceneOptions

        double max_iter;  // # Number of iterations 
        
        vector<double> FoV;           // # FoV

        vector<double> max_parallax;  // # max parallax
        
        double min_depth; 
        
        double max_depth; 
        
        vector<double> outlier_fraction;
        
        vector<double> focal_length; 
        
        vector<double> noise; 
        
        vector<int> number_correspondences; 
        */
        

        // aux. strings
        string line_i, line2;


        // try to open the file 
        // open file
        ifstream in_file(name_file);

        // check if file exists
        if (!in_file.good()) return false;
        
        
        /* READ PARAMETERS !!! */ 
        
        /* MAXIMUM NUMBER OF ITERATIONS */ 
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.max_iter = (double)atof(line_i.c_str());
        std::cout << "Maximum number of iterations: " << options.max_iter << std::endl;
        

        /* FOV */
        /*
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        // get nine components
        std::stringstream ss(line_i);
        
        getline(ss, line2, ',');
        int n_entries =  = (int)atof(line2.c_str());

        for (size_t j_col = 0; j_col < n_entries; j_col++)
            {
              getline(ss, line2, ',');
              FoV.push_back((double)atof(line2.c_str()));
            }
            
        // show
        std::cout << "FoV:\n";
        for (size_t j_col = 0; j_col < n_entries; j_col++)
            {
              std::cout << FoV[j_col] << " "; 
            }
        std::cout << std::endl; 
        */
        std::cout << "Field of view:\n";
        readSetParams<double>(in_file, options.n_fov, options.FoV);
        
        
        /* MAXIMUM PARALLAX */
        std::cout << "Maximum parallax:\n";
        readSetParams<double>(in_file, options.n_parallax, options.max_parallax);
            
        

        /* MINIMUM DEPTH POINTS */
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.min_depth = (double)atof(line_i.c_str());
        std::cout << "Minimum depth points: " << options.min_depth << std::endl;
        
        /* MAXIMUM DEPTH POINTS */
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.max_depth = (double)atof(line_i.c_str());
        std::cout << "Maximum depth points: " << options.max_depth << std::endl;
        
        
        /* Outlier fraction*/
        std::cout << "Outlier fraction:\n";
        readSetParams<double>(in_file, options.n_outliers, options.outlier_fraction);
        
        /* FOCAL LENGTH */
        std::cout << "Focal length:\n";
        readSetParams<double>(in_file, options.n_focal, options.focal_length);
        
        /* NOISE */
        std::cout << "Noise:\n";
        readSetParams<double>(in_file, options.n_noise, options.noise);
             
        
        /* NUMBER OF CORRESPONDENCES */
        std::cout << "Number of correspondences:\n";
        readSetParams<int>(in_file, options.n_corr, options.number_correspondences);
            
    
        /* TRANSLATION */ 
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.use_custom_t = (bool)atof(line_i.c_str());
        std::cout << "Use custom t:\n" << options.use_custom_t << std::endl;
        
        /* Read translation vector */ 
        std::cout << "Translation vector:\n";
        int dum_var;
        std::vector<double> temp_t; 
        readSetParams<double>(in_file, dum_var, temp_t);
        
        for (int i = 0; i < dum_var; i++)
                options.custom_t(i) = temp_t[i];
        
        
        /* ROTATION */ 
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.max_rotation = (double)atof(line_i.c_str());
        std::cout << "Max rotation:\n" << options.max_rotation << std::endl;
        
        /* METHOD INIT */ 
        // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.method_init = (double)atof(line_i.c_str());
        std::cout << "Method init:\n" << options.method_init << std::endl;
        
        
        /* USE MINIMAL FOR SOLVERS */ 
         // discard first line (comment)
        getline(in_file, line_i);
        // get param
        getline(in_file, line_i);
        options.use_all_corr = (bool)atof(line_i.c_str());
        std::cout << "Use all points for solvers?:\n" << options.use_all_corr << std::endl;
        
        
         // Close the file
        in_file.close();
        
        
        return true;
}


