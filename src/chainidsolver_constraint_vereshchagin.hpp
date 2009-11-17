// Copyright  (C)  2009 

// Version: 1.0
// Author: 
// Maintainer: 
// URL: http://www.orocos.org/kdl

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef KDL_CHAIN_IDSOLVER_VERESHCHAGIN_HPP
#define KDL_CHAIN_IDSOLVER_VERESHCHAGIN_HPP

#include "chainidsolver.hpp"

namespace KDL{
    /**
     * \brief Dynamics calculations by constraints based on Vereshchagin 1989.
     * Code also based on Roy Featherstone, Rigid body dynamics Algorithms page 132 and 127.
     * 
     * For a chain
     */
    class ChainIdSolver_Constraint_Vereshchagin : public ChainIdSolver{
    public:
        /**
         * Constructor for the solver, it will allocate all the necessary memory
         * \param chain The kinematic chain to calculate the inverse dynamics for, an internal copy will be made.
         * \param root_acc The acceleration vector of the root to use during the calculation.(most likely contains gravity)
         *
         */
        ChainIdSolver_Constraint_Vereshchagin(const Chain& chain,Vector root_acc);
        ~ChainIdSolver_Constraint_Vereshchagin(){};
        
        /**
         * Function to calculate from Cartesian forces to joint torques.
         * Input parameters;
         * \param q The current joint positions
         * \param q_dot The current joint velocities
         * \param f_ext The external forces (no gravity) on the segments
         * Output parameters:
         * \param q_dotdot The joint accelerations
         * \param torques the resulting torques for the joints
         */
        int CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques);
        
        //Function to calculate the determinant
        double determ(double[6][6],double);
        
        //Functions to calculate velocity, propagated inertia, propagated bias forces, constraint forces and accelerations
        void calc_one(int i, int j, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext);
        void calc_two(int i, int k_n);
        void constraint_calc(int k_n);
        void calc_three(int i, int j, int k_n, JntArray &q_dotdot, JntArray &torques);

    private:
        Chain chain;
        unsigned int nj;
        unsigned int ns;

        Twist acc_root;
        
        struct segment_info{
            Frame X;//Pose
            Twist S;//Unit twist
            Twist v;//twist
            Twist a;//acceleration twist
            Twist c;//constraint
            Wrench f;//wrench
            Wrench p_A_wrench;//dunno
            Eigen::Matrix<6,1,double> p_A;//wrench p of the bias forces in matrix form
            Eigen::Matrix<6,6,double> I_A;//I (expressed in 6*6 matrix) 
            Eigen::Matrix<6,1,double> U;//vector U[i] = I_A[i]*S[i]
            Eigen::Matrix<6,1,double> E;//vector E[i] = I_A[i]*c[i]
            double D;//vector D[i] = S[i]^T*U[i]
            Eigen::Matrix<6,1,double> p_a;//vector of wrench p of the bias forces (new) in matrix form
            Eigen::Matrix<6,6,double> I_a;//vector of RigidBodyInertia (new)
            Eigen::Matrix<6,6,double> X_matrix;//Transformation matrix
            Eigen::Matrix<6,6,double> X_matrix_inv;//Inverse Transformation matrix
            Eigen::Matrix<6,6,double> G;//G[i] = I_a*X_matrix_inv[i]
            Eigen::Matrix<6,1,double> F;//F[i]
            Eigen::Matrix<6,6,double> E_constr_A;//matrix with virtual unit constraint force due to acceleration constraints
            Eigen::Matrix<6,6,double> E_constr_a;
            Eigen::Matrix<6,6,double> M_constr;//acceleration energy already generated at link i
            Eigen::Matrix<6,1,double> G_constr;//magnitude of the constraint forces already generated at link i
            Eigen::Matrix<1,6,double> U_E;//U_E[i] = S[i]^T*E_constr_A[i];
            Eigen::Matrix<6,1,double> K;//K[i] = E_constr_A[i]^T*S[i]
            Eigen::Matrix<6,1,double> acc;
            double u;//vector u[i] = torques(i) - S[i]^T*(p_A[i] + I_A[i]*C[i])
            
        };
        double alfa_N[6][6];
        double v_constr[6][1];//v = M_0^-1*(beta_N - E_constr_a^T*a[0] - G_constr[0])
        double M_0_inverse[6][6];
        double beta_N[6][1];
        double v_constr_sum[6][1];
        double c_sum_a[6][1];
        double qdotdot_sum[6][1];
    };
}

#endif
