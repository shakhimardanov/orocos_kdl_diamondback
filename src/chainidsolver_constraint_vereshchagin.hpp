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
#include "articulatedbodyinertia.hpp"

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
        ChainIdSolver_Constraint_Vereshchagin(const Chain& chain,Vector root_acc,unsigned int nr_of_constraints, const Jacobian& alfa, const JntArray& beta);
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
        
    private:
        //Functions to calculate velocity, propagated inertia, propagated bias forces, constraint forces and accelerations
        void initial_upwards_sweep(const JntArray &q, const JntArray &q_dot, const Wrenches& f_ext);
        void downwards_sweep();
        void constraint_calculation();
        void final_upwards_sweep(JntArray &q_dotdot, JntArray &torques);

        Chain chain;
        unsigned int nj;
        unsigned int ns;
        unsigned int nc;
        Twist acc_root;

        typedef Eigen::Matrix<double,6,1> Vector6d;
        typedef Eigen::Matrix<double,6,6> Matrix6d;
        typedef Eigen::Matrix<double,6,Eigen::Dynamic> Matrix6Xd;
        
        struct segment_info{
            Frame F;//Pose
            Frame F_matrix;//Transformation matrix
            Twist Z;//Unit twist
            Twist v;//twist
            Twist acc;//acceleration twist
            Wrench f;//wrench
            Wrench U;//wrench p of the bias forces in matrix form
            Wrench R;//wrench p of the bias forces in matrix form
            Wrench R_tilde;//vector of wrench p of the bias forces (new) in matrix form
            Twist C;//constraint
            ArticulatedBodyInertia H;//I (expressed in 6*6 matrix) 
            ArticulatedBodyInertia P;//I (expressed in 6*6 matrix) 
            ArticulatedBodyInertia P_tilde;//I (expressed in 6*6 matrix) 
            Wrench PZ;//vector U[i] = I_A[i]*S[i]
            Wrench PC;//vector E[i] = I_A[i]*c[i]
            double D;//vector D[i] = S[i]^T*U[i]
            Matrix6Xd E;//matrix with virtual unit constraint force due to acceleration constraints
            Matrix6Xd E_tilde;
            Eigen::MatrixXd M;//acceleration energy already generated at link i
            Eigen::VectorXd G;//magnitude of the constraint forces already generated at link i
            Matrix6d UDS;//UDS =I_A*S*D_inv*S[i]^T;
            Vector6d K;//K[i] = E_constr_A[i]^T*S[i]
            double u;//vector u[i] = torques(i) - S[i]^T*(p_A[i] + I_A[i]*C[i])
            segment_info(unsigned int nc){
                E.resize(6,nc);
                E_tilde.resize(6,nc);
                G.resize(nc,nc);
                M.resize(nc,nc);
            };
        };

        std::vector<segment_info> results;
        
        Jacobian alfa_N;
        Eigen::VectorXd v_constr;//v = M_0^-1*(beta_N - E_constr_a^T*a[0] - G_constr[0])
        Eigen::MatrixXd M_0_inverse;
        JntArray beta_N;
        Eigen::VectorXd nu,nu_sum;
        Wrench qdotdot_sum;
    };
}

#endif
