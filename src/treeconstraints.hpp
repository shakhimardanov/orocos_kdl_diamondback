// Copyright  (C)  2009

// Version: 1.0
// Author: Carl Wouters
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

#ifndef KDL_CHAIN_IKSOLVER_RECURSIVE_NEWTON_EULER_HPP
#define KDL_CHAIN_IKSOLVER_RECURSIVE_NEWTON_EULER_HPP

#include "treeidsolver.hpp"

namespace KDL{
    /**
     * \brief 
     * 
     * Dynamics calculations by constraints based on Vereshchagin 1989 (Modelling and Control of Motion of Manipulational Robots).
     * Code also based on Roy Featherstone, Rigid body dynamics Algorithms page 132 and 127.
     * 
     * Expanded to Tree structures.
     */
    class Tree_c : public TreeIdSolver{
    public:
        /**
         * Constructor for the solver, it will allocate all the necessary memory
         * \param chain The kinematic chain to calculate the inverse dynamics for, an internal copy will be made.
         * \param root_acc The acceleration vector of the root to use during the calculation.(most likely contains gravity)
         *
         */
        Tree_c(const Tree& tree,Vector root_acc);
        ~Tree_c(){};
        
        /**
         * Function to calculate from Cartesian forces to joint torques.
         * Input parameters;
         * \param q The current joint positions
         * \param q_dot The current joint velocities
         * \param f_ext The external forces (no gravity) on the segments
         * Output parameters:
         * \param q_dotdot The joint accelerations
         * \param torques the resulting torques for the joints, this can be an Input parameter for forward dynamics
         */
        int CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques);
        
        //Function to calculate the determinant of a matrix
        //For calculation of inverse matrix (Cramer`s rule), will be changed to an svd calculation
        double determ(double num[6][6],double k);

		//Functions to calculate velocity, propagated inertia, propagated bias forces, constraint forces and accelerations
        void calc_one(const SegmentMap::const_iterator& it, int j, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext);
        void calc_two(const SegmentMap::const_iterator& it, int k_n[10]);
        void constraint_calc(int k_total, int k_n[10]);
        void calc_three(const SegmentMap::const_iterator& it, int j, int k_n[10], int k_total, JntArray &q_dotdot, JntArray &torques);
        
        //Functions for the recursive calculation for a tree structure (depth-first and Up recursion)
        void recursiveDepthFirst(const SegmentMap::const_iterator& it, int j, int k_n[10], int k_total, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques);//calculation of propagation inertia`s and constraints
        void recursiveUP(const SegmentMap::const_iterator& it, int j, int k_n[10], int k_total, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques);//calculation of accelerations


    private:
        Tree tree;
        unsigned int nj;
        unsigned int ns;
        std::vector<Frame> X;
        std::vector<Twist> S;
        std::vector<Twist> v;
        std::vector<Twist> a;
		std::vector<Twist> c;
        std::vector<Wrench> f;
		std::vector<Wrench> p_A_wrench;
        Twist acc_root;
        
        /*
         * Changes to TreeSegment:
         * 
         * 
         * variables in treesegment: p_A (articulated bias forces), I_A (articulated inertia), X_matrix (transformation matrix, for transf),
         * 							x_matrix_inv (transformation matrix, after transf), E_constr_A (Vereshchagin 1989, point c), M_constr (Vereshchagin 1989, point d),
         * 							G_constr (Vereshchagin 1989, point e), D
         * 
         * local variables: p_a, U, E, I_a, G, E_constr_a, U_E, K
         * 
         * 
         * Constraints class:
         * 
         * alfa_N, beta_N, k_n
         * 
         * 
         */   
        
                
        
        //Better solution: make vectors from matrices
        
        double p_A[6][1][10];//vector of wrench p of the bias forces in matrix form
        double I_A[6][6][10];//vector of I (expressed in 6*6 matrix) 
        double U[6][1][10];//vector U[i] = I_A[i]*S[i]
        double E[6][1][10];//vector E[i] = I_A[i]*c[i]
        double D[10];//vector D[i] = S[i]^T*U[i]
        double p_a[6][1];//vector of wrench p of the bias forces (new) in matrix form
        double I_a[6][6];//vector of RigidBodyInertia (new)
        double X_matrix[6][6][10];//Transformation matrix
        double X_matrix_inv[6][6][10];//Inverse Transformation matrix
        double G[6][6][10];//G[i] = I_a*X_matrix_inv[i]
        double F[6][1][10];//F[i]
		double E_constr_A[6][6][10][10];//matrix with virtual unit constraint force due to acceleration constraints (4th dimension is for different constraints on other segments)
		double E_constr_a[6][6][10][10];//(4th dimension is for different constraints on other segments)
		double M_constr[6][6][10][10];//acceleration energy already generated at link i (4th dimension is for different constraints on other segments)
		double G_constr[6][1][10][10];//magnitude of the constraint forces already generated at link i (4th dimension is for different constraints on other segments)
		double alfa_N[6][6][10];//(extra dimension is for different constraints on other segments)
		double U_E[1][6][10];//U_E[i] = S[i]^T*E_constr_A[i];
		double K[6][1][10];//K[i] = E_constr_A[i]^T*S[i]
		double u[10];//vector u[i] = torques(i) - S[i]^T*(p_A[i] + I_A[i]*C[i])
		double v_constr[6][1][10];//v = M_0^-1*(beta_N - E_constr_a^T*a[0] - G_constr[0])(extra dimension is for different constraints on other segments)
		double M_0_inverse[6][6];
		double beta_N[6][1][10];//(extra dimension is for different constraints on other segments)
		double v_constr_sum[6][1];
		double c_sum_a[6][1];
		double acc[6][1][10];
		double qdotdot_sum[6][1];
		
    };
}

#endif
