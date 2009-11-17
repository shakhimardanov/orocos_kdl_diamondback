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

#include "treeconstraints.hpp"
#include "frames_io.hpp"

namespace KDL{
    
    Tree_c::Tree_c(const Tree& tree_,Vector root_acc):
        tree(tree_),nj(tree.getNrOfJoints()),ns(tree.getNrOfSegments()),
        X(ns),S(ns),v(ns),a(ns),f(ns)
    {
        acc_root=-Twist(root_acc,Vector::Zero());//acceleration of the root, gravity is added here.
    }


    int Tree_c::CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques)
    {
    	
    	/**
     	* \brief 
     	* 
     	* Dynamics calculations by constraints based on Vereshchagin 1989 (Modelling and Control of Motion of Manipulational Robots).
     	* Code also based on Roy Featherstone, Rigid body dynamics Algorithms page 132 and 127.
     	* 
     	* Expanded to Tree structures.
     	*/
    	
    	
    	
        //Check sizes when in debug mode
        if(q.rows()!=nj || q_dot.rows()!=nj || q_dotdot.rows()!=nj || torques.rows()!=nj || f_ext.size()!=ns)
            return -1;
        unsigned int j=0;
                
        //		EXAMPLE
        //two constraints, both may only accelerate in vertical direction
        
        //All constraint segments need a defined string
        SegmentMap::const_iterator constraint_A_Iterator = tree.getSegment("constraint_A");
        SegmentMap::const_iterator constraint_B_Iterator = tree.getSegment("constraint_B");
        const TreeElement& constraint_A_Element = constraint_A_Iterator->second;
        const TreeElement& constraint_B_Element = constraint_B_Iterator->second;
        //If there is a new TreeElement with all the variables (header file), this isn`t needed
        int q_constraint_A = constraint_A_Element.q_nr;
        int q_constraint_B = constraint_B_Element.q_nr;
   
        
        
        //constraints dependent
        int k_total=10;//The total amount of constraints on all segments
        
        
        //constraints dependent
        int k_n[10];
        for(int i=0;i<10;i++)
        {
        	k_n[i]=0;
        }
        
        //Constraint A
        k_n[q_constraint_A] = 5;//example

		//Example: only vertical acceleration !!!!
		//Declaration of alfa_N
		for(int l=0;l<k_n[q_constraint_A];l++)
		{
			for(int k=0;k<k_n[q_constraint_A];k++)
			{
				alfa_N[l][k][q_constraint_A] = 0;
			}
		}
		alfa_N[0][0][q_constraint_A] = 1;
		alfa_N[1][1][q_constraint_A] = 1;
		alfa_N[3][2][q_constraint_A] = 1;
		alfa_N[4][3][q_constraint_A] = 1;
		alfa_N[5][4][q_constraint_A] = 1;

		//Declaration of beta_N, depends on k_n
		for(int l=0;l<k_n[q_constraint_A];l++)
		{
	    	beta_N[l][0][q_constraint_A] = 0;
		}
		
		
		//Constraint B
        k_n[q_constraint_B] = 5;//example

		//Example: only vertical acceleration !!!!
		//Declaration of alfa_N
		for(int l=0;l<k_n[q_constraint_B];l++)
		{
			for(int k=0;k<k_n[q_constraint_B];k++)
			{
				alfa_N[l][k][q_constraint_B] = 0;
			}
		}
		alfa_N[0][0][q_constraint_B] = 1;
		alfa_N[1][1][q_constraint_B] = 1;
		alfa_N[3][2][q_constraint_B] = 1;
		alfa_N[4][3][q_constraint_B] = 1;
		alfa_N[5][4][q_constraint_B] = 1;

		//Declaration of beta_N, depends on k_n
		for(int l=0;l<k_n[q_constraint_B];l++)
		{
	    	beta_N[l][0][q_constraint_B] = 0;
		}
             
        
        /*
         * All variables are stored in arrays, the location of the right variable of the children/parent is done by an iterator of the tree structure for each segment.
         */
        
        
        /*
         * For the calculation of different constraints of a tree, every segments needs to be passed three times.
         * 1.a)calculation of velocity`s (Up direction, from root to leaf)
         * 1.b)calculation of propagated inertia`s, propagated bias forces and constraint forces
         * 2)calculation of joint and segment accelerations
         */
         
        // Recursive algorithm, from root to leaf and back
        // 1)
        SegmentMap::const_iterator rootIterator = tree.getSegment("root");
		Tree_c::recursiveDepthFirst(rootIterator, j, k_n, k_total, q, q_dot, q_dotdot, f_ext, torques);//calculation of velocity, Propagated inertia, Propagated bias forces and constraint forces
		
		// Recursive algorithm, from root to leaf
		// 2)
		Tree_c::recursiveUP(rootIterator, j, k_n, k_total, q, q_dot, q_dotdot, f_ext, torques);//calculation of accelerations


    }

	// 1)
	void Tree_c::recursiveDepthFirst(const SegmentMap::const_iterator& it, int j, int k_n[10], int k_total, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques)//calculation of propagation inertia`s and constraints
	{
		const TreeElement& currentElement = it->second;		
		int nc = currentElement.children.size();//Number of Children
		
		calc_one(it, j, q, q_dot, q_dotdot, f_ext);//Calculation root to leaf
		for(int i=0;i<nc;i++)
		{
			SegmentMap::const_iterator childrenIt = currentElement.children[i];
			recursiveDepthFirst(childrenIt, j, k_n, k_total, q, q_dot, q_dotdot, f_ext, torques);
		}
		calc_two(it, k_n);//Calculation leaf to root
		SegmentMap::const_iterator rootIterator = tree.getSegment("root");
		if(it == rootIterator)
		{
			constraint_calc(k_total, k_n);//Calculation of v(constraint magnitude)
		}
		
	}
	
	// 2)
	void Tree_c::recursiveUP(const SegmentMap::const_iterator& it, int j, int k_n[10], int k_total, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques)//calculation of accelerations
	{
		const TreeElement& currentElement = it->second;		
		int nc = currentElement.children.size();//Number of Children
		
		calc_three(it, j, k_n, k_total, q_dotdot, torques);
		for(int i=0;i<nc;i++)
		{
			SegmentMap::const_iterator childrenIt = currentElement.children[i];
			recursiveUP(childrenIt, j, k_n, k_total, q, q_dot, q_dotdot, f_ext, torques);
		}
	}


    //Function for matrix inverse calculation, Determinant
    //Will be deleted when svd is used
    double Tree_c::determ(double num[6][6],double k)
    {
		double s=1;
		double det=0;
		double b[6][6];

		if(k==1)
		{
	    	return(num[0][0]);
		}
		else
    	{
            det=0;
            for(int p=0;p<k;p++)
            {
                int m=0;
            	int n=0;
            	for(int l=0;l<k;l++)
            	{
                    for(int r=0;r<k;r++)
                    {
                    	b[l][r]=0;
                    	if(l!=0&&r!=p)
                    	{
                            b[m][n]=num[l][r];
                            if(n<(k-2))
						    {
                            	n++;
			    			}
                            else
                            {
                            	n=0;
                            	m++;
                            }
                    	}
                    }
            	}
            	det=det+s*(num[0][p]*Tree_c::determ(b,k-1));
            	s=-1*s;
            }
    	}
    	return det;
    }
    
    //Calculate velocity, transformation matrices and segment inertia matrices
    void Tree_c::calc_one(const SegmentMap::const_iterator& it, int j, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext)
	{
		double q_,qdot_;
		const TreeElement& currentElement = it->second;
    	if(currentElement.segment.getJoint().getType()!=Joint::None)
    	{
        	q_=q(currentElement.q_nr);
        	qdot_=q_dot(currentElement.q_nr);
        	j++;
    	}
    	else
        	q_=qdot_=0.0;
        	
        //i = the place in the array of the current segment	
        int i = currentElement.q_nr;
        //i_1 = the place in the array of the parent segment
        SegmentMap::const_iterator parentIt = currentElement.parent;
        const TreeElement& parentElement = parentIt->second;
        int i_1 = parentElement.q_nr;	
            
    	//Calculate segment properties: X,S,vj,cj
    	//(Roy Featherstone p132)
    	X[i]=currentElement.segment.pose(q_);//Remark this is the inverse of the 
                                        //frame for transformations from 
                                        //the parent to the current coord frame
    	//Transform velocity and unit velocity to segment frame
    	Twist vj=X[i].M.Inverse(currentElement.segment.twist(q_,qdot_));
    	S[i]=X[i].M.Inverse(currentElement.segment.twist(q_,1.0));
    	//We can take cj=0, see remark section 3.5, page 55 (Featherstone) since the unit velocity vector S of our joints is always time constant

 		//Transformation matrix X_matrix[i] and X_matrix_inv[i] (X_matrix[i]*I[i]*X_matrix_inv[i]), for and after transformation matrices
 		//Optimalisation possible !!
 		//Make an operator of this
   		for(int l=0;l<3;l++)
   		{
			for(int k=0;k<3;k++)
			{
    			X_matrix[l][k][i] = X[i].M(l,k);
   				X_matrix_inv[l][k][i] = X[i].M(k,l);
			}
			for(int k=3;k<6;k++)
			{
    			X_matrix[l][k][i] = 0;
			}
			X_matrix_inv[0][l+3][i] = -X[i].p(2)*X[i].M(l,1) + X[i].p(1)*X[i].M(l,2);
			X_matrix_inv[1][l+3][i] = X[i].p(2)*X[i].M(l,0) - X[i].p(0)*X[i].M(l,2);
			X_matrix_inv[2][l+3][i] = -X[i].p(1)*X[i].M(l,0) + X[i].p(0)*X[i].M(l,1);
   		}
   		for(int l=3;l<6;l++)
   		{
			for(int k=0;k<3;k++)
			{
    			X_matrix_inv[l][k][i] = 0;
			}
			for(int k=3;k<6;k++)
			{
    			X_matrix[l][k][i] = X[i].M(l-3,k-3);
    			X_matrix_inv[l][k][i] = X[i].M(k-3,l-3);
			}
			X_matrix[l][0][i] = -X[i].M(l-3,1)*X[i].p(2) + X[i].M(l-3,2)*X[i].p(1);
       		X_matrix[l][1][i] = X[i].M(l-3,0)*X[i].p(2) - X[i].M(l-3,2)*X[i].p(0);
       		X_matrix[l][2][i] = -X[i].M(l-3,0)*X[i].p(1) + X[i].M(l-3,1)*X[i].p(0);
   		}

    	//calculate velocity of the segment (in segment coordinates)
    	//(Roy Featherstone p132)	    
    	if(i==0){//root
        	v[i]=vj;
    	}else{
        	v[i]=X[i].Inverse(v[i_1])+vj;
    	}

    	//c[i] = cj + v[i]xvj (remark: cj=0)
    	//(Roy Featherstone p132)
   		c[i] = v[i]*vj;//This is a cross product
	    
    	//Calculate the force for the joint
    	//Collect RigidBodyInertia and external forces
    	//Possible: make a class, ArticulatedBodyInertia
    	RigidBodyInertia Ii=currentElement.segment.getInertia();
	    
   		//Make Inertia matrix out of RigidBodyInertia
   		//(Roy Featherstone p132)
   		//Optimalisation possible !! Inertia matrix is a symmetric matrix
   		for(int l=0;l<3;l++)
   		{
			for(int m=0;m<3;m++)
			{
    			if(l==m)
    			{
					I_A[l][m][i] = Ii.m;
    			}
    			else
    			{
					I_A[l][m][i] = 0;
				}
			}
   		}
   		int r=0;
   		for(int l=3;l<6;l++)
   		{
			for(int m=3;m<6;m++)
			{
    			I_A[l][m][i] = Ii.I.data[r];
    			r++;
			}
   		}
   		I_A[0][3][i] = I_A[3][0][i] = 0;
   		I_A[0][4][i] = I_A[4][0][i] = Ii.m*Ii.h(3);
   		I_A[0][5][i] = I_A[5][0][i] = -Ii.m*Ii.h(2);
   		I_A[1][3][i] = I_A[3][1][i] = -Ii.m*Ii.h(3);
   		I_A[1][4][i] = I_A[4][1][i] = 0;
	   	I_A[1][5][i] = I_A[5][1][i] = Ii.m*Ii.h(1);
	   	I_A[2][3][i] = I_A[3][2][i] = Ii.m*Ii.h(2);
	   	I_A[2][4][i] = I_A[4][2][i] = -Ii.m*Ii.h(1);
	   	I_A[2][5][i] = I_A[5][2][i] = 0;
	
	   	//wrench p of the bias forces
	   	//external forces on the segments (in body coordinates)
	   	p_A_wrench[i] = v[i]*(Ii*v[i]) - f_ext[i];
	   
	   	//Change p_A_wrench[i] to matrix form
	   	//Can be deleted and changed further down
	   	for(int l=0;l<6;l++)
	   	{
			p_A[l][0][i] = p_A_wrench[i](l);
	   	}
	}
	
	//Calculate propagated inertia, propagated bias forces and constraint forces
	void Tree_c::calc_two(const SegmentMap::const_iterator& it, int k_n[10])
	{
		//i = the place in the array of the current segment	
		const TreeElement& currentElement = it->second;
		int i = currentElement.q_nr;
		//i_1 = the place in the array of the parent segment
        SegmentMap::const_iterator parentIt = currentElement.parent;
        const TreeElement& parentElement = parentIt->second;
        int i_1 = parentElement.q_nr;
        
        
	    //U[i] = I_A[i]*S[i];
	    //(Roy Featherstone p132)
	    //local variable (no need to store it out of function)
    	for(int l=0;l<6;l++)
    	{
			U[l][0][i] = 0;
			for(int k=0;k<6;k++)
			{
	    		U[l][0][i] += I_A[l][k][i]*S[i](k);
			}
    	}

    	//D[i] = S[i]^T*U[i]    
    	//(Roy Featherstone p132)
    	D[i] = 0;//We can add the inertia of the motor here.
  	 	for(int k=0;k<6;k++)
  	 	{
			D[i] += S[i](k)*U[k][0][i];
    	}

    	//E[i] = I_A*c[i]
    	//(Roy Featherstone p132)
    	//local variable (no need to store it out of function)
    	for(int l=0;l<6;l++)
    	{
        	E[l][0][i] = 0;
        	for(int m=0;m<6;m++)
        	{
    	    	E[l][0][i] += I_A[l][m][i]*c[i](m);
        	}
    	}

    	//u[i] = torques(i) - S[i]^T*(p_A[i] + I_A[i]*C[i])
    	//(Roy Featherstone p127)
    	u[i] = 0;//For inverse dynamics, torques = 0
    	for(int k=0;k<6;k++)
    	{
			u[i] -= S[i](k)*(p_A[k][0][i] + E[k][0][i]);
    	}
   
    	if(i != 0)
    	{
			//I_a = I_A[i] - U[i]*D[i]^-1*U[i]^T
			//(Roy Featherstone p132)
			//local variable (no need to store it out of function)
			for(int l=0;l<6;l++)
			{
	    		for(int k=l;k<6;k++)
	    		{
					I_a[l][k] = I_A[l][k][i] - (U[l][0][i]*U[k][0][i])/D[i];
	    		}
			}
			//declaration of symmetry
			//optimalisation possible !!
			for(int l=0;l<6;l++)
			{
	    		for(int k=0;k<6;k++)
	    		{
					I_a[k][l] = I_a[l][k];
	    		}
			}

			//p_a = p_A[i] + I_A*c[i] + U[i]*D[i]^-1*u[i]
			//p_a = p_A[i] + E[i] + U[i]*D[i]^-1*u[i]
			//(Roy Featherstone p132)
			//local variable (no need to store it out of function)
			for(int l=0;l<6;l++)
			{
	    		p_a[l][0] = p_A[l][0][i] + E[l][0][i] + U[l][0][i]*(u[i]/D[i]);
			}

			//I_A[i-1] = I_A[i-1] + X[i]*I_a*X[i].inverse()
			//Transformation matrix X_matrix[i] and X_matrix_inv[i]

			//G[i] = I_a*X_matrix_inv[i]
			//Optimalisation possible !!
			for(int l=0;l<6;l++)
			{
	    		for(int k=0;k<6;k++)
	    		{
					G[l][k][i] = 0;
					for(int m=0;m<6;m++)//sum
					{
		    			G[l][k][i] += I_a[l][m]*X_matrix_inv[m][k][i];
					}
	    		}
			}

			//I_A[i-1] = I_A[i-1] + X[i]*G[i]
			//Optimalisation possible !!
			for(int l=0;l<6;l++)
			{
	    		for(int k=0;k<6;k++)
	    		{
					for(int m=0;m<6;m++)//sum
					{
		    			I_A[l][k][i_1] += X_matrix[l][m][i]*G[m][k][i];
					}
	    		}
			}
		
			//p_A[i-1] = p_A[i-1] + X[i]*p_a
			//(Roy Featherstone p132)
			for(int l=0;l<6;l++)
			{
				for(int m=0;m<6;m++)//sum
				{
		    		p_A[l][0][i_1] += X_matrix[l][m][i]*p_a[m][0];
				}
			}


			//             									CONSTRAINTS
			
			
			
			
			
			/* 
			 * 
			 * 
			 * 												NOT FINISHED:
			 * 
			 * 
			 * 
			 * 
			 * When two constraints join together in a joint, the propagated matrix E_constr_A needs to be updated to the new one.
			 * 
			 * 
			 * It needs to update to:
			 * 
			 * E_constr_A_new = [E_constr_A_old1 E_constr_A_old2] 
			 * M_constr_new = [ M_constr_old1     0
			 * 					 	 0         M_constr_old2 ]
			 * 
			 * 
			 * 
			 * Possible solution: check -> if(E_constr_A !=0 for 2 or more children) then make E_constr_A_new, M_constr_new and work with these ones further
			 * 								else just go further
			 * 
			 * 
			 * 
			 * 
			 */


			for(int constr=0;constr<ns;constr++)
			{
				if(k_n[constr] != 0)
				{
					//E_constr_A[ns-1] = X[ns-1]*alfa_N
					//(Vereshchagin 1989, point c)
					if(i=constr)
					{
			    		for(int l=0;l<6;l++)
			    		{
							for(int m=0;m<k_n[constr];m++)
							{
				    			E_constr_A[l][m][i][constr] = 0;
				    			for(int k=0;k<6;k++)//sum
				    			{
									E_constr_A[l][m][i][constr] += X_matrix[l][k][i]*alfa_N[k][m][constr];
				    			}
							}
			    		}
					}
		
					//U_E[i] = S[i]^T*E_constr_A[i];
					//(Vereshchagin 1989, point c)
					for(int l=0;l<k_n[constr];l++)
					{
			    		U_E[0][l][i] = 0;
			    		for(int k=0;k<6;k++)
			    		{
							U_E[0][l][i] += S[i](k)*E_constr_A[k][l][i][constr];
			    		}
					}
		
		
					//E_constr_a[i] = E_constr_A[i] - U[i]*D[i]^-1*U_E[i]
					//(Vereshchagin 1989, point c)
					for(int l=0;l<6;l++)
					{
			    		for(int k=0;k<k_n[constr];k++)
			    		{
							E_constr_a[l][k][i][constr] = E_constr_A[l][k][i][constr] - (U[l][0][i]*U_E[0][k][i])/D[i];
			    		}
					}
		
					//E_constr_A[i-1] = X[i]*E_constr_a[i];
					//(Vereshchagin 1989, point c)
					for(int l=0;l<6;l++)
					{
			    		for(int k=0;k<k_n[constr];k++)
			    		{
							E_constr_A[l][k][i_1][constr] = 0;
							for(int m=0;m<6;m++)//sum
							{
				    			E_constr_A[l][k][i_1][constr] += X_matrix[l][m][i]*E_constr_a[m][k][i][constr];
							}
			    		}
					}
		
					//M_constr[i-1] = M_constr[i] - E_constr_A[i]^T*S[i]*D[i]^-1*S[i]^T*E_constr_A[i]
					//(Vereshchagin 1989, point d)
				
					//K[i] = E_constr_A[i]^T*S[i]
					for(int l=0;l<k_n[constr];l++)
					{
			    		K[l][0][i] = 0;
			    		for(int k=0;k<6;k++)
			    		{
							K[l][0][i] += E_constr_A[l][k][i][constr]*S[i](k);
			    		}
					}
				
					//M_constr[i-1] = M_constr[i] - K[i]*D[i]^-1*K[i]^T
					//(Vereshchagin 1989, point d)
					for(int l=0;l<k_n[constr];l++)
					{
			    		for(int k=0;k<k_n[constr];k++)
			    		{
							if(i=constr)
							{
				    			M_constr[l][k][i][constr] = 0;
							}
							M_constr[l][k][i_1][constr] = M_constr[l][k][i][constr] - (K[l][0][i]*K[0][k][i])/D[i];
			    		}
					}
		
					//G_constr[i-1] = G_constr[i] + E_constr[i]^T*c[i] + E_constr[i]^T*S[i]*u[i]/D[i]
					//(Vereshchagin 1989, point e)
					for(int l=0;l<k_n[constr];l++)
					{
			    		if(i=constr)
			    		{
							G_constr[l][0][i][constr] = 0;
			    		}
			    		G_constr[l][0][i_1][constr] = G_constr[l][0][i][constr];
			    		for(int k=0;k<6;k++)
			    		{
							G_constr[l][0][i_1][constr] += E_constr_A[k][l][i][constr]*c[i](k);
			    		}
			    		G_constr[l][0][i_1][constr] += K[l][0][i]*(u[i]/D[i]);
					}
				}
			}
    	}
    }
    
    //Calculate the magnitude of the constraint forces
    //(Vereshchagin 1989, point f)
	void Tree_c::constraint_calc(int k_total, int k_n[10])
	{

		for(int constr=0;constr<ns;constr++)//This has to go when a new tree segment is used !!
		{
			if(k_n[constr] != 0)
			{
				//v_constr = M_0_inverse*(beta_N - E_constr_a^T*a[0] - G_constr[0])
		
				//M_0_inverse, always k_n*k_n matrix
				//Optimalisation possible !!
		
		  		double a_constr[6][6];
				double d;
		
				for(int l=0;l<k_total;l++)
				{
			    	for(int m=0;m<k_total;m++)
			    	{
						a_constr[l][m] = M_constr[l][m][0][constr];
			    	}
				}
		
		        d = Tree_c::determ(a_constr,k_total);
		   		if(d<=1/100000000)
				{
			    	//error !! Not inversible
				}
		  		else
				{
			    	double b[6][6];
			    	double fac[6][6];
			    	double pow;
		    	    double b_fac[6][6];
		
				    if(k_total == 1)
				    {
						M_0_inverse[0][0] = 1/a_constr[0][0];
				    }
				    else
				    {
		    	        for(int q=0;q<k_total;q++)
		    	        {
		        	    	for(int p=0;p<k_total;p++)
		        	    	{
		            	        int m=0;
		                        int n=0;
		            	        for(int l=0;l<k_total;l++)
		            	        {
		                	    	for(int r=0;r<k_total;r++)
		                	    	{
		                    	        b[l][r]=0;
		                                if(l!=q&&r!=p)
		                    	        {
		                        	    	b[m][n]=a_constr[l][r];
		                        	    	if(n<(k_total-2))
		                            	        n++;
		                        	    	else
		                        	    	{
		                            	        n=0;
		                            	        m++;
		                        	    	}
		                    	    	}
		                	    	}
		            	        }
						        pow = 1;
						        for(int l=0;l<q+p;l++)
						        {
								    pow = pow*-1;
					        	}
			            	        fac[q][p]=pow*Tree_c::determ(b,k_total-1);
			        	    }
			    	    }
			
			    	    for(int l=0;l<k_total;l++)
			    	    {
			        	    for(int r=0;r<k_total;r++)
			        	    {
			            	        b_fac[l][r]=fac[r][l];
			        	    }
			    	    }
			
				        for(int l=0;l<k_total;l++)
				        {
					    	for(int r=0;r<k_total;r++)
					    	{
					        	M_0_inverse[l][r] = b_fac[l][r]/d;
					    	}
				        }
				    }
				}

		
				//v_constr_sum = beta_N - E_constr_a^T*acc[0] - G_constr[0]
				for(int l=0;l<k_total;l++)
				{
			    	v_constr_sum[l][0] = beta_N[l][0][constr] - G_constr[l][0][0][constr];
			    	for(int k;k<6;k++)
			    	{
						v_constr_sum[l][0] -= E_constr_a[k][l][0][constr]*acc_root(k);
			    	}
				}
		
				//v_constr = M_0_inverse*v_constr_sum
				for(int l=0;l<k_total;l++)
				{
			    	v_constr[l][0][constr] = 0;
			    	for(int k=0;k<k_total;k++)
			    	{
						v_constr[l][0][constr] += M_0_inverse[l][k]*v_constr_sum[k][0];
			    	}
				}
			}
		}
	}
   
	//Calculate the accelerations of the joints and segments and the torques needed in the joints
	void Tree_c::calc_three(const SegmentMap::const_iterator& it, int j, int k_n[10], int k_total, JntArray &q_dotdot, JntArray &torques)
	{
		//i = the place in the array of the current segment	
		const TreeElement& currentElement = it->second;
		int i = currentElement.q_nr;
		//i_1 = the place in the array of the parent segment
        SegmentMap::const_iterator parentIt = currentElement.parent;
        const TreeElement& parentElement = parentIt->second;
        int i_1 = parentElement.q_nr;
		
		//Calculation of joint and segment accelerations
    	//qdotdot(i) = D[i]-1(Q[i] - S[i]^T(p_A[i] + I_A[i](c[i] + a[i-1]) + E_constr_A[i]*v_constr))
    	//(Vereshchagin 1989, point g)
    	//c_sum_a = c[i] + a[i-1]
    	if(i=0)
    	{
			for(int l=0;l<6;l++)
			{
	    		c_sum_a[l][0] = c[i](l) + acc_root(l);
			}
    	}
    	else
    	{
			for(int l=0;l<6;l++)
			{
	    		c_sum_a[l][0] = c[i](l) + acc[l][0][i_1];
			}
    	}
    	//qdotdot_sum = p_A[i] + I_A[i](c_sum_a) + E_constr_A[i]*v_constr
    	for(int l=0;l<6;l++)
    	{
			qdotdot_sum[l][0] = p_A[l][0][i];
			for(int m=0;m<6;m++)
			{
	    		qdotdot_sum[l][0] += I_A[l][m][i]*c_sum_a[m][0];
			}
			for(int constr=0;constr<ns;constr++)
			{
				if(k_n[constr] != 0)
				{
					for(int k=0;k<k_total;k++)//total matrix with all the constraints and constraint multipliers in
					{
	    				qdotdot_sum[l][0] += E_constr_A[l][k][i][constr]*v_constr[k][0][constr];
					}
				}
			}
    	}
    	
    	//joint accelerations
    	//qdotdot(i) = (Q[i] - S[i]^T*qdotdot_sum)/D[i]
    	//(Vereshchagin 1989, point g)
    	q_dotdot(i) = 0/D[i];//Q[i] = 0
    	for(int l=0;l<6;l++)
    	{
			q_dotdot(i) -= (S[i](l)*qdotdot_sum[l][0])/D[i];
    	}
    
    	//accelerations of the segments
    	//acc[i] = X_matrix_inv[i]*acc[i-1] + S[i]*qdotdot(i) + c[i]
    	//(Vereshchagin 1989, point h)
    	if(i=0)
    	{
			for(int l=0;l<6;l++)
			{
	    		acc[l][0][i] = S[i](l)*q_dotdot(i) + c[i](l);
	    		for(int m=0;m<6;m++)
	    		{
					acc[l][0][i] += X_matrix_inv[l][m][i]*acc_root(m);
	    		}
			}
    	}
    	else
    	{
			for(int l=0;l<6;l++)
			{
	    		acc[l][0][i] = S[i](l)*q_dotdot(i) + c[i](l);
	    		for(int m=0;m<6;m++)
	    		{
					acc[l][0][i] += X_matrix_inv[l][m][i]*acc[m][0][i_1];
	    		}
			}
    	}
	
    	//Torques needed for constraint
    	//Forces in the joints
    	for(int k=0;k<6;k++)
    	{
			F[k][0][i] = 0;
			for(int constr=0;constr<ns;constr++)
			{
				if(k_n[constr] != 0)
				{
					for(int m=0;m<k_total;m++)//total matrix of all the constraints
					{
	    				F[k][0][i] += E_constr_A[k][m][i][constr]*v_constr[m][0][constr];
					}
    			}
			}
    	}

    	//Torques
    	if(currentElement.segment.getJoint().getType()!=Joint::None)
    	{
        	torques(i) = 0;
        	for(int k=0;k<6;k++)
        	{
            	torques(i) += S[i](k)*F[k][0][i];
				j--;//could check if j=0
        	}
    	}	
	}

}//namespace
