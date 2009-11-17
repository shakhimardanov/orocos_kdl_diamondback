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

#include "chainidsolver_constraint_vereshchagin.hpp"
#include "frames_io.hpp"

namespace KDL{
    
    ChainIdSolver_Constraint_Vereshchagin::ChainIdSolver_Constraint_Vereshchagin(const Chain& chain_,Vector root_acc):
        chain(chain_),nj(chain.getNrOfJoints()),ns(chain.getNrOfSegments()),
        X(ns),S(ns),v(ns),a(ns),f(ns)
    {
        acc_root=-Twist(root_acc,Vector::Zero());
    }


    int ChainIdSolver_Constraint_Vereshchagin::CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques)
    {
        //Check sizes when in debug mode
        if(q.rows()!=nj || q_dot.rows()!=nj || q_dotdot.rows()!=nj || torques.rows()!=nj || f_ext.size()!=ns)
            return -1;
        unsigned int j=0;
        
        
        
		//    EXAMPLE

		int k_n = 5;//example

		//Example: only vertical acceleration !!!!
		//Declaration of alfa_N
		for(int l=0;l<k_n;l++)
		{
			for(int k=0;k<k_n;k++)
			{
				alfa_N[l][k] = 0;
			}
		}
		alfa_N[0][0] = 1;
		alfa_N[1][1] = 1;
		alfa_N[3][2] = 1;
		alfa_N[4][3] = 1;
		alfa_N[5][4] = 1;

		//Declaration of beta_N, depends on k_n
		for(int l=0;l<k_n;l++)
		{
	    	beta_N[l][0] = 0;
		}
		
		
		
		
		//Sweep from root to leaf
        for(unsigned int i=0;i<ns;i++)
        {
			ChainIdSolver_Constraint_Vereshchagin::calc_one(i, j, q, q_dot, q_dotdot, f_ext);
        }
	

        for(int i=ns-1;i>=0;i--)
        {
        	ChainIdSolver_Constraint_Vereshchagin::calc_two(i, k_n);
        }

		ChainIdSolver_Constraint_Vereshchagin::constraint_calc(k_n);

		for(int i=0;i<ns;i++)
		{
	    	ChainIdSolver_Constraint_Vereshchagin::calc_three(i, j, k_n, q_dotdot, torques);
		}

    }

    //Function for matrix inverse calculation, Determinant

    double ChainIdSolver_Constraint_Vereshchagin::determ(double num[6][6],double k)
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
            	det=det+s*(num[0][p]*ChainIdSolver_Constraint_Vereshchagin::determ(b,k-1));
            	s=-1*s;
            }
    	}
    	return det;
    }
    
    void ChainIdSolver_Constraint_Vereshchagin::calc_one(int i, int j, const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext)
	{
		double q_,qdot_,qdotdot_;
    	if(chain.getSegment(i).getJoint().getType()!=Joint::None)
    	{
        	q_=q(j);
        	qdot_=q_dot(j);
        	j++;
    	}
    	else
        	q_=qdot_=0.0;
            
    	//Calculate segment properties: X,S,vj,cj
    	X[i]=chain.getSegment(i).pose(q_);//Remark this is the inverse of the 
                                        //frame for transformations from 
                                        //the parent to the current coord frame
    	//Transform velocity and unit velocity to segment frame
    	Twist vj=X[i].M.Inverse(chain.getSegment(i).twist(q_,qdot_));
    	S[i]=X[i].M.Inverse(chain.getSegment(i).twist(q_,1.0));
    	//We can take cj=0, see remark section 3.5, page 55 since the unit velocity vector S of our joints is always time constant

 		//Transformation matrix X_matrix[i] and X_matrix_inv[i] (X_matrix[i]*I[i]*X_matrix_inv[i])
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
    	if(i==0){
        	v[i]=vj;
    	}else{
        	v[i]=X[i].Inverse(v[i-1])+vj;
    	}

    	//c[i] = cj + v[i]xvj (remark: cj=0)
   		c[i] = v[i]*vj;//This is a cross product
	    
    	//Calculate the force for the joint
    	//Collect RigidBodyInertia and external forces
    	RigidBodyInertia Ii=chain.getSegment(i).getInertia();
	    
   		//Make Inertia matrix out of RigidBodyInertia
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
	   	for(int l=0;l<6;l++)
	   	{
			p_A[l][0][i] = p_A_wrench[i](l);
	   	}
	}
	
	void ChainIdSolver_Constraint_Vereshchagin::calc_two(int i, int k_n)
	{
	    //U[i] = I_A[i]*S[i];
    	for(int l=0;l<6;l++)
    	{
			U[l][0][i] = 0;
			for(int k=0;k<6;k++)
			{
	    		U[l][0][i] += I_A[l][k][i]*S[i](k);
			}
    	}

    	//D[i] = S[i]^T*U[i]    
    	D[i] = 0;//We can add the inertia of the motor here.
  	 	for(int k=0;k<6;k++)
  	 	{
			D[i] += S[i](k)*U[k][0][i];
    	}

    	//E[i] = I_A*c[i]
    	for(int l=0;l<6;l++)
    	{
        	E[l][0][i] = 0;
        	for(int m=0;m<6;m++)
        	{
    	    	E[l][0][i] += I_A[l][m][i]*c[i](m);
        	}
    	}

    	//u[i] = torques(i) - S[i]^T*(p_A[i] + I_A[i]*C[i])
    	u[i] = 0;//For inverse dynamics, torques = 0
    	for(int k=0;k<6;k++)
    	{
			u[i] -= S[i](k)*(p_A[k][0][i] + E[k][0][i]);
    	}
   
    	if(i != 0)
    	{
			//I_a = I_A[i] - U[i]*D[i]^-1*U[i]^T
			//I_a is a usable variable, it isn`t needed out of for loop(i)
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
			//p_a is a usable variable, it isn`t needed out of for loop(i)

			//p_a = p_A[i] + E[i] + U[i]*D[i]^-1*u[i]
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
		    			I_A[l][k][i-1] += X_matrix[l][m][i]*G[m][k][i];
					}
	    		}
			}
		
			//p_A[i-1] = p_A[i-1] + X[i]*p_a
			for(int l=0;l<6;l++)
			{
				for(int m=0;m<6;m++)//sum
				{
		    		p_A[l][0][i-1] += X_matrix[l][m][i]*p_a[m][0];
				}
			}


			//             CONSTRAINTS


			//E_constr_A[ns-1] = X[ns-1]*alfa_N
			if(i=ns-1)
			{
	    		for(int l=0;l<6;l++)
	    		{
					for(int m=0;m<k_n;m++)
					{
		    			E_constr_A[l][m][i] = 0;
		    			for(int k=0;k<6;k++)//sum
		    			{
							E_constr_A[l][m][i] += X_matrix[l][k][i]*alfa_N[k][m];
		    			}
					}
	    		}
			}

			//U_E[i] = S[i]^T*E_constr_A[i];
			for(int l=0;l<k_n;l++)
			{
	    		U_E[0][l][i] = 0;
	    		for(int k=0;k<6;k++)
	    		{
					U_E[0][l][i] += S[i](k)*E_constr_A[k][l][i];
	    		}
			}


			//E_constr_a[i] = E_constr_A[i] - U[i]*D[i]^-1*U_E[i]
			for(int l=0;l<6;l++)
			{
	    		for(int k=0;k<k_n;k++)
	    		{
					E_constr_a[l][k][i] = E_constr_A[l][k][i] - (U[l][0][i]*U_E[0][k][i])/D[i];
	    		}
			}

			//E_constr_A[i-1] = X[i]*E_constr_a[i];
			for(int l=0;l<6;l++)
			{
	    		for(int k=0;k<k_n;k++)
	    		{
					E_constr_A[l][k][i-1] = 0;
					for(int m=0;m<6;m++)//sum
					{
		    			E_constr_A[l][k][i-1] += X_matrix[l][m][i]*E_constr_a[m][k][i];
					}
	    		}
			}

			//M_constr[i-1] = M_constr[i] - E_constr_A[i]^T*S[i]*D[i]^-1*S[i]^T*E_constr_A[i]
		
			//K[i] = E_constr_A[i]^T*S[i]
			for(int l=0;l<k_n;l++)
			{
	    		K[l][0][i] = 0;
	    		for(int k=0;k<6;k++)
	    		{
					K[l][0][i] += E_constr_A[l][k][i]*S[i](k);
	    		}
			}
		
			//M_constr[i-1] = M_constr[i] - K[i]*D[i]^-1*K[i]^T
			for(int l=0;l<k_n;l++)
			{
	    		for(int k=0;k<k_n;k++)
	    		{
					if(i=ns-1)
					{
		    			M_constr[l][k][i] = 0;
					}
					M_constr[l][k][i-1] = M_constr[l][k][i] - (K[l][0][i]*K[0][k][i])/D[i];
	    		}
			}

			//G_constr[i-1] = G_constr[i] + E_constr[i]^T*c[i] + E_constr[i]^T*S[i]*u[i]/D[i]
			for(int l=0;l<k_n;l++)
			{
	    		if(i=ns-1)
	    		{
					G_constr[l][0][i] = 0;
	    		}
	    		G_constr[l][0][i-1] = G_constr[l][0][i];
	    		for(int k=0;k<6;k++)
	    		{
					G_constr[l][0][i-1] += E_constr_A[k][l][i]*c[i](k);
	    		}
	    		G_constr[l][0][i-1] += K[l][0][i]*(u[i]/D[i]);
			}
    	}
    }
    
	void ChainIdSolver_Constraint_Vereshchagin::constraint_calc(int k_n)
	{
		//v_constr = M_0_inverse*(beta_N - E_constr_a^T*a[0] - G_constr[0])

		//M_0_inverse, always k_n*k_n matrix
		//Optimalisation possible !!

  		double a_constr[6][6];
		double d;

		for(int l=0;l<k_n;l++)
		{
	    	for(int m=0;m<k_n;m++)
	    	{
				a_constr[l][m] = M_constr[l][m][0];
	    	}
		}

        d = ChainIdSolver_Constraint_Vereshchagin::determ(a_constr,k_n);
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

		    if(k_n == 1)
		    {
				M_0_inverse[0][0] = 1/a_constr[0][0];
		    }
		    else
		    {
    	        for(int q=0;q<k_n;q++)
    	        {
        	    	for(int p=0;p<k_n;p++)
        	    	{
            	        int m=0;
                        int n=0;
            	        for(int l=0;l<k_n;l++)
            	        {
                	    	for(int r=0;r<k_n;r++)
                	    	{
                    	        b[l][r]=0;
                                if(l!=q&&r!=p)
                    	        {
                        	    	b[m][n]=a_constr[l][r];
                        	    	if(n<(k_n-2))
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
	            	        fac[q][p]=pow*ChainIdSolver_Constraint_Vereshchagin::determ(b,k_n-1);
	        	    }
	    	    }
	
	    	    for(int l=0;l<k_n;l++)
	    	    {
	        	    for(int r=0;r<k_n;r++)
	        	    {
	            	        b_fac[l][r]=fac[r][l];
	        	    }
	    	    }
	
		        for(int l=0;l<k_n;l++)
		        {
			    	for(int r=0;r<k_n;r++)
			    	{
			        	M_0_inverse[l][r] = b_fac[l][r]/d;
			    	}
		        }
		    }
		}

		
		//v_constr_sum = beta_N - E_constr_a^T*acc[0] - G_constr[0]
		for(int l=0;l<k_n;l++)
		{
	    	v_constr_sum[l][0] = beta_N[l][0] - G_constr[l][0][0];
	    	for(int k;k<6;k++)
	    	{
				v_constr_sum[l][0] -= E_constr_a[k][l][0]*acc_root(k);
	    	}
		}

		//v_constr = M_0_inverse*v_constr_sum
		for(int l=0;l<k_n;l++)
		{
	    	v_constr[l][0] = 0;
	    	for(int k=0;k<k_n;k++)
	    	{
				v_constr[l][0] += M_0_inverse[l][k]*v_constr_sum[k][0];
	    	}
		}	
	}    

	void ChainIdSolver_Constraint_Vereshchagin::calc_three(int i, int j, int k_n, JntArray &q_dotdot, JntArray &torques)
	{
		//Calculation of joint and segment accelerations
    	//qdotdot[i] = D[i]-1(Q[i] - S[i]^T(p_A[i] + I_A[i](c[i] + a[i-1]) + E_constr_A[i]*v_constr))
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
	    		c_sum_a[l][0] = c[i](l) + acc[l][0][i-1];
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
			for(int k=0;k<k_n;k++)
			{
	    		qdotdot_sum[l][0] += E_constr_A[l][k][i]*v_constr[k][0];
			}
    	}
    	//qdotdot[i] = (Q[i] - S[i]^T*qdotdot_sum)/D[i]
    	q_dotdot(i) = 0/D[i];//Q[i] = 0
    	for(int l=0;l<6;l++)
    	{
			q_dotdot(i) -= (S[i](l)*qdotdot_sum[l][0])/D[i];
    	}
    
    	//acc[i] = X_matrix_inv[i]*acc[i-1] + S[i]*qdotdot[i] + c[i]
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
					acc[l][0][i] += X_matrix_inv[l][m][i]*acc[m][0][i-1];
	    		}
			}
    	}
	
    	//Torques needed for constraint
    	//Forces in the joints
    	for(int k=0;k<6;k++)
    	{
			F[k][0][i] = 0;
			for(int m=0;m<k_n;m++)
			{
	    		F[k][0][i] += E_constr_A[k][m][i]*v_constr[m][0];
			}
    	}

    	//Torques
    	if(chain.getSegment(i).getJoint().getType()!=Joint::None)
    	{
        	torques(j) = 0;
        	for(int k=0;k<6;k++)
        	{
            	torques(j) += S[i](k)*F[k][0][i];
				j--;
        	}
    	}	
	}

}//namespace
