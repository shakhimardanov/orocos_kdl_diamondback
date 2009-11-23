b// Copyright  (C)  2009

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
#include <Eigen/SVD>

namespace KDL{
    using namespace Eigen;

    ChainIdSolver_Constraint_Vereshchagin::ChainIdSolver_Constraint_Vereshchagin(const Chain& chain_,Vector root_acc, unsigned int nr_of_constraints,const Jacobian& alfa,const JntArray& beta):
        chain(chain_),nj(chain.getNrOfJoints()),ns(chain.getNrOfSegments()),nc(nr_of_constraints),
        results(ns,segment_info(nc)),alfa_N(alfa),beta_N(beta)
    {
        acc_root=-Twist(root_acc,Vector::Zero());
        M_0_inverse.resize(nc,nc);
        assert(alfa_N.columns()!=nc);
        assert(beta_N.rows()!=nc);
       
    }


    int ChainIdSolver_Constraint_Vereshchagin::CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Wrenches& f_ext,JntArray &torques)
    {
        //Check sizes when in debug mode
        if(q.rows()!=nj || q_dot.rows()!=nj || q_dotdot.rows()!=nj || torques.rows()!=nj || f_ext.size()!=ns)
            return -1;
                
        this->initial_upwards_sweep(q, q_dot, f_ext);
        this->downwards_sweep();
		this->constraint_calculation();
        this->final_upwards_sweep(q_dotdot, torques);
    }

    void ChainIdSolver_Constraint_Vereshchagin::initial_upwards_sweep(const JntArray &q, const JntArray &qdot, const Wrenches& f_ext)
	{
        unsigned int j=0;
        for(unsigned int i=0;i<ns;++i){
            //Calculate segment properties: X,S,vj,cj
            const Segment& segment=chain.getSegment(i);
            segment_info& s=results[i];
            s.F=segment.pose(q(j));
            //frame for transformations from 
            //the parent to the current coord frame
            //Transform velocity and unit velocity to segment frame
            Twist vj=s.F.M.Inverse(segment.twist(q(j),qdot(j)));
            s.Z=s.F.M.Inverse(segment.twist(q(j),1.0));
            //We can take cj=0, see remark section 3.5, page 55 since the unit velocity vector S of our joints is always time constant
            
          
            //calculate velocity of the segment (in segment coordinates)	    
            if(i!=0)
                s.v=s.F.Inverse(results[i-1].v)+vj;
            else
                s.v=vj;
            
            //c[i] = cj + v[i]xvj (remark: cj=0, since our S is not time dependent in local coordinates)
            s.C = s.v*vj;//This is a cross product
        
            s.H=segment.getInertia();
        
            //wrench p of the bias forces
            //external forces on the segments (in body coordinates)
            s.U = s.v*(s.H*s.v) - f_ext[i];

            if(segment.getJoint().getType()!=Joint::None)
                {
                    j++;
                }
        }
    }

	void ChainIdSolver_Constraint_Vereshchagin::downwards_sweep()
	{
        for(int i=ns-1;i>0;i--)
            {
                //Get a handle for the segment we are working on.
                segment_info& s=results[i];
                //For segment N,
                if(i=ns-1){
                    s.P_tilde=s.H;
                    s.R_tilde=s.U;
                    s.M.setZero();
                    s.G.setZero(); 
                    for(unsigned int r=0;r<6;r++)
                        for(unsigned int c=0;c<nc;c++)
                            s.E_tilde(r,c)=alfa_N(r,c);
                }else{
                    //a)Pi_tilde=Hi+P(i+1) -P(i+1)*Zi/D(i+1)*Zi'*P(i+1)'
                    segment_info& parent=results[i+1];
                    Vector6d vPZ;
                    vPZ<<Vector3d::Map(parent.PZ.force.data),Vector3d::Map(parent.PZ.torque.data);
                    Matrix6d PZDZtPt;
                    PZDZtPt.part<SelfAdjoint>()= (vPZ*vPZ.transpose()).lazy()/parent.D;
                    s.P_tilde = ((s.H + parent.P) - ArticulatedBodyInertia(PZDZtPt.corner<3,3>(TopLeft),PZDZtPt.corner<3,3>(TopRight),PZDZtPt.corner<3,3>(BottomRight)));
                    
                    //b) Ri_tilde = Ui+R(i+1)+P(i+1)*C(i+1)(Q(i+1)-Zi'((Ri+1)+P(i+1)*C(i+1))
                    double torque = -dot(parent.Z,parent.R+parent.PC);//torque=Q(i+1)-Zi'(R(i+1)+P(i+1)*C(i+1)), Qi=external joint torque, 0 in our case.
                    s.R_tilde = ((s.U+parent.R)+parent.PC)+parent.PZ/parent.D*torque;
                    
                    //c) Ei_tilde = E(i+1) - P(i+1)*Zi/D(i+1)*Zi'*E(i+1)
                    Vector6d Zi;
                    Zi<<Vector3d::Map(parent.Z.rot.data),Vector3d::Map(parent.Z.vel.data);
                    Matrix6d tmp = -(vPZ*Zi.transpose()).lazy();
                    tmp/=parent.D;
                    s.E_tilde=(tmp*parent.E).lazy();
                    s.E_tilde+=parent.E;
                    
                    //d) Mi = M(i+1) - E(i+1)'*Zi/D(i+1)*Zi'*E(i+1)
                    //Ki = E(i+1)'*Zi
                    Vector6d Ki = (parent.E.transpose()*Zi).lazy();
                    s.M.part<SelfAdjoint>()=-(Ki*Ki.transpose()).lazy();
                    s.M/=parent.D;
                    s.M+=parent.M;

                    //e) Gi = G(i+1) + E(i+1)'(C(i+1) + Zi/D(i+1)*(Q(i+1)-Zi'(R(i+1)+P(i+1)*C(i+1))
                    Twist CiminZetc = parent.C+parent.Z/parent.D*torque;
                    Vector6d Ci;
                    Ci<<Vector3d::Map(CiminZetc.rot.data),Vector3d::Map(CiminZetc.vel.data);
                    s.G=(parent.E.transpose()*Ci).lazy();
                    s.G+=parent.G;
                }
                //a) Pi=F'*Pi_tilde*F
                s.P=s.F*s.P_tilde;//This is a full transform F*P*F'
                s.PC = s.P*s.C;
                s.PZ = s.P*s.Z;
                //a)Di=di+Z(i-1)'*Pi*Z(i-1)
                s.D = dot(s.Z,s.PZ);//We can add the inertia of the motor, di here.

                //b) Ri = F'Ri_tilde
                s.R=s.F*s.R_tilde;
                //c) Ei = F'Ei_tilde
                //Every column of E is a Wrench, so use Wrench representation for transformation
                for(unsigned int c=0;c<nc;c++){
                    Wrench col(Vector(s.E_tilde(0,c),s.E_tilde(1,c),s.E_tilde(2,c)),
                               Vector(s.E_tilde(3,c),s.E_tilde(4,c),s.E_tilde(5,c)));
                    col=s.F*col;
                    s.E.col(c)<<Vector3d::Map(col.force.data),Vector3d::Map(col.torque.data);
                }                                
                
            }
    }

    
	void ChainIdSolver_Constraint_Vereshchagin::constraint_calculation()
	{
		//f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)

		//M_0_inverse, always nc*nc matrix
		//Optimalisation possible !!
        results[0].M.computeInverse(&M_0_inverse);
        Vector6d acc;
        acc<<Vector3d::Map(acc_root.rot.data),Vector3d::Map(acc_root.vel.data);
        nu_sum=-(results[0].E_tilde.transpose()*acc).lazy();
        nu_sum+=beta_N.data;
        nu_sum -= results[0].G;
        nu = (M_0_inverse * nu_sum).lazy();
    }    


	void ChainIdSolver_Constraint_Vereshchagin::final_upwards_sweep(JntArray &q_dotdot, JntArray &torques){
        {
            unsigned int j=0;
            
            for(unsigned int i=0;i<ns;i++){
                segment_info& s=results[i];
                //Calculation of joint and segment accelerations
                //g) qdotdot[i] = D[i]-1(Q[i] - Z[i]^T(R[i] + P[i](C[i] + acc[i-1]) + E[i]*nu))
                Twist CAcc;
                if(i=0)
                    CAcc=s.C+acc_root;
                else
                    CAcc=s.C+results[i-1].acc;
                
                Vector6d c_f = s.E*v_constr;
                Wrench constraint_force(Vector(c_f(0),c_f(1),c_f(2)),Vector(c_f(3),c_f(4),c_f(5)));
                Wrench total_force=(s.R+s.P*CAcc)+constraint_force;
                
                q_dotdot(j) = -dot(s.Z,total_force)/s.D;//Q[i] = 0

                //h) acc[i] = Fi*(acc[i-1] + Z[i]*qdotdot[i] + c[i]
                if(i!=0)
                    s.acc=s.F*(results[i-1].acc+s.Z*q_dotdot(j)+s.C);
                
                if(chain.getSegment(i).getJoint().getType()!=Joint::None)
                    {
                        j++;
                    }
            }
        }
    }
}//namespace
