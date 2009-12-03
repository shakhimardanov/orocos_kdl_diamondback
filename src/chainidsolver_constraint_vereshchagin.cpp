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
#include "utilities/svd_eigen_HH.hpp"
#include<Eigen/Cholesky>

namespace KDL{
    using namespace Eigen;

    ChainIdSolver_Constraint_Vereshchagin::ChainIdSolver_Constraint_Vereshchagin(const Chain& chain_,Twist root_acc,unsigned int _nc):
        chain(chain_),nj(chain.getNrOfJoints()),ns(chain.getNrOfSegments()),nc(_nc),
        results(ns+1,segment_info(nc))
    {
        acc_root=root_acc;

        //Provide the necessary memory for computing the inverse of M0
        M_0_inverse.resize(nc,nc);
        Um=MatrixXd::Identity(nc,nc);
        Vm=MatrixXd::Identity(nc,nc);
        Sm=VectorXd::Ones(nc);
        tmpm=VectorXd::Ones(nc);
    }


    int ChainIdSolver_Constraint_Vereshchagin::CartToJnt(const JntArray &q, const JntArray &q_dot, JntArray &q_dotdot, const Jacobian& alfa, const JntArray& beta,const Wrenches& f_ext,JntArray &torques)
    {
        //Check sizes always
        if(q.rows()!=nj || q_dot.rows()!=nj || q_dotdot.rows()!=nj || torques.rows()!=nj || f_ext.size()!=ns)
            return -1;
        if(alfa.columns()!=nc||beta.rows()!=nc)
            return -2;

        //do an upward recursion for position and velocities
        this->initial_upwards_sweep(q, q_dot, f_ext);
        //do an inward recursion for inertia, forces and constraints
        this->downwards_sweep(alfa, torques);
        //Solve for the constraint forces
		this->constraint_calculation(beta);
        //do an upward recursion to propagate the result
        this->final_upwards_sweep(q_dotdot, torques);
    }

    void ChainIdSolver_Constraint_Vereshchagin::initial_upwards_sweep(const JntArray &q, const JntArray &qdot, const Wrenches& f_ext)
	{
        unsigned int j=0;
        F_total=Frame::Identity();
        for(unsigned int i=0;i<ns;i++){
            //Express everything in the segments reference frame (body coordinates)
            //which is at the segments tip, i.e. where the next joint is attached.
            
            //Calculate segment properties: X,S,vj,cj
            const Segment& segment=chain.getSegment(i);
            segment_info& s=results[i+1];
            //The pose between the joint root and the segment tip (tip expressed in joint root coordinates)
            s.F=segment.pose(q(j));
            F_total=F_total*s.F;
            //The velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
            Twist vj=s.F.M.Inverse(segment.twist(q(j),qdot(j)));
            //The unit velocity due to the joint motion of the segment expressed in the segments reference frame (tip)
            s.Z=s.F.M.Inverse(segment.twist(q(j),1.0));
            //Put Z in the joint root reference frame:
            s.Z=s.F*s.Z;
            //std::cout<<"q: "<<q(j)<<std::endl;
            //The total velocity of the segment expressed in the the segments reference frame (tip)
            if(i!=1)
                s.v=s.F.Inverse(results[i].v)+vj;
            else
                s.v=vj;
            
            //c[i] = cj + v[i]xvj (remark: cj=0, since our S is not time dependent in local coordinates)
            //The velocity product acceleration
            s.C = s.v*vj;//This is a cross product
            //Put C in the joint root reference frame
            s.C=s.F*s.C;
            //The rigid body inertia of the segment, expressed in the segments reference frame (tip)
            s.H=segment.getInertia();
        
            //wrench of the rigid body bias forces and the external forces on the segment (in body coordinates, tip)
            s.U = s.v*(s.H*s.v) - f_ext[i];
            
            if(segment.getJoint().getType()!=Joint::None)
                {
                    j++;
                }
        }
    }

	void ChainIdSolver_Constraint_Vereshchagin::downwards_sweep(const Jacobian& alfa,const JntArray &torques)
	{
        for(int i=ns;i>=0;i--)
            {
                //Get a handle for the segment we are working on.
                segment_info& s=results[i];
                //For segment N,
                //tilde is in the segment refframe (tip, not joint root)
                //without tilde is at the joint root (the childs tip!!!)
                //P_tilde is the articulated body inertia
                //R_tilde is the sum of external and coriolis/centrifugal forces
                //M is the (unit) acceleration energy already generated at link i
                //G is the (unit) magnitude of the constraint forces at link i
                //E are the (unit) constraint forces due to the constraints
                if(i==(ns)){
                    s.P_tilde=s.H;
                    s.R_tilde=s.U;
                    s.M.setZero();
                    s.G.setZero();
                    //changeBase(alfa_N,F_total.M.Inverse(),alfa_N2);
                    for(unsigned int r=0;r<3;r++)
                        for(unsigned int c=0;c<nc;c++){
                            s.E_tilde(r,c)=alfa(r+3,c);
                            s.E_tilde(r+3,c)=alfa(r,c);
                        }
                    //Change the reference frame of alfa to the segmentN tip frame
                    Rotation base_to_end=F_total.M.Inverse();
                    for(unsigned int c=0;c<nc;c++){
                        Wrench col(Vector(s.E_tilde(3,c),s.E_tilde(4,c),s.E_tilde(5,c)),
                                   Vector(s.E_tilde(0,c),s.E_tilde(1,c),s.E_tilde(2,c)));
                        col=base_to_end*col;
                        s.E_tilde.col(c)<<Vector3d::Map(col.torque.data),Vector3d::Map(col.force.data);
                    }
                }else{
                    //For all others:
                    //Everything should expressed in the body coordinates of segment i
                    segment_info& child=results[i+1];
                    //Copy PZ into a vector so we can do matrix manipulations, put torques above forces
                    Vector6d vPZ;
                    vPZ<<Vector3d::Map(child.PZ.torque.data),Vector3d::Map(child.PZ.force.data);
                    Matrix6d PZDPZt=(vPZ*vPZ.transpose()).lazy();
                    PZDPZt/=child.D;
                    //equation a) (see Vereshchagin89) PZDPZt=[I,H;H',M]
                    s.P_tilde=s.H+child.P-ArticulatedBodyInertia(PZDPZt.corner<3,3>(BottomRight),PZDPZt.corner<3,3>(TopRight),PZDPZt.corner<3,3>(TopLeft));
                    //equation b) (see Vereshchagin89)
                    s.R_tilde=s.U+child.R+child.PC+child.PZ/child.D*child.u;
                    //equation c) (see Vereshchagin89)
                    s.E_tilde=child.E;
                    s.E_tilde-=(vPZ*child.EZ.transpose()).lazy()/child.D;

                    //equation d) (see Vereshchagin89)
                    s.M=child.M;
                    s.M-=(child.EZ*child.EZ.transpose()).lazy()/child.D;
                    
                    //equation e) (see Vereshchagin89)
                    s.G=child.G;
                    Twist CiZDu=child.C+(child.Z/child.D)*child.u;
                    Vector6d vCiZDu;
                    vCiZDu<<Vector3d::Map(CiZDu.rot.data),Vector3d::Map(CiZDu.vel.data);
                    s.G+=(child.E.transpose()*vCiZDu).lazy();

                }
                if(i!=0){
                    //Transform all results to joint root coordinates of segment i (== body coordinates segment i-1)
                    //equation a)
                    s.P=s.F*s.P_tilde;
                    //equation b)
                    s.R=s.F*s.R_tilde;
                    //equation c), in matrix: torques above forces, so switch and switch back
                    for(unsigned int c=0;c<nc;c++){
                        Wrench col(Vector(s.E_tilde(3,c),s.E_tilde(4,c),s.E_tilde(5,c)),
                                   Vector(s.E_tilde(0,c),s.E_tilde(1,c),s.E_tilde(2,c)));
                        col=s.F*col;
                        s.E.col(c)<<Vector3d::Map(col.torque.data),Vector3d::Map(col.force.data);
                    }
                    
                    //needed for next recursion
                    s.PZ=s.P*s.Z;
                    s.D=dot(s.Z,s.PZ);
                    s.PC=s.P*s.C;
                    
                    //u=(Q-Z(R+PC)=sum of external forces along the joint axes, 
                    //R are the forces comming from the children, 
                    //Q is taken zero (do we need to take the previous calculated torques?
                    s.u=torques(i-1)-dot(s.Z,s.R+s.PC);
                    
                    //Matrix form of Z, put rotations above translations
                    Vector6d vZ;
                    vZ<<Vector3d::Map(s.Z.rot.data),Vector3d::Map(s.Z.vel.data);
                    s.EZ=(s.E.transpose()*vZ).lazy();
                }

                /*        
                std::cout<<"For segment "<<i<<std::endl;
                std::cout<<"D: "<<s.D<<std::endl;
                std::cout<<"Z: "<<s.Z<<std::endl;
                std::cout<<"PZ: "<<s.PZ<<std::endl;
                std::cout<<"E'Z: "<<s.EZ<<std::endl;
                std::cout<<"E~: "<<s.E_tilde<<std::endl;
                std::cout<<"E: "<<s.E<<std::endl;

                std::cout<<"E: "<<s.E<<std::endl;
                std::cout<<"Z: "<<s.Z.rot<<s.Z.vel<<std::endl;
                std::cout<<"G: "<<s.G<<std::endl;
                std::cout<<"M: "<<s.M<<std::endl;
                Matrix6d tmp;
                tmp<<s.P_tilde.I,s.P_tilde.H,s.P_tilde.H.transpose(),s.P_tilde.M;
                std::cout<<"P~: \n"<<tmp<<std::endl;
                tmp<<s.P.I,s.P.H,s.P.H.transpose(),s.P.M;
                std::cout<<"P: \n"<<tmp<<std::endl;
                */
          
            }
    }
    
    
	void ChainIdSolver_Constraint_Vereshchagin::constraint_calculation(const JntArray& beta)
	{
		//equation f) nu = M_0_inverse*(beta_N - E0_tilde`*acc0 - G0)
        //M_0_inverse, always nc*nc symmetric matrix
        //std::cout<<"M0: "<<results[0].M<<std::endl;
        //results[0].M-=MatrixXd::Identity(nc,nc);
        //std::cout<<"augmented M0: "<<results[0].M<<std::endl;
        
        //M_0_inverse=results[0].M.inverse();
        svd_eigen_HH(results[0].M,Um,Sm,Vm,tmpm);
        /*
        std::cout<<"U: "<<Um<<std::endl;
        std::cout<<"S: "<<Sm<<std::endl;
        std::cout<<"V: "<<Vm<<std::endl;
        */  
        //truncated svd, what would sdls, dls physically mean?
        for(unsigned int i=0;i<nc;i++)
            if(Sm(i)<1e-4)
                Sm(i)=0.0;
            else 
                Sm(i)=1/Sm(i);
        
        results[0].M=(Vm*Sm.asDiagonal()).lazy();
        M_0_inverse=(results[0].M*Um.transpose()).lazy();
        //results[0].M.ldlt().solve(MatrixXd::Identity(nc,nc),&M_0_inverse);
        //results[0].M.computeInverse(&M_0_inverse);
        //std::cout<<"inv(M0): "<<M_0_inverse<<std::endl;
        
        Vector6d acc;
        acc<<Vector3d::Map(acc_root.rot.data),Vector3d::Map(acc_root.vel.data);
        nu_sum=-(results[0].E_tilde.transpose()*acc).lazy();
        nu_sum+=beta.data;
        nu_sum -= results[0].G;
        nu = (M_0_inverse * nu_sum).lazy();
        /*
        std::cout<<"nu: "<<nu<<std::endl;
        std::cout<<"beta: "<<beta.data<<std::endl;
        std::cout<<"G: "<<results[0].G<<std::endl;
        */
    }    
    

	void ChainIdSolver_Constraint_Vereshchagin::final_upwards_sweep(JntArray &q_dotdot, JntArray &torques){
        {
            unsigned int j=0;
            
            for(unsigned int i=1;i<=ns;i++){
                segment_info& s=results[i];
                //Calculation of joint and segment accelerations
                //equation g) qdotdot[i] = D^-1*(Q - Z'(R + P(C + acc[i-1]) + E*nu))
                // = D^-1(u - Z'(P*acc[i-1] + E*nu)
                Twist a_p;
                if(i==1)
                    a_p=-acc_root;
                else
                    a_p=results[i-1].acc;
                
                //std::cout<<"a': "<<at<<std::endl;

                //The contribution of the constraint forces at segment i
                Vector6d tmp = s.E*nu;
                Wrench constraint_force = Wrench(Vector(tmp(3),tmp(4),tmp(5)),
                                                 Vector(tmp(0),tmp(1),tmp(2)));
                //std::cout<<"constraint force"<<i<<": "<<constraint_force<<std::endl;
                //Contribution of the acceleration of the parent (i-1)
                Wrench parent_force = s.P*a_p;
                //std::cout<<"parent force"<<i<<": "<<parent_force<<std::endl;
                //The constraint force and acceleration force projected on the joint axes -> axis torque/force
                double constraint_torque = dot(s.Z,parent_force+constraint_force);
                //The result should be the torque at this joint
                torques(j) = (s.u-constraint_torque);
                q_dotdot(j) = torques(j)/s.D;
                
                //equation h) acc[i] = Fi*(acc[i-1] + Z[i]*qdotdot[i] + c[i]
                s.acc=s.F.Inverse(a_p+s.Z*q_dotdot(j)+s.C);
                
                //std::cout<<"a: "<<s.acc<<std::endl;

                if(chain.getSegment(i).getJoint().getType()!=Joint::None)
                    {
                        j++;
                    }
            }
        }
    }
}//namespace
