#include "solvertest.hpp"
#include <frames_io.hpp>
#include <framevel_io.hpp>
#include <kinfam_io.hpp>
#include <time.h>
#include <fstream>
#include <chainidsolver_constraint_vereshchagin.hpp>

CPPUNIT_TEST_SUITE_REGISTRATION( SolverTest );

using namespace KDL;

void SolverTest::setUp()
{
    srand( (unsigned)time( NULL ));

    chain1.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotZ),
                             Frame(Vector(0.0,0.0,0.0))));
    chain1.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotX),
                             Frame(Vector(0.0,0.0,0.9))));
    chain1.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::None),
                             Frame(Vector(-0.4,0.0,0.0))));
    chain1.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotX),
                             Frame(Vector(0.0,0.0,1.2))));
    chain1.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::None),
                             Frame(Vector(0.4,0.0,0.0))));
    chain1.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),
                             Frame(Vector(0.0,0.0,1.4))));
    chain1.addSegment(Segment("Segment 7", Joint("Joint 7", Joint::RotX),
                             Frame(Vector(0.0,0.0,0.0))));
    chain1.addSegment(Segment("Segment 8", Joint("Joint 8", Joint::RotZ),
                             Frame(Vector(0.0,0.0,0.4))));
    chain1.addSegment(Segment("Segment 9", Joint("Joint 9", Joint::None),
                             Frame(Vector(0.0,0.0,0.0))));

    chain2.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotZ),
                              Frame(Vector(0.0,0.0,0.5))));
    chain2.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotX),
                              Frame(Vector(0.0,0.0,0.4))));
    chain2.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),
                              Frame(Vector(0.0,0.0,0.3))));
    chain2.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotX),
                              Frame(Vector(0.0,0.0,0.2))));
    chain2.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotZ),
                              Frame(Vector(0.0,0.0,0.1))));


    chain3.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotZ),
                             Frame(Vector(0.0,0.0,0.0))));
    chain3.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotX),
                             Frame(Vector(0.0,0.0,0.9))));
    chain3.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotZ),
                             Frame(Vector(-0.4,0.0,0.0))));
    chain3.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotX),
                             Frame(Vector(0.0,0.0,1.2))));
    chain3.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::None),
                             Frame(Vector(0.4,0.0,0.0))));
    chain3.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),
                             Frame(Vector(0.0,0.0,1.4))));
    chain3.addSegment(Segment("Segment 7", Joint("Joint 7", Joint::RotX),
                             Frame(Vector(0.0,0.0,0.0))));
    chain3.addSegment(Segment("Segment 8", Joint("Joint 8", Joint::RotZ),
                             Frame(Vector(0.0,0.0,0.4))));
    chain3.addSegment(Segment("Segment 9", Joint("Joint 9", Joint::RotY),
                             Frame(Vector(0.0,0.0,0.0))));


    chain4.addSegment(Segment("Segment 1", Joint("Joint 1", Vector(10,0,0), Vector(1,0,1),Joint::RotAxis),
			       Frame(Vector(0.0,0.0,0.5))));
    chain4.addSegment(Segment("Segment 2", Joint("Joint 2", Vector(0,5,0), Vector(1,0,0),Joint::RotAxis),
			       Frame(Vector(0.0,0.0,0.4))));
    chain4.addSegment(Segment("Segment 3", Joint("Joint 3", Vector(0,0,5), Vector(1,0,4),Joint::RotAxis),
                              Frame(Vector(0.0,0.0,0.3))));
    chain4.addSegment(Segment("Segment 4", Joint("Joint 4", Vector(0,0,0), Vector(1,0,0),Joint::RotAxis),
                              Frame(Vector(0.0,0.0,0.2))));
    chain4.addSegment(Segment("Segment 5", Joint("Joint 5", Vector(0,0,0), Vector(0,0,1),Joint::RotAxis),
                              Frame(Vector(0.0,0.0,0.1))));


   //chain definition for vereshchagin's dynamic solver
    Joint rotJoint0 = Joint(Joint::RotZ, 1, 0, 0.01);
    Joint rotJoint1 = Joint(Joint::RotZ, 1, 0, 0.01);

    Frame refFrame(Rotation::RPY(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0));
    Frame frame1(Rotation::RPY(0.0, 0.0, 0.0), Vector(0.0, -0.4, 0.0));
    Frame frame2(Rotation::RPY(0.0, 0.0, 0.0), Vector(0.0, -0.4, 0.0));

    Segment segment1 = Segment(rotJoint0, frame1);
    Segment segment2 = Segment(rotJoint1, frame2);

    //rotational inertia around symmetry axis of rotation
    RotationalInertia rotInerSeg1(0.0, 0.0, 0.0, 0.0, 0.0, 0.0); 

    //spatial inertia
    RigidBodyInertia inerSegment1(0.3, Vector(0.0, -0.4, 0.0), rotInerSeg1);
    RigidBodyInertia inerSegment2(0.3, Vector(0.0, -0.4, 0.0), rotInerSeg1);
    segment1.setInertia(inerSegment1);
    segment2.setInertia(inerSegment2);

    //chain segments
    chaindyn.addSegment(segment1);
    chaindyn.addSegment(segment2);
    
}

void SolverTest::tearDown()
{
//     delete fksolverpos;
//     delete jacsolver;
//     delete fksolvervel;
//     delete iksolvervel;
//     delete iksolverpos;
}

void SolverTest::FkPosAndJacTest()
{
    ChainFkSolverPos_recursive fksolver1(chain1);
    ChainJntToJacSolver jacsolver1(chain1);
    FkPosAndJacLocal(chain1,fksolver1,jacsolver1);
    ChainFkSolverPos_recursive fksolver2(chain2);
    ChainJntToJacSolver jacsolver2(chain2);
    FkPosAndJacLocal(chain2,fksolver2,jacsolver2);
    ChainFkSolverPos_recursive fksolver3(chain3);
    ChainJntToJacSolver jacsolver3(chain3);
    FkPosAndJacLocal(chain3,fksolver3,jacsolver3);
    ChainFkSolverPos_recursive fksolver4(chain4);
    ChainJntToJacSolver jacsolver4(chain4);
    FkPosAndJacLocal(chain4,fksolver4,jacsolver4);
}

void SolverTest::FkVelAndJacTest()
{
    ChainFkSolverVel_recursive fksolver1(chain1);
    ChainJntToJacSolver jacsolver1(chain1);
    FkVelAndJacLocal(chain1,fksolver1,jacsolver1);
    ChainFkSolverVel_recursive fksolver2(chain2);
    ChainJntToJacSolver jacsolver2(chain2);
    FkVelAndJacLocal(chain2,fksolver2,jacsolver2);
    ChainFkSolverVel_recursive fksolver3(chain3);
    ChainJntToJacSolver jacsolver3(chain3);
    FkVelAndJacLocal(chain3,fksolver3,jacsolver3);
    ChainFkSolverVel_recursive fksolver4(chain4);
    ChainJntToJacSolver jacsolver4(chain4);
    FkVelAndJacLocal(chain4,fksolver4,jacsolver4);
}

void SolverTest::FkVelAndIkVelTest()
{
    //Chain1
    std::cout<<"square problem"<<std::endl;
    ChainFkSolverVel_recursive fksolver1(chain1);
    ChainIkSolverVel_pinv iksolver1(chain1);
    ChainIkSolverVel_pinv_givens iksolver_pinv_givens1(chain1);
    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkVelAndIkVelLocal(chain1,fksolver1,iksolver1);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkVelAndIkVelLocal(chain1,fksolver1,iksolver_pinv_givens1);

    //Chain2
    std::cout<<"underdetermined problem"<<std::endl;
    ChainFkSolverVel_recursive fksolver2(chain2);
    ChainIkSolverVel_pinv iksolver2(chain2);
    ChainIkSolverVel_pinv_givens iksolver_pinv_givens2(chain2);
    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkVelAndIkVelLocal(chain2,fksolver2,iksolver2);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkVelAndIkVelLocal(chain2,fksolver2,iksolver_pinv_givens2);

    //Chain3
    std::cout<<"overdetermined problem"<<std::endl;
    ChainFkSolverVel_recursive fksolver3(chain3);
    ChainIkSolverVel_pinv iksolver3(chain3);
    ChainIkSolverVel_pinv_givens iksolver_pinv_givens3(chain3);
    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkVelAndIkVelLocal(chain3,fksolver3,iksolver3);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkVelAndIkVelLocal(chain3,fksolver3,iksolver_pinv_givens3);

    //Chain4
    std::cout<<"overdetermined problem"<<std::endl;
    ChainFkSolverVel_recursive fksolver4(chain4);
    ChainIkSolverVel_pinv iksolver4(chain4);
    ChainIkSolverVel_pinv_givens iksolver_pinv_givens4(chain4);
    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkVelAndIkVelLocal(chain4,fksolver4,iksolver4);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkVelAndIkVelLocal(chain4,fksolver4,iksolver_pinv_givens4);
}

void SolverTest::FkPosAndIkPosTest()
{
    std::cout<<"square problem"<<std::endl;
    ChainFkSolverPos_recursive fksolver1(chain1);
    ChainIkSolverVel_pinv iksolver1v(chain1);
    ChainIkSolverVel_pinv_givens iksolverv_pinv_givens1(chain1);
    ChainIkSolverPos_NR iksolver1(chain1,fksolver1,iksolver1v);
    ChainIkSolverPos_NR iksolver1_givens(chain1,fksolver1,iksolverv_pinv_givens1,1000);

    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkPosAndIkPosLocal(chain1,fksolver1,iksolver1);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkPosAndIkPosLocal(chain1,fksolver1,iksolver1_givens);

    std::cout<<"underdetermined problem"<<std::endl;
    ChainFkSolverPos_recursive fksolver2(chain2);
    ChainIkSolverVel_pinv iksolverv2(chain2);
    ChainIkSolverVel_pinv_givens iksolverv_pinv_givens2(chain2);
    ChainIkSolverPos_NR iksolver2(chain2,fksolver2,iksolverv2);
    ChainIkSolverPos_NR iksolver2_givens(chain2,fksolver2,iksolverv_pinv_givens2,1000);

    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkPosAndIkPosLocal(chain2,fksolver2,iksolver2);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkPosAndIkPosLocal(chain2,fksolver2,iksolver2_givens);

    std::cout<<"overdetermined problem"<<std::endl;
    ChainFkSolverPos_recursive fksolver3(chain3);
    ChainIkSolverVel_pinv iksolverv3(chain3);
    ChainIkSolverVel_pinv_givens iksolverv_pinv_givens3(chain3);
    ChainIkSolverPos_NR iksolver3(chain3,fksolver3,iksolverv3);
    ChainIkSolverPos_NR iksolver3_givens(chain3,fksolver3,iksolverv_pinv_givens3,1000);

    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkPosAndIkPosLocal(chain3,fksolver3,iksolver3);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkPosAndIkPosLocal(chain3,fksolver3,iksolver3_givens);

    std::cout<<"underdetermined problem with WGs segment constructor"<<std::endl;
    ChainFkSolverPos_recursive fksolver4(chain4);
    ChainIkSolverVel_pinv iksolverv4(chain4);
    ChainIkSolverVel_pinv_givens iksolverv_pinv_givens4(chain4);
    ChainIkSolverPos_NR iksolver4(chain4,fksolver4,iksolverv4,1000);
    ChainIkSolverPos_NR iksolver4_givens(chain4,fksolver4,iksolverv_pinv_givens4,1000);

    std::cout<<"KDL-SVD-HouseHolder"<<std::endl;
    FkPosAndIkPosLocal(chain4,fksolver4,iksolver4);
    std::cout<<"KDL-SVD-Givens"<<std::endl;
    FkPosAndIkPosLocal(chain4,fksolver4,iksolver4_givens);
}

void SolverTest::VereshchaginTest(){

    Vector constrainXLinear(0.0, 0.0, 0.0);
    Vector constrainXAngular(0.0, 0.0, 0.0);
    Vector constrainYLinear(0.0, 1.0, 0.0);
    Vector constrainYAngular(0.0, 0.0, 0.0);
    // Vector constrainZLinear(0.0, 0.0, 0.0);
    //Vector constrainZAngular(0.0, 0.0, 0.0);
    Twist constraintForcesX(constrainXLinear, constrainXAngular);
    Twist constraintForcesY(constrainYLinear, constrainYAngular);
    //Twist constraintForcesZ(constrainZLinear, constrainZAngular);
    Jacobian alpha(1);
    //alpha.setColumn(0, constraintForcesX);
    alpha.setColumn(0, constraintForcesY);
    //alpha.setColumn(0, constraintForcesZ);

    //Acceleration energy at  the end-effector
    JntArray betha(1); //set to zero
    betha(0) = 0.0;
    //betha(1) = 0.0;
    //betha(2) = 0.0;

    //arm root acceleration
    Vector linearAcc(0.0, 10, 0.0); //gravitational acceleration along Y
    Vector angularAcc(0.0, 0.0, 0.0);
    Twist twist1(linearAcc, angularAcc);

    //external forces on the arm
    Vector externalForce1(0.0, 0.0, 0.0);
    Vector externalTorque1(0.0, 0.0, 0.0);
    Vector externalForce2(0.0, 0.0, 0.0);
    Vector externalTorque2(0.0, 0.0, 0.0);
    Wrench externalNetForce1(externalForce1, externalTorque1);
    Wrench externalNetForce2(externalForce2, externalTorque2);
    Wrenches externalNetForce;
    externalNetForce.push_back(externalNetForce1);
    externalNetForce.push_back(externalNetForce2);
    //~Definition of constraints and external disturbances
    //-------------------------------------------------------------------------------------//


    //Definition of solver and initial configuration
    //-------------------------------------------------------------------------------------//
    int numberOfConstraints = 1;
    ChainIdSolver_Vereshchagin constraintSolver(chain, twist1, numberOfConstraints);
    
    //These arrays of joint values contain actual and desired values
    //actual is generated by a solver and integrator
    //desired is given by an interpolator
    //error is the difference between desired-actual
    int k = 3;
    JntArray jointPoses[k];
    JntArray jointRates[k];
    JntArray jointAccelerations[k];
    JntArray jointTorques[k];
    JntArray biasqDotDot(chain.getNrOfJoints());
    for (int i = 0; i < k; i++)
    {
        JntArray jointValues(chain.getNrOfJoints());
        jointPoses[i] = jointValues;
        jointRates[i] = jointValues;
        jointAccelerations[i] = jointValues;
        jointTorques[i] = jointValues;
    }

    //cartesian space/link values
    k = 4;
    Frames cartX[k];
    Twists cartXDot[k];
    Twists cartXDotDot[k];
    Twist accLink;
    for (int i = 0; i < k; i++) //i is number of variables (actual, desired, error)
    {
        for (int j = 0; j < k-1; j++) //j is number of links
        {
            cartX[i].push_back(frame1);
            cartXDot[i].push_back(accLink);
            cartXDotDot[i].push_back(accLink);
        }

    }


    // Initial arm position configuration/constraint
    JntArray jointInitialPose(chain.getNrOfJoints());
    jointInitialPose(0) = 0.0; // initial joint0 pose
    jointInitialPose(1) = M_PI/6.0; //initial joint1 pose, negative in clockwise
    //j0=0.0, j1=pi/6.0 correspond to x = 0.2, y = -0.7464
    //j0=2*pi/3.0, j1=pi/4.0 correspond to x = 0.44992, y = 0.58636

    //actual
    jointPoses[0](0) = jointInitialPose(0);
    jointPoses[0](1) = jointInitialPose(1);
    //desired
    jointPoses[1](0) = jointInitialPose(0);
    jointPoses[1](1) = jointInitialPose(1);

   
    //~Definition of solver and initial configuration
    //-------------------------------------------------------------------------------------//



    //Definition of process main loop
    //-------------------------------------------------------------------------------------//
    double taskTimeConstant = 2; //Time required to complete the task move(frameinitialPose, framefinalPose) default T=10.0
    double simulationTime = 2.5*taskTimeConstant;
    double timeDelta = 0.001;
    bool status;
    
    double ksi = 1; //damping factor
    double Kp = 0.02/(taskTimeConstant*taskTimeConstant);
    double Kv = 1*ksi/taskTimeConstant;
    double Ki = 1.0;
    double Ka = 0.0;


    //Interpolator parameters:
    double b0_y = -0.7464102; //should come from initial joint configuration
    double b1_y = 0.0;
    double b2_y = 0.0;//((0.5 + 0.7464)*3.0 / TimeConstant * TimeConstant);
    double b3_y = 0.0;//-((0.5 + 0.7464)*2.0 / TimeConstant * TimeConstant * TimeConstant);

    double b0_x = 0.2; //should come from initial joint configuration
    double b1_x = 0.0;
    double b2_x = 0.0;//((0.5 + 0.7464)*3.0 / TimeConstant * TimeConstant);
    double b3_x = 0.0;//-((0.5 + 0.7464)*2.0 / TimeConstant * TimeConstant * TimeConstant);

    /*
    double a0_q0 = jointPoses[1](0);
    double a1_q0 = 0.0;
    double a2_q0 = ((qFinalPose(0) - jointPoses[1](0))*3.0 / TimeConstant * TimeConstant);
    double a3_q0 = -((qFinalPose(0) - jointPoses[1](0))*2.0 / TimeConstant * TimeConstant * TimeConstant);

    double a0_q1 = jointPoses[1](1);
    double a1_q1 = 0.0;
    double a2_q1 = ((qFinalPose(1) - jointPoses[1](1))*3.0 / TimeConstant * TimeConstant);
    double a3_q1 = -((qFinalPose(1) - jointPoses[1](1))*2.0 / TimeConstant * TimeConstant * TimeConstant);
     */

    
    for (double t = 0.0; t <=simulationTime; t = t + timeDelta)
    {

        //Interpolation (Desired) q = a0+a1t+a2t^2+a3t^3
        //Do we feed these values? then         jointPoses[0] = jointPoses[1];         jointRates[0] = jointRates[1];
        // But I don't think so, in control desired trajectory plays role of the reference that the controller should strive for from its actual (previous, current) state.

        cartX[1][1].p[0] = b0_x + b1_x * t + b2_x * t * t + b3_x * t * t*t;
        cartX[1][1].p[1] = b0_y + b1_y * t + b2_y * t * t + b3_y * t * t*t;
        cartXDot[1][1].vel[0] = b1_x + 2 * b2_x * t + 3 * b3_x * t*t;
        cartXDot[1][1].vel[1] = b1_y + 2 * b2_y * t + 3 * b3_y * t*t;
        cartXDotDot[1][1].vel[0] = 2 * b2_x + 6 * b3_x*t;
        cartXDotDot[1][1].vel[1] = 2 * b2_y + 6 * b3_y*t;
       // printf("Desired Cartesian values: %f          %f      %f     %f         %f        %f      %f\n", t, cartX[1][1].p[0], cartX[1][1].p[1], cartXDot[1][1].vel[0], cartXDot[1][1].vel[1], cartXDotDot[1][1].vel[0], cartXDotDot[1][1].vel[1]);


        status = constraintSolver.CartToJnt(jointPoses[0], jointRates[0], jointAccelerations[0], alpha, betha, externalNetForce, jointTorques[0]);
        /*
        constraintSolver.initial_upwards_sweep(jointPoses[0], jointRates[0], jointAccelerations[0], externalNetForce);
        constraintSolver.downwards_sweep(alpha, jointTorques[0]);
        constraintSolver.constraint_calculation(betha);
        constraintSolver.final_upwards_sweep(jointAccelerations[0], cartXDotDot[0], jointTorques[0]);
       */
        //Integration(robot joint values for rates and poses; actual) at the given "instanteneous" interval for joint position and velocity.
        jointRates[0](0) = jointRates[0](0) + jointAccelerations[0](0) * timeDelta; //Euler Forward
        jointPoses[0](0) = jointPoses[0](0) + (jointRates[0](0) - jointAccelerations[0](0) * timeDelta / 2.0) * timeDelta; //Trapezoidal rule
        jointRates[0](1) = jointRates[0](1) + jointAccelerations[0](1) * timeDelta; //Euler Forward
        jointPoses[0](1) = jointPoses[0](1) + (jointRates[0](1) - jointAccelerations[0](1) * timeDelta / 2.0) * timeDelta;
        //printf("Joint %f          %f      %f       %f     %f       %f      %f     %f      %f\n", t, jointPoses[0](0), jointPoses[0](1), jointRates[0](0), jointRates[0](1), jointAccelerations[0](0), jointAccelerations[0](1), jointTorques[0](0), jointTorques[0](1));
        
        //constraintSolver.initial_upwards_sweep(jointPoses[0], jointRates[0], jointAccelerations[0], externalNetForce);
        //constraintSolver.getLinkCartesianPose(cartX[0]);
        //constraintSolver.getLinkCartesianVelocity(cartXDot[0]);
        //constraintSolver.getLinkCartesianAcceleration(cartXDotDot[0]);
        //printf("Cartesian actual %f          %f      %f     %f         %f        %f      %f     %f\n", t, cartX[0][1].p.x(), cartX[0][1].p.y(), cartXDot[0][1].vel[0], cartXDot[0][1].vel[1], cartXDotDot[0][1].vel[0], cartXDotDot[0][1].vel[1], cartXDotDot[0][1].rot[2]);

        //Error
        /*
        jointPoses[2](0) = jointPoses[1](0) - jointPoses[0](0);
        jointPoses[2](1) = jointPoses[1](1) - jointPoses[0](1);
        jointRates[2](0) = jointRates[1](0) - jointRates[0](0);
        jointRates[2](1) = jointRates[1](1) - jointRates[0](1);
        jointAccelerations[2](0) = jointAccelerations[1](0) - jointAccelerations[0](0);
        jointAccelerations[2](1) = jointAccelerations[1](1) - jointAccelerations[0](1);
        printf("Errors: %f          %f      %f          %f     %f       %f      %f\n", t, jointPoses[2](0), jointPoses[2](1), jointRates[2](0), jointRates[2](1), jointAccelerations[2](0), jointAccelerations[2](1));
        */
        //e = d - a
        //cartX[2][1].p[0] = cartX[1][1].p[0] - cartX[0][1].p[0];
        //cartX[2][1].p[1] = cartX[1][1].p[1] - cartX[0][1].p[1];
        //cartXDot[2][1].vel[0] = cartXDot[1][1].vel[0] - cartXDot[0][1].vel[0];
        //cartXDot[2][1].vel[1] = cartXDot[1][1].vel[1] - cartXDot[0][1].vel[1];
        //cartXDotDot[2][1].vel[0] = cartXDotDot[1][1].vel[0] - cartXDotDot[0][1].vel[0];
        //cartXDotDot[2][1].vel[1] = cartXDotDot[1][1].vel[1] - cartXDotDot[0][1].vel[1];
        //cartXDotDot[2][1].rot[2] = cartXDotDot[1][1].rot[2] - cartXDotDot[0][1].rot[2];
        //cartX[3][1].p[0] += timeDelta*cartX[2][1].p[0]; //for integral term;
        //cartX[3][1].p[1] += timeDelta*cartX[2][1].p[1];
        //printf("Cartesian error  %f          %f      %f     %f         %f        %f      %f     %f\n", t, cartX[2][1].p[0], cartX[2][1].p[1], cartXDot[2][1].vel[0], cartXDot[2][1].vel[1], cartXDotDot[2][1].vel[1], cartXDotDot[2][1].vel[1], cartXDotDot[2][1].rot[2]);
        //Regulator or controller
        //externalNetForce[1].force[0] =  ( (Kv)*cartXDot[2][1].vel[0] + (Kp)*cartX[2][1].p[0] + Ki*cartX[3][1].p[0]);

        
        //betha(0) = ( (Kv)*cartXDot[2][1].vel[1] - (Ka)*cartXDot[2][1].rot[2] + (Kp)*cartX[2][1].p[1]);// + Ki*cartX[3][1].p[0]);
        //betha(1) = betha(1) +(  (Kv)*cartXDot[2][1].vel[1] + (Kp)*cartX[2][1].p[1]);// + Ki*cartX[3][1].p[0]);
        //printf("betha   %f      %f\n", betha(0), betha(1));
    


/*
    unsigned int nc=1;
    Jacobian alfa(nc);
    JntArray beta(nc);
    alfa.data.transpose()<<0,0,1,0,0,0;
    //0,0,0,0,0,0;
//        0,1,0,0,0,0,
//        0,0,0,1,0,0,
//        0,0,0,0,1,0,
//        0,0,0,0,0,1;
        
    ChainIdSolver_Constraint_Vereshchagin ikacc(chaindyn,Twist(Vector(0.0,0,-9.81),Vector::Zero()),nc);
    ChainFkSolverVel_recursive fkposvel(chaindyn);
    JntArrayVel q(chaindyn.getNrOfJoints());
    JntArray qdotdot(chaindyn.getNrOfJoints());
    JntArray torques(chaindyn.getNrOfJoints());
    Wrenches f_ext(chaindyn.getNrOfSegments());
    SetToZero(q);
    SetToZero(torques);
    for(unsigned int i=0;i<chaindyn.getNrOfSegments();i++)
        SetToZero(f_ext[i]);
    q.q(0)=0.0;
    q.q(1)=M_PI_4;
    //q.q(0)=0;
    //q.q(1)=0.5;
    FrameVel ee_start,ee;
    std::cout<<"q: "<<q.q<<std::endl;
    fkposvel.JntToCart(q,ee_start);
    fkposvel.JntToCart(q,ee);
    std::cout<<"ee: "<<ee.p.p.y()<<","<<ee.p.p.z()<<std::endl;
    double T=2.;
    double dt=0.001;
    std::ofstream fout;
    fout.open("vereshchagin_result.txt");
    fout<<"time"<<"  "<<"torque_1"<<"  "<<"torque_2"<<"  "<<"q_1"<<"  "<<"q_2"<<"  "<<"qdot_1"<<"  "<<"qdot_2"<<"  "<<"qdotdot_1"<<"  "<<"qdotdot_2"<<"  "<<"y"<<"  "<<"z"<<std::endl;
    clock_t start, finish;
    start = clock();        
    int times=100;
    double prev_error=0.0;
    double integral=0.0;
    for(double t=0.0;t<10;t+=dt){
        //f_ext[1]=ee.M.R.Inverse(Wrench(Vector(0.0,-1.0,1.0 ),Vector::Zero()));
        //beta.data<<1;//,0;
        //double intbeta0=0.0;
        double error=0.5-ee.p.p.z();
        integral+=error*dt;
        double deriv = (error-prev_error)/dt;
        beta.data(0)=20.0*error+15.*integral+0.*deriv-8.0*ee.p.v.z();
        //beta.data(1)=2000*(ee_start.p.p.y()-ee.p.p.y())-20.0*ee.p.v.y();//0.0;//(ee_start.p.p.y()-ee.p.p.y())/(T*T);
        //std::cout<<"beta: "<<beta<<std::endl;
        torques(0)=-0.5*q.qdot(0);
        torques(1)=-0.7*q.qdot(1);
        ikacc.CartToJnt(q.q, q.qdot, qdotdot,alfa,beta,f_ext,torques);
        q.q.data+=(q.qdot.data+qdotdot.data*dt/2)*dt;
        q.qdot.data+=qdotdot.data*dt;
        /*
        std::cout<<"qdotdot:"<<qdotdot<<std::endl;
        std::cout<<"qdot:"<<q.qdot<<std::endl;
        std::cout<<"q: "<<q.q<<std::endl;
        */
        fout<<t<<"   ";
        fout<<torques(0)<<"  "<<torques(1)<<"  ";
        fout<<q.q(0)<<"  "<<q.q(1)<<"   ";
        fout<<q.qdot(0)<<"  "<<q.qdot(1)<<"   ";
        fout<<qdotdot(0)<<"  "<<qdotdot(1)<<"   ";
        
        fkposvel.JntToCart(q,ee);
        
        fout<<ee.p.p.y()<<"  "<<ee.p.p.z()<<std::endl;
        /*
        std::cout<<"ee.p: "<<ee.p.p.y()<<","<<ee.p.p.z()<<std::endl;
        std::cout<<"ee.v: "<<ee.p.v.y()<<","<<ee.p.v.z()<<std::endl;
        */
    }
    fout.close();
    finish = clock();
        
    //std::cout<<(double(finish - start))<<std::endl;
    
    std::cout<<(double(finish - start)/CLOCKS_PER_SEC )<<std::endl;
    
    //std::cout<<"q: "<<q.q<<std::endl;
    //std::cout<<"u: "<<torques<<std::endl;
    fkposvel.JntToCart(q,ee);
    std::cout<<"ee_start: "<<ee_start.p.p.y()<<","<<ee_start.p.p.z()<<std::endl;
    std::cout<<"ee.p: "<<ee.p.p.y()<<","<<ee.p.p.z()<<std::endl;
    //std::cout<<"ee.v: "<<ee.p.v.y()<<","<<ee.p.v.z()<<std::endl;
    
*/
}


void SolverTest::FkPosAndJacLocal(Chain& chain,ChainFkSolverPos& fksolverpos,ChainJntToJacSolver& jacsolver)
{
    double deltaq = 1E-4;
    double epsJ   = 1E-4;

    Frame F1,F2;

    JntArray q(chain.getNrOfJoints());
    Jacobian jac(chain.getNrOfJoints());

    for(unsigned int i=0;i<chain.getNrOfJoints();i++){
        random(q(i));
    }

    jacsolver.JntToJac(q,jac);

    for (int i=0; i< q.rows() ;i++) {
        // test the derivative of J towards qi
        double oldqi = q(i);
        q(i) = oldqi+deltaq;
        CPPUNIT_ASSERT(0==fksolverpos.JntToCart(q,F2));
        q(i) = oldqi-deltaq;
        CPPUNIT_ASSERT(0==fksolverpos.JntToCart(q,F1));
        q(i) = oldqi;
        // check Jacobian :
        Twist Jcol1 = diff(F1,F2,2*deltaq);
        Twist Jcol2(Vector(jac(0,i),jac(1,i),jac(2,i)),
                    Vector(jac(3,i),jac(4,i),jac(5,i)));

        //CPPUNIT_ASSERT_EQUAL(true,Equal(Jcol1,Jcol2,epsJ));
        CPPUNIT_ASSERT_EQUAL(Jcol1,Jcol2);
    }
}

void SolverTest::FkVelAndJacLocal(Chain& chain, ChainFkSolverVel& fksolvervel, ChainJntToJacSolver& jacsolver)
{
    double deltaq = 1E-4;
    double epsJ   = 1E-4;

    JntArray q(chain.getNrOfJoints());
    JntArray qdot(chain.getNrOfJoints());

    for(unsigned int i=0;i<chain.getNrOfJoints();i++){
        random(q(i));
        random(qdot(i));
    }
    JntArrayVel qvel(q,qdot);
    Jacobian jac(chain.getNrOfJoints());

    FrameVel cart;
    Twist t;

    jacsolver.JntToJac(qvel.q,jac);
    CPPUNIT_ASSERT(fksolvervel.JntToCart(qvel,cart)==0);
    MultiplyJacobian(jac,qvel.qdot,t);
    CPPUNIT_ASSERT_EQUAL(cart.deriv(),t);
}

void SolverTest::FkVelAndIkVelLocal(Chain& chain, ChainFkSolverVel& fksolvervel, ChainIkSolverVel& iksolvervel)
{
    double epsJ   = 1E-7;

    JntArray q(chain.getNrOfJoints());
    JntArray qdot(chain.getNrOfJoints());

    for(unsigned int i=0;i<chain.getNrOfJoints();i++){
        random(q(i));
        random(qdot(i));
    }
    JntArrayVel qvel(q,qdot);
    JntArray qdot_solved(chain.getNrOfJoints());
        
    FrameVel cart;
    
    CPPUNIT_ASSERT(0==fksolvervel.JntToCart(qvel,cart));
    
    int ret = iksolvervel.CartToJnt(qvel.q,cart.deriv(),qdot_solved);
    CPPUNIT_ASSERT(0<=ret);
    
    qvel.deriv()=qdot_solved;
    
    if(chain.getNrOfJoints()<=6)
        CPPUNIT_ASSERT(Equal(qvel.qdot,qdot_solved,1e-5));
    else{
        FrameVel cart_solved;
        CPPUNIT_ASSERT(0==fksolvervel.JntToCart(qvel,cart_solved));
        CPPUNIT_ASSERT(Equal(cart.deriv(),cart_solved.deriv(),1e-5));
    }
}


void SolverTest::FkPosAndIkPosLocal(Chain& chain,ChainFkSolverPos& fksolverpos, ChainIkSolverPos& iksolverpos)
{
    JntArray q(chain.getNrOfJoints());
    for(unsigned int i=0;i<chain.getNrOfJoints();i++){
        random(q(i));
    }
    JntArray q_init(chain.getNrOfJoints());
    double tmp;
    for(unsigned int i=0;i<chain.getNrOfJoints();i++){
        random(tmp);
        q_init(i)=q(i)+0.1*tmp;
    }
    JntArray q_solved(q);

    Frame F1,F2;

    CPPUNIT_ASSERT(0==fksolverpos.JntToCart(q,F1));
    CPPUNIT_ASSERT(0==iksolverpos.CartToJnt(q_init,F1,q_solved));
    CPPUNIT_ASSERT(0==fksolverpos.JntToCart(q_solved,F2));

    CPPUNIT_ASSERT_EQUAL(F1,F2);
    //CPPUNIT_ASSERT_EQUAL(q,q_solved);

}


