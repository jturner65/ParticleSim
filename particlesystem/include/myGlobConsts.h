#ifndef myGlobConsts_h
#define myGlobConsts_h
#include <map>
#include <list>
#include <vector>
#include <deque>
#include <Eigen/Dense>
#include <Eigen/StdVector>
using namespace std;
namespace particleSystem{
    
	static const double partRad = .04;

	static const double distFromSnowGlobe = 10;
	static const double shakingCoefficient = .1;

	static const double snowGlobStartRad = 1;			//internal snowglobe radius where snowflakes can start
	static const double snowGlobRad = 10;
	static const double partSqColDist = .1;				//potential particle collision distance
	static const double globDeltaT = .01;
	static const double epsVal = .000001;

	static const double dragTTFrcCoef = -.2;			//general damping for tinker toys
	static const double dragMSFrcCoef = -.3;			//general damping for mass spring
	static const double dragFluidFrcCoef = -.15;		//general damping for fluids

	static const double PI = atan(1.0) * 4;

	static const double twoPI = 2 * PI;
	static const double PI_2 = PI * .5;
	static const double VERLET1mDAMP = .99999;          //1 minus some tiny damping term for verlet stability
	static const int CONJGRAD_MAXITERS = 30;
	static const int  numPosCnstIters = 300;     //# of iterations for position constraint satisfaction


	/** Helper types for STL containers with fixed size vectorizable Eigen memory allocators - from MRPT. */
	template <class TYPE1, class TYPE2 = TYPE1>
	struct aligned_containers {
		typedef std::pair<TYPE1, TYPE2> pair_t;
		typedef std::vector<TYPE1, Eigen::aligned_allocator<TYPE1> > vector_t;
		typedef std::vector<std::vector<TYPE1, Eigen::aligned_allocator<TYPE1>>> vec_vec_t;
		typedef std::deque<TYPE1, Eigen::aligned_allocator<TYPE1> > deque_t;
		typedef std::list<TYPE1, Eigen::aligned_allocator<TYPE1> > list_t;
		typedef std::map<TYPE1, TYPE2, std::less<TYPE1>, Eigen::aligned_allocator<std::pair<const TYPE1, TYPE2> > > map_t;
		typedef std::multimap<TYPE1, TYPE2, std::less<TYPE1>, Eigen::aligned_allocator<std::pair<const TYPE1, TYPE2> > > multimap_t;
	};


    //10 when inv pend working      //number of scenes - 0:galileo, 1:snowglobe, 2:single constraint, 3:clacker, 4:4 circles, 5:4 circle jump-tracks, 6: 4 part inverted pendulum, 7:9 bead-on-wires with varying integrators, 8 : 9 bead on wires with varying lambda rk4 gen, 9 : seaweed field

	//particle - has current values for position, velocity and forceAcc in idx 0 of appropriate deques
	enum ForceType {F_NONE, S_SCALAR, S_VECTOR, ATTR, REPL, DAMPSPRING, DSPR_THETABAR};
	static const char* ForceType2str[] = {"None", "Gravity-type force (scalar particle quantity)", "Air Drag-type force (vector particle quantity)", "Attractor", "Repulsor", "Damped Spring", "Force back to ThetaBar"};

	enum ConstraintType {C_NONE, C_Circular};
	static const char* ConstraintType2str[] = {"None", "Circular/Bar-type constraint"};

	enum CollisionType {CL_NONE, FLAT, PARTICLE, SPHERE};
	static const char* CollisionType2str[] = {"None", "Flat surface", "Particle to particle", "Inside sphere"};

	enum SolverType {GROUND, EXP_E, MIDPOINT, RK3, RK4, IMP_E, TRAP, VERLET, RK4_G};
	static const char* SolverType2str[] = {"Ground_Truth", "Explicit_Euler", "Midpoint", "RK3", "RK4", "Implicit_Euler", "Trapezoidal", "Verlet", "Gen_RK4"};

	enum SceneType {BALL_DROP, SNOW_GLOBE, MPM_FLUID, CNSTR_1, NEWT_CRDL, CNSTR_4,CNSTR_4_JMP, INV_PEND, SOVLER_9, RK4_LAM_9, SEAWEED, MSPR_MTN_PROJ,MSPR_MTN_PRO2};
    static const char* SceneType2str[] = {"Ball Drop", "Snow Globe", "MPM Fluid", "1 Initial Constraint", "Newton's Cradle", "4 Init Constraints", "4 Cnstr Indy Jones", "Inverted Pendulum", "9 Solvers", "RK4 Gen Form Comp", "Seaweed", "Mass Spring RT", "Mass Spring 2" };

	static const int numCurrSystems = 13;

	static const Eigen::Vector3d gravVec = Eigen::Vector3d(0, -9.8, 0);		//used in ground truth eval

    static const int numSolvers =  9; //number of numeric solvers uses for simulation of falling particles : 9 + none

}//namespace
#endif