#ifndef MyParticleWorld_h
#define MyParticleWorld_h

#include <vector>
#include <memory>
#include <deque>
#include "myGlobConsts.h"
#include <Eigen/Dense>
#include "ControlPanel.h"

namespace particleSystem{
    class myParticle;
    class myForce;
    class mySolver;
    class mySystem;
    class mySpring;
    class myCollider;
    class myConstraint;
}


using namespace std;
using namespace particleSystem;

class MyParticleWorld{
	public:
		//per system per frame per particle
		static vector<vector<vector<std::shared_ptr<myParticle>>>> systemParticleHolder;
		static vector<vector<vector<std::shared_ptr<myConstraint>>>> systemConstraintHolder;
		static vector<vector<vector<std::shared_ptr<mySpring>>>> systemSpringHolder;

		//systems for this project  - should be # of different scenes	
		static vector<std::shared_ptr<mySystem>> systems;

        static deque<Eigen::Vector3d> msClickVal;        
		static deque<Eigen::Vector3d> msRelVal;
        static deque<Eigen::Vector3d> msDragVal;
		
		static vector<bool> flags;
		static const int debug						= 0;					//init to false
		static const int msClicked					= 1;					//init to false			
		static const int msDragged					= 2;					//init to false			
		static const int msReleased					= 3;					//init to false		

		static const int numFlags = 4; 

		//static int curSolverIDX;
		static int shakeCountdown;
		static unsigned int numSnowFlakes;
		static int clickOnPartIDX;

        //for mass-spring rhombic triacontahedron
        static const double rt_fi;//golden ratio
        static vector<Eigen::Vector3d> rt_crnrPts;                         //locations  of corner points
        static const int numRTparts = 33;
        static const int numRTEdges = 90;
        static const int numRTMsclSpr = 32;                     //these all radiate out from the center
        static const int rt_edgePtIdxs[numRTEdges][2];                   //idxs in crnrPts vector of endpoints for every edge
        static const int cnstOppPts[numRTMsclSpr][2];                        //idxs of pairs of points that are opposite each other - has dupes, idx 0 is in ascending order

        //for mass spring creature 2
        static vector<Eigen::Vector3d> mass2_crnrPts;
        static const int numMS2Parts = 13;
        static const int numMS2Edges = 34;
        static const int numMSMsclSpr = 14;
        static const int mass2_edgePtIdxs[numMS2Edges][2];                      //idxs in crnrPts vector of endpoints for every edge
        static const int mass2_springPtIdxs[numMSMsclSpr][2];                      //idxs in crnrPts vector of endpoints for every edge

	    //use _usePartSolver means use a different solver for each particle
		static void buildRhTrHdrn(Eigen::Vector3d& sLoc, int curSystemIDX); //rhombic triacontahedron (30 sided die) - for mass spring motion project
		static void buildMsSprMtn2(Eigen::Vector3d& sLoc, int curSystemIDX); //some other mass-spring object

		static void initScene(int sceneNum, double deltaT, SolverType _solveType, int numP, int numF, int numC);
		static void handlePause(int curSystemIDX);
		static void handleTimeStep(int curSystemIDX);
};

#endif