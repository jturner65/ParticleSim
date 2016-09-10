#include <string>
#include <vector>
#include <sstream>
#include <ctime>
#include <math.h>
#include <memory>
#include "MyParticleWorld.h"
#include "myParticle.h"
#include "mySystem.h"
#include "myForce.h"
#include "mySolver.h"
#include "myConstraint.h"
#include "mySpring.h"
using namespace particleSystem;
using namespace std;

extern int curSystemIDX;

//static vars
//int MyParticleWorld::curSolverIDX = -1;
vector<bool> MyParticleWorld::flags(MyParticleWorld::numFlags, false);

int MyParticleWorld::shakeCountdown = 0;
unsigned int MyParticleWorld::numSnowFlakes;
int MyParticleWorld::clickOnPartIDX = -1;
//Eigen::Vector3d MyParticleWorld::partCOM = Eigen::Vector3d(0,0,0);
//double MyParticleWorld::tsCntr = 0;
deque<Eigen::Vector3d>  MyParticleWorld::msDragVal(3, Eigen::Vector3d(0, 0, 0));
deque<Eigen::Vector3d>  MyParticleWorld::msRelVal(3, Eigen::Vector3d(0, 0, 0));
deque<Eigen::Vector3d>  MyParticleWorld::msClickVal(3, Eigen::Vector3d(0, 0, 0));

const double MyParticleWorld::rt_fi = ((1 + sqrt(5.0))/2.0);
//these are relative locations of all particles
vector<Eigen::Vector3d> MyParticleWorld::rt_crnrPts = [](double _fi)->vector<Eigen::Vector3d>{
    double c = .6;
    double fi = c*_fi;
    double fiSq = c*_fi*_fi;
    double fiCu = c*_fi*_fi*_fi;	
	vector<Eigen::Vector3d> pts(MyParticleWorld::numRTparts, Eigen::Vector3d(0, 0, 0));
    pts[0] =  Eigen::Vector3d( fiSq, 0, fiCu);
    pts[1] =  Eigen::Vector3d( 0, fi, fiCu);
    pts[2] =  Eigen::Vector3d(-fiSq, 0, fiCu);
    pts[3] =  Eigen::Vector3d( 0,-fi, fiCu);
    pts[4] =  Eigen::Vector3d ( fiSq, fiSq, fiSq);
    pts[5] =  Eigen::Vector3d ( 0, fiCu, fiSq);
    pts[6] =  Eigen::Vector3d (-fiSq, fiSq, fiSq);   
    pts[7] =  Eigen::Vector3d (-fiSq,-fiSq, fiSq);  
    pts[8] =  Eigen::Vector3d ( 0,-fiCu, fiSq); 
    pts[9] =  Eigen::Vector3d ( fiSq,-fiSq, fiSq);    
    pts[10] = Eigen::Vector3d ( fiCu, 0, fi);  
    pts[11] = Eigen::Vector3d (-fiCu, 0, fi);    
    pts[12] = Eigen::Vector3d ( fiCu, fiSq,0);   
    pts[13] = Eigen::Vector3d ( fi , fiCu, 0); 
    pts[14] = Eigen::Vector3d (-fi , fiCu, 0);   
    pts[15] = Eigen::Vector3d (-fiCu, fiSq, 0);   
    pts[16] = Eigen::Vector3d (-fiCu,-fiSq, 0);   
    pts[17] = Eigen::Vector3d (-fi ,-fiCu, 0); 
    pts[18] = Eigen::Vector3d ( fi ,-fiCu, 0);    
    pts[19] = Eigen::Vector3d ( fiCu,-fiSq, 0);   
    pts[20] = Eigen::Vector3d ( fiCu, 0,-fi);   
    pts[21] = Eigen::Vector3d (-fiCu, 0,-fi);  
    pts[22] = Eigen::Vector3d (fiSq, fiSq,-fiSq);   
    pts[23] = Eigen::Vector3d (0, fiCu,-fiSq); 
    pts[24] = Eigen::Vector3d (-fiSq, fiSq,-fiSq);  
    pts[25] = Eigen::Vector3d (-fiSq,-fiSq,-fiSq);   
    pts[26] = Eigen::Vector3d (0,-fiCu,-fiSq);  
    pts[27] = Eigen::Vector3d ( fiSq,-fiSq,-fiSq);  
    pts[28] = Eigen::Vector3d ( fiSq, 0,-fiCu);  
    pts[29] = Eigen::Vector3d (0, fi,-fiCu);  
    pts[30] = Eigen::Vector3d (-fiSq, 0,-fiCu); 
    pts[31] = Eigen::Vector3d (0,-fi,-fiCu);	
    pts[32] = Eigen::Vector3d (0,0,0);	            //center point for springs
    return pts;}(MyParticleWorld::rt_fi);

//these edges are constraints - we want springs attaching center particle to all outside particles
const int MyParticleWorld::rt_edgePtIdxs[numRTEdges][2] = { 
		{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 4}, {4, 5}, {1, 5},						
	{5, 6}, {6, 2}, {2, 11}, {2, 7}, {7, 8}, {8, 3}, {8, 9}, {0, 9},
	{0, 10}, {10, 12}, {12, 13}, {5, 13}, {5, 14}, {6, 15}, {11, 15}, {11, 16},
	{16, 7}, {9, 19}, {19, 10}, {4, 12}, {28, 29}, {29, 30}, {30, 31}, {31, 28},
	{28, 22}, {22, 23}, {23, 29}, {30, 24}, {30, 21}, {30, 25}, {31, 26}, {27, 28},
	{28, 20}, {20, 12}, {12, 22}, {23, 24}, {24, 15}, {15, 21}, {21, 16}, {16, 25},
	{25, 26}, {26, 27}, {27, 19}, {19, 20}, {13, 23}, {14, 23}, {14, 15}, {8, 18},
	{18, 26}, {26, 17}, {17, 8}, {18, 19}, {16, 17}
	////30 more edges that bisect romboids - adds 30 more
		,{1, 3}, {1, 4}, {1,6},  {3, 9}, {3, 7},
		{4,10}, {4,13},	{6,11}, {6,14}, 
		{7,11}, {7,17},	{9,10}, {9,18},
		{10,20},{11,21},{13,14},{13,22},
		{14,24},{17,18},{17,25},{18,27},
		{20,27},{20,22},{21,24},{21,25},
		{22,29},{24,29},{25,31},{27,31},{29,31}		
	};
const int MyParticleWorld::cnstOppPts[numRTMsclSpr][2] = {
    {0, 30},  {1, 31},  {2, 28},  {3, 29},  {4, 25},  {5, 26},  {6, 27},  {7, 22},  {8, 23},
    {9, 24},  {10, 21}, {11, 20}, {12, 16}, {13, 17}, {14, 18}, {15, 19}, 
    {16, 12}, {17, 13}, {18, 14}, {19, 15}, {20, 11}, {21, 10}, {22, 7},  {23, 8},  {24, 9}, 
    {25, 4},  {26, 5},  {27, 6},  {28, 2},  {29, 3},  {30, 0},  {31, 1}
};

vector<Eigen::Vector3d> MyParticleWorld::mass2_crnrPts = []()->vector<Eigen::Vector3d>{
    double mult = .5, 
        minX = -mult, minY = -mult, minZ = -mult,minX2 = 2*minX,minX3 = minX*3,
        maxX = mult, maxY = mult, maxZ = mult, maxX2 = 2*maxX, maxX3 = maxX*3;
    
	vector<Eigen::Vector3d> pts(MyParticleWorld::numMS2Parts, Eigen::Vector3d(0, 0, 0));
        pts[0] =  Eigen::Vector3d(minX, minY, minZ);//base of body
        pts[1] =  Eigen::Vector3d(maxX, minY, minZ);
        pts[2] =  Eigen::Vector3d(minX, minY, maxZ);
        pts[3] =  Eigen::Vector3d(maxX, minY, maxZ);
        pts[4] =  Eigen::Vector3d(0,maxY,0);          //point of body
        //left side
        pts[5] = Eigen::Vector3d(minX2, maxY, minZ);
        pts[6] = Eigen::Vector3d(minX2, maxY, maxZ);
        pts[7] = Eigen::Vector3d(minX3, minY, minZ);      //feet
        pts[8] = Eigen::Vector3d(minX3, minY, maxZ);
        //right side
        pts[9]  = Eigen::Vector3d(maxX2, maxY, minZ);
        pts[10] = Eigen::Vector3d(maxX2, maxY, maxZ);
        pts[11] = Eigen::Vector3d(maxX3, minY, minZ);     //feet
        pts[12] = Eigen::Vector3d(maxX3, minY, maxZ);

    return pts;}();

const int MyParticleWorld::mass2_edgePtIdxs[numMS2Edges][2] = { //constraints - 28 edges
    {0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,4}, {2,3}, {2,4}, {3,4}, //body pyramid
    {0,5}, {0,6}, {2,5}, {2,6}, {5,6}, {5,7}, {5,8}, {6,7}, {6,8}, {7,8}, //left side
    {1,9}, {1,10}, {3,9}, {3,10}, {9,10}, {9,11}, {9,12}, {10,11}, {10,12}, {11,12}, //right side
    {4,5}, {4,6}, {4,9}, {4,10}
};

const int MyParticleWorld::mass2_springPtIdxs[numMSMsclSpr][2] = {//springs - 8
    {0,7}, {0,8}, {2,7}, {2,8}, {4,7}, {4,8}, //left side - idx 4&5 are "up" springs
    {1,11}, {1,12}, {3,11}, {3,12},{4,11}, {4,12},     //right side springs idx 10&11 are "up" springs
    {7,11},  {8,12}                //knees together
};
 
//systems for each part
vector<std::shared_ptr<mySystem>> MyParticleWorld::systems(numCurrSystems);
vector<vector<vector<std::shared_ptr<myParticle>>>> MyParticleWorld::systemParticleHolder(numCurrSystems);
vector<vector<vector<std::shared_ptr<myConstraint>>>>  MyParticleWorld::systemConstraintHolder(numCurrSystems);
vector<vector<vector<std::shared_ptr<mySpring>>>>  MyParticleWorld::systemSpringHolder(numCurrSystems);

////copy code below

//calculate correct z value given 2 other vals between +/- snowGlobStartRad
//double MyParticleWorld::calcZ(double x, double y){return sqrt((snowGlobStartRad*snowGlobStartRad) - (x*x) - (y*y)) * ((rand() % 2 ==1) ? 1 : -1);}
//mass spring motion - final project for cs7492
void MyParticleWorld::buildRhTrHdrn(Eigen::Vector3d& sLoc, int curSystemIDX) {
	double ks = 500, kd = 40;
	stringstream ss;
	//forces       
	systems[curSystemIDX]->buildDefForces(std::string("MassSpring"), -.05);
	//particles
	for (unsigned int i = 0; i < MyParticleWorld::numRTparts; ++i) {
		//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, MyParticleWorld::rt_crnrPts[i] - Eigen::Vector3d(0, 1, 15), Eigen::Vector3d(0, 0, 0),RK4_G));
		systems[curSystemIDX]->addParticle(1.0, MyParticleWorld::rt_crnrPts[i] - Eigen::Vector3d(0, 1, 15), Eigen::Vector3d(0, -2, 0), RK4_G);
		systems[curSystemIDX]->p.back()->setOrigMass(((i == MyParticleWorld::numRTparts - 1)) ? 100 : 1);     //only last one is heavier
		//systems[curSystemIDX]->p.back()->setVelocity(Eigen::Vector3d(0, -2, 0));          //some initial velocity
	}
	//constraints
	for (int i = 0; i<MyParticleWorld::numRTEdges; ++i) {
		systems[curSystemIDX]->buildAndSetCnstrnts(MyParticleWorld::rt_edgePtIdxs[i][0], MyParticleWorld::rt_edgePtIdxs[i][1], -1, systems[curSystemIDX]->p[MyParticleWorld::rt_edgePtIdxs[i][1]]->getPosition());
	}
	int numXBodyCnstr = (int)MyParticleWorld::numRTMsclSpr / 2;
	for (int i = 0; i<numXBodyCnstr; ++i) {
		systems[curSystemIDX]->buildAndSetCnstrnts(MyParticleWorld::cnstOppPts[i][0], MyParticleWorld::cnstOppPts[i][1], -1, systems[curSystemIDX]->p[MyParticleWorld::cnstOppPts[i][1]]->getPosition());
	}
	////springs
	for (int i = 0; i<MyParticleWorld::numRTMsclSpr; ++i) {
		//  buildAndSetNonPathCnstrntVals(i,numRTMsclSpr, p[i]->position[0]);
		////make springs from last particle to all others here
		ss.str("");
		ss << "RT Spring_" << i << "_a_" << systems[curSystemIDX]->p[i]->ID << "__b_" << systems[curSystemIDX]->p[MyParticleWorld::numRTparts - 1]->ID;
		//cout<<ss.str()<<endl;
		double rLen = (systems[curSystemIDX]->p[i]->getPosition() - systems[curSystemIDX]->p[MyParticleWorld::numRTparts - 1]->getPosition()).norm();
		systems[curSystemIDX]->spr.push_back(std::make_shared<mySpring>(ss.str(), systems[curSystemIDX]->p[i], systems[curSystemIDX]->p[MyParticleWorld::numRTparts - 1], ks, kd, rLen));
	}
	//collider ground
	systems[curSystemIDX]->buildGndCollider(.4, .9, 500, -5, 500, std::string("Ground_MassSpring"), Eigen::Vector3d(0, -5, -12));// , Eigen::Vector3d(0, -8, -12));
}

void MyParticleWorld::buildMsSprMtn2(Eigen::Vector3d& sLoc, int curSystemIDX) {
	double ks = 100, kd = 150;
	stringstream ss;
	//forces       
	std::shared_ptr<myForce>tmpForce1(new myForce("MassSpring_Grav", (gravVec), S_SCALAR)); //add gravity force
	systems[curSystemIDX]->f.push_back(tmpForce1);
	std::shared_ptr<myForce>tmpForce2(new myForce("MassSpring_Drag", -.05, S_VECTOR));		//add fluid drag 
	systems[curSystemIDX]->f.push_back(tmpForce2);
	//particles
	double mass = 1;
	for (unsigned int i = 0; i < MyParticleWorld::numMS2Parts; ++i) {
		//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, MyParticleWorld::mass2_crnrPts[i] - Eigen::Vector3d(0, 1, 10), Eigen::Vector3d(0, 0, 0),RK4));
		systems[curSystemIDX]->addParticle(1, MyParticleWorld::mass2_crnrPts[i] - Eigen::Vector3d(0, 1, 10), Eigen::Vector3d(0, -2, 0),RK4);
		//systems[curSystemIDX]->p.back()->setOrigMass(mass);        //only last one is heavier
		//systems[curSystemIDX]->p.back()->setVelocity(Eigen::Vector3d(0, -2, 0));          //some initial velocity
	}
	//constraints
	for (int i = 0; i<MyParticleWorld::numMS2Edges; ++i) {
		systems[curSystemIDX]->buildAndSetCnstrnts(MyParticleWorld::mass2_edgePtIdxs[i][0], MyParticleWorld::mass2_edgePtIdxs[i][1], -1, systems[curSystemIDX]->p[MyParticleWorld::mass2_edgePtIdxs[i][1]]->getPosition());
	}
	////springs //mass2_springPtIdxs
	for (int i = 0; i<MyParticleWorld::numMSMsclSpr; ++i) {
		ss.str("");
		ss << "MassSpring2 Spring_" << i << "_a_" << systems[curSystemIDX]->p[MyParticleWorld::mass2_springPtIdxs[i][0]]->ID << "__b_" << systems[curSystemIDX]->p[MyParticleWorld::mass2_springPtIdxs[i][1]]->ID;
		//cout<<ss.str()<<endl;
		double rLen = (systems[curSystemIDX]->p[MyParticleWorld::mass2_springPtIdxs[i][0]]->getPosition() - systems[curSystemIDX]->p[MyParticleWorld::mass2_springPtIdxs[i][1]]->getPosition()).norm();
		systems[curSystemIDX]->spr.push_back(std::make_shared<mySpring>(ss.str(), systems[curSystemIDX]->p[MyParticleWorld::mass2_springPtIdxs[i][0]], systems[curSystemIDX]->p[MyParticleWorld::mass2_springPtIdxs[i][1]], ks, kd, rLen));
	}
	for (int i = 4; i<6; ++i) { systems[curSystemIDX]->spr[i]->scaleSpring(2);		systems[curSystemIDX]->spr[i + 6]->scaleSpring(2); }
	for (int i = 12; i<MyParticleWorld::numMSMsclSpr; ++i) { systems[curSystemIDX]->spr[i]->scaleSpring(1.5); }
	//collider ground
	systems[curSystemIDX]->buildGndCollider(0, .9, 500, -3, 500, std::string("Ground_MassSpring"), Eigen::Vector3d(0, -3, -12));//, Eigen::Vector3d(0, -8, -12));
}

void MyParticleWorld::initScene(int curSystemIDX, double _deltaT, SolverType _solveType, int numP, int numF, int numC) {
	srand((unsigned int)(time(NULL)));
    myParticle::ID_gen = 0;     //reinit all particle numbering so always start at 0 for each system
	myConstraint::ID_gen = 0;
	myCollider::ID_gen = 0;
    //tsCntr = 0;             //re-init animation control
    stringstream ss;
    ss<<"System : "<<SceneType2str[curSystemIDX];
	//initialize scene with starting conditions	
	systems[curSystemIDX] = std::make_shared<mySystem>(ss.str(), _deltaT, numP, numF, numC);

	switch(curSystemIDX){
		case BALL_DROP : {		//ball drop
			for (unsigned int i = 0; i < numSolvers; ++i){//1 particle per solver
				//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, Eigen::Vector3d(-3.0 + i, 8, -12), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), particleSystem::SolverType(i)));
				systems[curSystemIDX]->addParticle(1, Eigen::Vector3d(-3.0 + i, 8, -12), particleSystem::SolverType(i));
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), 0);//no drag
			//collider ground
			systems[curSystemIDX]->buildGndCollider(.7, 0, 10, -8, 10, std::string("GroundPart1"), Eigen::Vector3d(0, -8, -12));//, Eigen::Vector3d(0, -8, -12));
			break;}

		case SNOW_GLOBE : {		//snow globe
			double x, y, z, theta, phi, gndLocY = -6;
			numSnowFlakes = 6000;// mUI->mControl->getCurPartCount();
			double startRad = (.9 * snowGlobRad);// (snowGlobStartRad);
			for (unsigned int i = 0; i < numSnowFlakes; ++i){
				z = ((2.0*(rand() / (1.0*RAND_MAX))) - 1.0) * startRad;
				phi = (rand()/(1.0*RAND_MAX)) * twoPI * 1.0;
				theta = asin(z / startRad);
				x = startRad * cos(theta) * cos(phi);
				y = std::max(startRad * cos(theta) * sin(phi), gndLocY+.1);
				//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1.0, Eigen::Vector3d(x, y, z - (snowGlobRad + distFromSnowGlobe)), Eigen::Vector3d(0, 0, 0),RK4));
				systems[curSystemIDX]->addParticle(1.0, Eigen::Vector3d(x, y, z - (snowGlobRad + distFromSnowGlobe)),RK4);
				systems[curSystemIDX]->p.back()->setColor(1, 1, 1, 1);
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -4.0);
			systems[curSystemIDX]->buildGndCollider(.00001, 0, 10, gndLocY, 10, std::string("GroundPart1"), Eigen::Vector3d(0, gndLocY, -12));//, Eigen::Vector3d(0, gndLocY, -12));
			systems[curSystemIDX]->buildGlobeCollider(.00001, 0, snowGlobRad, distFromSnowGlobe);
			int numCells = (snowGlobRad + 1) * 4;
			double cellDim = (snowGlobRad + 1) / (.5*numCells);
			systems[curSystemIDX]->buildFluidBox(numCells, numCells, numCells, .0006, 0.0001, systems[curSystemIDX]->colliders.back()->center, Eigen::Vector3d(cellDim, cellDim, cellDim));
			systems[curSystemIDX]->fluidBox->radSq = snowGlobRad * snowGlobRad;
			systems[curSystemIDX]->fluidBox->setIsMesh(true);
			systems[curSystemIDX]->flags[systems[curSystemIDX]->useMassMat] = false;
			cout << "Fluid globe : " << *(systems[curSystemIDX]->fluidBox) << "\n";
			break;}
		case MPM_FLUID: {//MPM fluid TODO
			double x, y, z, theta, phi, gndLocY = -6;
			numSnowFlakes = 6000;// mUI->mControl->getCurPartCount();
			double startRad = (.9 * snowGlobRad);// (snowGlobStartRad);
			for (unsigned int i = 0; i < numSnowFlakes; ++i) {
				z = ((2.0*(rand() / (1.0*RAND_MAX))) - 1.0) * startRad;
				phi = (rand() / (1.0*RAND_MAX)) * twoPI * 1.0;
				theta = asin(z / startRad);
				x = startRad * cos(theta) * cos(phi);
				y = std::max(startRad * cos(theta) * sin(phi), gndLocY + .1);
				//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1.0, Eigen::Vector3d(x, y, z - (snowGlobRad + distFromSnowGlobe)), Eigen::Vector3d(0, 0, 0),RK4));
				systems[curSystemIDX]->addParticle(1.0, Eigen::Vector3d(x, y, z - (snowGlobRad + distFromSnowGlobe)),RK4);
				systems[curSystemIDX]->p.back()->setColor(1, 1, 1, 1);
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -4.0);
			systems[curSystemIDX]->buildGndCollider(.00001, 0, 10, gndLocY, 10, std::string("GroundPart1"), Eigen::Vector3d(0, gndLocY, -12));//, Eigen::Vector3d(0, gndLocY, -12));
			systems[curSystemIDX]->buildGlobeCollider(.00001, 0, snowGlobRad, distFromSnowGlobe);
			int numCells = (snowGlobRad + 1) * 4;
			double cellDim = (snowGlobRad + 1) / (.5*numCells);
			systems[curSystemIDX]->buildFluidBox(numCells, numCells, numCells, .0006, 0.0001, systems[curSystemIDX]->colliders.back()->center, Eigen::Vector3d(cellDim, cellDim, cellDim));
			systems[curSystemIDX]->fluidBox->radSq = snowGlobRad * snowGlobRad;
			systems[curSystemIDX]->fluidBox->setIsMesh(true);
			systems[curSystemIDX]->flags[systems[curSystemIDX]->useMassMat] = false;
			cout << "MPM Fluid globe : " << *(systems[curSystemIDX]->fluidBox) << "\n";
			break; }
		case CNSTR_1 : {		//tinkertoy/bead on wire assignment
			//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, Eigen::Vector3d(2, 0, -5), Eigen::Vector3d(0, 0, 0),RK4));
			systems[curSystemIDX]->addParticle(1, Eigen::Vector3d(2, 0, -5),RK4);
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.3);
			//initialize first contraint - circle wire
			systems[curSystemIDX]->buildAndSetCnstrnts(0, 0, 2, Eigen::Vector3d(0, 0, -5));

			break;}
		case NEWT_CRDL : {//multiple circular constraints - newton cradle
			//initialize first particle
			for (unsigned int i = 0; i < 4; ++i){
				systems[curSystemIDX]->addParticle(1, Eigen::Vector3d(2.5 + i / 2.0, 0, -5),RK4);
				systems[curSystemIDX]->addParticle(1, Eigen::Vector3d(-2.5 - i / 2.0, 0, -5),RK4);
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);

			for (unsigned int i = 0; i < 4; ++i){
				systems[curSystemIDX]->buildAndSetCnstrnts((i * 2), (i * 2), 2, Eigen::Vector3d(.5 + i / 2.0, 0, -5));
				systems[curSystemIDX]->buildAndSetCnstrnts((i * 2) + 1, (i * 2) + 1, 2, Eigen::Vector3d(-.5 - i / 2.0, 0, -5));
			}
			//build global structures needed for constraint calculations
			systems[curSystemIDX]->buildCnstrntStruct(false);
			break;}
		case CNSTR_4 : {//perpetual motion - particles cause each other to keep moving, or speed up
			std::vector<Eigen::Vector3d> clocList(4);	//init part positions
			std::vector<Eigen::Vector3d> locList(4);	//init part positions
			clocList[0] = Eigen::Vector3d(2, 2, -5);
			clocList[1] = Eigen::Vector3d(-2, 2, -5);
			clocList[2] = Eigen::Vector3d(2, -2, -5);
			clocList[3] = Eigen::Vector3d(-2, -2, -5);
			locList[0] = clocList[0] + Eigen::Vector3d((2 * sin(PI / 6.0)), (2 * cos(PI / 6.05)), 0);
			locList[1] = clocList[1] + Eigen::Vector3d((2 * sin(11 * PI / 6.0)), (2 * cos(11 * PI / 6.0)), 0);
			locList[2] = clocList[2] + Eigen::Vector3d((2 * sin(11 * PI / 6.0)), (2 * cos(11 * PI / 6.0)), 0);
			locList[3] = clocList[3] + Eigen::Vector3d((2 * sin(PI / 6.0)), (2 * cos(PI / 4.0)), 0);
			for (int i = 0; i < 4; ++i) {
				systems[curSystemIDX]->addParticle(1, locList[i],RK4);
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);
			for (unsigned int i = 0; i < 4; ++i) {
				systems[curSystemIDX]->buildAndSetCnstrnts(i, i, 2, clocList[i]);
			}
			//build global structures needed for constraint calculations
			systems[curSystemIDX]->buildCnstrntStruct(false);
			break;}
		case CNSTR_4_JMP :	{//multi particle rollercoaster - particles should jump from constraint to constraint, based on the results of the constraint calculations
			for(int i = 0; i<10; ++i){
				double angle = twoPI * (rand()/(1.0 * RAND_MAX));
				//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, Eigen::Vector3d((2 * sin(angle)) + (2 * ((rand() % 2 == 0) ? 1 : -1)), (2 * cos(angle)) + 2, -5), Eigen::Vector3d(0, 0, 0),RK4));
				systems[curSystemIDX]->addParticle(1, Eigen::Vector3d((2 * sin(angle)) + (2 * ((rand() % 2 == 0) ? 1 : -1)), (2 * cos(angle)) + 2, -5),RK4);
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);

			for (int i = 0; i<10; ++i) { systems[curSystemIDX]->buildRollerCoasterConstraints(i, 2); }
			//build global structures needed for constraint calculations - use part 5 version
			systems[curSystemIDX]->buildCnstrntStruct(true);
			break;}
        case INV_PEND : {//  inverted pendulum
			systems[curSystemIDX]->buildInvPend(Eigen::Vector3d(0, -1, -5));   //location of particles
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);
			systems[curSystemIDX]->buildGndCollider(.7, 0, 20, -2, 20, std::string("GroundPart7"), Eigen::Vector3d(0, -2, -5));//, Eigen::Vector3d(0, -2, -5));
			systems[curSystemIDX]->partCOM = systems[curSystemIDX]->calcAndSetCOM(0, systems[curSystemIDX]->p.size()); 

            break;}
		case SOVLER_9 : {//varying integrators on 9 beads on wires
			for (unsigned int i = 0; i < 3; ++i){
				for (unsigned int j = 0; j < 3; ++j) {
					systems[curSystemIDX]->addParticle(1, Eigen::Vector3d(-4 + (i*3.0), (j*-3.0) + 3.0, -5), particleSystem::SolverType((i * 3 + j) % (numSolvers)));
					//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, Eigen::Vector3d(-4 + (i*3.0), (j*-3.0) + 3.0, -5), Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 0), particleSystem::SolverType((i * 3 + j) % (numSolvers))));
				}
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);
			for (unsigned int i = 0; i < 3; ++i){
				for (unsigned int j = 0; j < 3; ++j) {
					systems[curSystemIDX]->buildAndSetCnstrnts((i * 3) + j, (i * 3) + j, 1, Eigen::Vector3d(-3 + (i*3.0), (j*-3.0) + 3.0, -5));
					
				}
			}
			systems[curSystemIDX]->buildCnstrntStruct(false);
			break;}
        case RK4_LAM_9 : {//general rk4 bead on wire w/varying lambda
			for (unsigned int i = 0; i < 3; ++i) {
				for (unsigned int j = 0; j < 3; ++j) {
					systems[curSystemIDX]->addParticle(1, Eigen::Vector3d(-4 + (i*3.0), (j*-3.0) + 3.0, -5),RK4_G);
					//systems[curSystemIDX]->p.push_back(std::make_shared<myParticle>(1, Eigen::Vector3d(-4 + (i*3.0), (j*-3.0) + 3.0, -5), Eigen::Vector3d(0, 0, 0),RK4_G));
				}
			}
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);
			for (unsigned int i = 0; i < 3; ++i) {
				for (unsigned int j = 0; j < 3; ++j) {
					systems[curSystemIDX]->buildAndSetCnstrnts((i * 3) + j, (i * 3) + j, 1, Eigen::Vector3d(-3 + (i*3.0), (j*-3.0) + 3.0, -5));

				}
			}
			systems[curSystemIDX]->buildCnstrntStruct(false);
            break;}
        case SEAWEED : {//  seaweed field - 3 x 4 array of particle "poles" inv pendulum

			for(int xLoc = -2; xLoc < 3; ++xLoc){
                for (int zLoc = -6; zLoc < -3; ++zLoc){
					systems[curSystemIDX]->buildInvPend(Eigen::Vector3d(xLoc, -1, zLoc));   //location of first particle
                }
            }
			//add gravity force
			systems[curSystemIDX]->buildDefForces(std::string(SceneType2str[curSystemIDX]), -.1);
			//collider ground
			systems[curSystemIDX]->buildGndCollider(.7, 0, 20, -2, 20, std::string("GroundPart10"), Eigen::Vector3d(0, -2, -5));//, Eigen::Vector3d(0, -2, -5));
			int numCells = 14;				//#of cells on each side of fluid box, including bounds
			int seaWeedWidth = 6;			//width of "seaweed" patch
			double cellDim = (seaWeedWidth+1) / (.5*numCells);
			systems[curSystemIDX]->buildFluidBox(numCells, numCells, numCells, .0012, 0.01, Eigen::Vector3d(0, 2, -5), Eigen::Vector3d(cellDim, cellDim, cellDim));
			systems[curSystemIDX]->fluidBox->radSq = 6 * 6;
            break;}
        case MSPR_MTN_PROJ : {//        mass spring motion - final project for cs7492
			buildRhTrHdrn(Eigen::Vector3d(0, -1, -15), curSystemIDX);
			systemSpringHolder[curSystemIDX].push_back(systems[curSystemIDX]->spr);
            break;}
        case MSPR_MTN_PRO2 : {//        mass spring motion - final project for cs7492
			buildMsSprMtn2(Eigen::Vector3d(0, -1, -15), curSystemIDX);
			systemSpringHolder[curSystemIDX].push_back(systems[curSystemIDX]->spr);
            break;}
		default : {break;}
	}//switch
    //Make sure to call this after any particles have been added to a system
	systems[curSystemIDX]->initMassSystem();
	//display-related structures
	systemParticleHolder[curSystemIDX].push_back(systems[curSystemIDX]->p);
	systemConstraintHolder[curSystemIDX].push_back(systems[curSystemIDX]->c);
}//initscene

void MyParticleWorld::handleTimeStep(int curSystemIDX) {
	auto sys = systems[curSystemIDX];
	switch (curSystemIDX) {
		case BALL_DROP : {//galileo
			sys->applyForcesToSystem();
			if (MyParticleWorld::flags[MyParticleWorld::msDragged]) { 
				Eigen::Vector3d msDiff = MyParticleWorld::msDragVal[0] - MyParticleWorld::msDragVal[1];
				sys->addForcesToTinkerToys(20, MyParticleWorld::flags[MyParticleWorld::msDragged], msDiff);
			}
			sys->handlePartCldrCollision();//checks for, and handles, collisions for all particles to see if they hit the ground
			sys->handlePartPartCollision();		//handle particle-particle interraction for balldrop
			sys->invokeSolverDerivEval();//moves updated state back to particles - ready to be drawn
			break;}
		case SNOW_GLOBE : {//snowglobe
			bool tmpBool = sys->handleSnowGlobeTimeStep(.1, MyParticleWorld::flags[MyParticleWorld::msDragged], MyParticleWorld::msDragVal[0], MyParticleWorld::msDragVal[1]);
			MyParticleWorld::flags[MyParticleWorld::msDragged] = tmpBool;
			break;}
		case MPM_FLUID: {//MPM fluid TODO



			break; }
		case CNSTR_1://single bead onwire	
		case NEWT_CRDL  ://clackers newton cradle
        case CNSTR_4_JMP ://4 beads jump constraints
		case SOVLER_9://varying integrators with bead on wire 
		case RK4_LAM_9://varying lambda general rk4 with bead on wire
			{
				Eigen::Vector3d msDiff = MyParticleWorld::msDragVal[0] - MyParticleWorld::msDragVal[1];
				sys->handleTTTimeStep(20, MyParticleWorld::flags[MyParticleWorld::msDragged], curSystemIDX == CNSTR_4_JMP, msDiff);
			break;}

		case CNSTR_4://4 beads on wire - each tries to influence others
			{
			Eigen::Vector3d msDiff = MyParticleWorld::msDragVal[0] - MyParticleWorld::msDragVal[1];
			sys->handleTTTimeStepRepel(20, MyParticleWorld::flags[MyParticleWorld::msDragged], false, msDiff);
			break; }

		case MSPR_MTN_PROJ :    // mass spring motion
		case MSPR_MTN_PRO2:	{
			Eigen::Vector3d msDiff = MyParticleWorld::msDragVal[0] - MyParticleWorld::msDragVal[1];
			sys->handleMassSpringTimeStep(curSystemIDX == MSPR_MTN_PRO2, MyParticleWorld::flags[MyParticleWorld::msClicked], MyParticleWorld::flags[MyParticleWorld::msDragged], numMSMsclSpr, msDiff);
			systemSpringHolder[curSystemIDX].push_back(sys->spr);
			break;}
        case SEAWEED :
        case INV_PEND ://inv pendulum
            { 				
				double fmult = ((curSystemIDX == SEAWEED) ? .1 : 20000);
				MyParticleWorld::flags[MyParticleWorld::msDragged] = sys->handleInvPendTimeStep(fmult, (curSystemIDX == SEAWEED), false, MyParticleWorld::flags[MyParticleWorld::msDragged], (curSystemIDX != SEAWEED), MyParticleWorld::msDragVal[0], MyParticleWorld::msDragVal[1]);
			break;}
 	    default : {return;}//exit method if no appropriate current scene number
	}//switch

	systemParticleHolder[curSystemIDX].push_back(sys->p);
	systemConstraintHolder[curSystemIDX].push_back(sys->c);
}//handleTimeStep

//set sceneModded == true if changed - called during pause/non-simulation
void MyParticleWorld::handlePause(int curSystemIDX) {
	bool sceneModded = false;           //set true only if particle or constraint added
	bool sceneModdedClk = false;           //set true only if particle or constraint added
	bool sceneModdedRel = false;           //set true only if particle or constraint added
	auto sys = systems[curSystemIDX];
	switch (curSystemIDX) {
		case BALL_DROP : {		//ball drop
			if (MyParticleWorld::flags[MyParticleWorld::msClicked]) {
				systems[curSystemIDX]->addParticle(1.0, Eigen::Vector3d(msClickVal[0][0] / 46, msClickVal[0][1] / 35, -12),RK4);
				
				MyParticleWorld::flags[MyParticleWorld::msClicked] = false;
                sceneModded = true;
			}//MyParticleWorld::flags[MyParticleWorld::msClicked]
			break;	 }
		case CNSTR_1  : 
		case NEWT_CRDL: 
		case CNSTR_4  : 
        case INV_PEND : 		//bead on wire, inv pendulum
            {      
			if (MyParticleWorld::flags[MyParticleWorld::msClicked]) {
				Eigen::Vector3d mClickLoc = Eigen::Vector3d(msClickVal[0][0] / 86, msClickVal[0][1] / 86, -5);
				sceneModdedClk = sys->handlePauseDrawClick(clickOnPartIDX, mClickLoc, false);
				MyParticleWorld::flags[MyParticleWorld::msClicked] = false;
			}//MyParticleWorld::flags[MyParticleWorld::msClicked]				
			MyParticleWorld::flags[MyParticleWorld::msDragged] = false;
			if(MyParticleWorld::flags[MyParticleWorld::msReleased]){
				MyParticleWorld::flags[MyParticleWorld::msReleased] = false;
				Eigen::Vector3d mRelLoc = Eigen::Vector3d(msRelVal[0][0] / 86, msRelVal[0][1] / 86, -5);
				sceneModdedRel = sys->handlePauseDrawRel(clickOnPartIDX, mRelLoc, false);
			}
			sceneModded = (sceneModdedClk || sceneModdedRel);
			break;	 }
		case MPM_FLUID: {//MPM fluid TODO



			break; }
		case MSPR_MTN_PROJ: {break; }//do nothing for now
		case MSPR_MTN_PRO2: {break; }
		case SNOW_GLOBE: {break; }          //no snow globe handling
    	case CNSTR_4_JMP : {break;	 }          //don't add particles to constraint-jumpers
        case SOVLER_9 : {break;    }           //9 wire integrator exp 
        case RK4_LAM_9 : {break;    }           //9 wire rk4 exp : don't add particles to 9-wire experiments
        case SEAWEED : {break;}
  	}//switch
    if(sceneModded){
		systemParticleHolder[curSystemIDX].push_back(sys->p);
		systemConstraintHolder[curSystemIDX].push_back(sys->c);
        if((curSystemIDX == MSPR_MTN_PROJ) || (curSystemIDX == MSPR_MTN_PRO2)){            systemSpringHolder[curSystemIDX].push_back(sys->spr);}
		sys->initMassSystem();				//call after particles added to/removed from system - recalcs mass matrix
	}

}//handlePause
