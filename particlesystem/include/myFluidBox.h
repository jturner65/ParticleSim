#ifndef myFluidBox_h
#define myFluidBox_h

#include <memory>
#include <array>
#include <map>
#include <Eigen/Dense>
#include "myParticle.h"
#include "myGlobConsts.h"

namespace particleSystem {
	//boundary cell
	class myFluidBndObj {
	public:
		myFluidBndObj(int _idx, const Eigen::Ref<const Eigen::Vector3d>& _normToCtr, double _mag) : idx(_idx), norm(_normToCtr), mag(_mag) {}
		~myFluidBndObj() {}

	public: //fields
		const int idx;
		Eigen::Vector3d norm;
		double mag;					//0-1, pct of cell within boundary

	};

	//template<std::size_t numCellX, std::size_t numCellY, std::size_t numCellZ>
	class myFluidBox {
	public:
		typedef particleSystem::aligned_containers<Eigen::Vector3d>::vector_t eignVecTyp;
		eignVecTyp vortVec, accelVecX, accelVecY, accelVecZ;


		myFluidBox();
		myFluidBox(int _numCellX, int _numCellY, int _numCellZ, double _diffusion, double _viscosity, double _deltaT, int _sviters, const Eigen::Ref<const Eigen::Vector3d>& _cellSz);
		virtual ~myFluidBox();
		inline void initVecs() {
			oldDensity = new double[numCells];
			density = new double[numCells];
			Vx = new double[numCells];
			Vy = new double[numCells];
			Vz = new double[numCells];
			Vx0 = new double[numCells];
			Vy0 = new double[numCells];
			Vz0 = new double[numCells];
			isOOB = new bool[numCells];
			vortVec.reserve(numCells);
			vortVec.resize(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { vortVec[i].setZero(); }
			accelVecX.reserve(numCells);
			accelVecX.resize(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { accelVecX[i].setZero(); }
			accelVecY.reserve(numCells);
			accelVecY.resize(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { accelVecY[i].setZero(); }
			accelVecZ.reserve(numCells);
			accelVecZ.resize(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { accelVecZ[i].setZero(); }
		}

		//if this fluidbox uses a mesh (like a sphere or some other mesh), set up boundary 
		//structure to hold idx of upper fwd left corner of cell containing bound and vector 
		//from center of cell to center of mesh, whO2Se magnitude corresponds 
		//amount of cube inside mesh.
		inline void initMeshBounds() {
			double* sq_dists = new double[numCells];
			double amtInside=0;
			for (unsigned int k = 0; k < numCellZ; ++k) {
				for (unsigned int j = 0; j < numCellY; ++j) {
					for (unsigned int i = 0; i < numCellX; ++i) {
						sq_dists[IX(i, j, k)] = getXYZSqDistFromCtr(i, j, k);
					}
				}
			}
			double minA = 88888888888, maxA = -888888888;
			for (unsigned int k = 0; k < numCellZ; ++k) {
				for (unsigned int j = 0; j < numCellY; ++j) {
					for (unsigned int i = 0; i < numCellX; ++i) {
						int idx = IX(i, j, k);
						if (checkIDXAtThreshold(radSq*.95, idx, i, j, k, sq_dists, amtInside)) {
							minA = (minA > amtInside ? amtInside : minA);
							maxA = (maxA < amtInside ? amtInside : maxA);
							sphereBnds[idx] = make_shared<myFluidBndObj>(idx, (center - getIdxLoc(idx)).normalized(), amtInside);			//set idx's in-pointing normal
						}
					}
				}
			}
			delete[] sq_dists;
			cout << "# boundary cells for mesh : " << sphereBnds.size() << " min/max amt inside : " << minA << "," << maxA << " center of mesh : " << E3VecToStr(center)<<"\n";
		}//initSphereBounds

		//returns true if cell given by x,y,z encapsulates a given threshold value from the center of the fluid box
		inline bool checkIDXAtThreshold(double thresh, unsigned int _idx, unsigned int _x, unsigned int _y, unsigned int _z, double* sq_dists, double& amtInside) {
			amtInside = 0;
			double minVal = 9999999999, maxVal = -9999999999, d;
			unsigned int idx = 0, x, y, z;
			for (unsigned int k = 0; k < 2; ++k) {
				z = (_z >= numCellZ-1 ? _z : _z + k);
				for (unsigned int j = 0; j < 2; ++j) {
					y = (_y >= numCellY-1 ? _y : _y + j);
					for (unsigned int i = 0; i < 2; ++i) {
						x = (_x >= numCellX-1 ? _x : _x + i);
						d = sq_dists[IX(x, y, z)];
						minVal = (minVal > d ? d : minVal);
						maxVal = (maxVal < d ? d : maxVal);
					}
				}
			}
			if (minVal >= thresh) {
				isOOB[_idx] = true;
				//sphereExtBnds.push_back(_idx);
				amtInside = 0;
				return false;
			}
			isOOB[_idx] = false;
			if ((maxVal < thresh) || (maxVal - minVal == 0)) { amtInside = 0;	return false; }
			amtInside = (thresh - minVal)/(maxVal - minVal);
			return true;
		}//checkIDXAtThreshold

		inline int getInBndsIDX(int _ckIdx, int idx) {	return isOOB[_ckIdx] ? idx : _ckIdx;}

		inline double getXYZSqDistFromCtr(int x, int y, int z) {
			Eigen::Vector3d res;
			res << cellSz(0)*(x - halfNmCellX), cellSz(1)*(y - halfNmCellY), cellSz(2)*(z - halfNmCellZ);			//res here is dist vals from center
			return res.dot(res);
		}
		//returns Squared distance from center of fluid box
		inline double getIdxSqDistFromCtr(int idx) {
			Eigen::Vector3d res(0, 0, 0);
			int x = 0, y = 0, z = 0;			//x,y,z idx's in 3d version of array
			invIX(idx, x, y, z);
			return getXYZSqDistFromCtr(x, y, z);
		}
		inline Eigen::Vector3d getVecToCenter(int x, int y, int z) {
			Eigen::Vector3d res(0, 0, 0);
			res << cellSz(0)*(x - halfNmCellX), cellSz(1)*(y - halfNmCellY), cellSz(2)*(z - halfNmCellZ);			//res here is dist vals from center
			return res;
		}

		//x,y,z coords in the world of cell
		inline Eigen::Vector3d getIdxLoc(int idx) {
			//Eigen::Vector3d res(0, 0, 0);
			int x = 0, y = 0, z = 0;			//x,y,z idx's in 3d version of array
			invIX(idx, x, y, z);
			//cout << "IDX : " << idx << " -> " << x << "," << y << "," << z << "\t halfNumCells : " << halfNmCellX << "," << halfNmCellY << "," << halfNmCellZ << "\n";
			//res << cellSz(0)*(x - halfNmCellX), cellSz(1)*(y - halfNmCellY), cellSz(2)*(z - halfNmCellZ);			//res here is dist vals from center
			Eigen::Vector3d res = getVecToCenter(x, y, z);		//res here is dist vals from center
			res += center;
			return res;
		}
		//normalized vector from cell @idx to center - length is 1 + difference in distance between min dist vert of cube and threshold/boundary
		inline Eigen::Vector3d getIdxCtrVec(int idx) {
			Eigen::Vector3d res = (center - getIdxLoc(idx)).normalized();
			return res;
		}

		inline void setStartLoc() { startLoc << center[0] - cellSz(0) * (halfNmCellX - .5), center[1] - cellSz(1) * (halfNmCellY - .5), center[2] - cellSz(2) * (halfNmCellZ - .5); }
		//inline void swap(vector<double>& a, vector<double>& b) { auto t = a;	a = b;	b = t; }
		//inline void swap(double* a, double* b) { double* t = a;	a = b;	b = t; }
		//set location of center of box in world coords
		void setCenter(const Eigen::Ref<const Eigen::Vector3d>& _ctr) {
			center<<_ctr; 
			setStartLoc();
			Eigen::Vector3d tmpHalfCell(halfNmCellX, halfNmCellY, halfNmCellZ);
			//tmpHalfCell << halfNmCellX, halfNmCellY, halfNmCellZ;
			ctrSzHalfNC = (center - tmpHalfCell.cwiseProduct(cellSz)).cwiseQuotient(cellSz);			//precalculated (center - sz*halfNumCell)/sz for force query for particles
			//cout << "Fluid Box Center : " << center(0) << "," << center(1) << "," << center(2) << "\tSize : " << cellSz(0) << "," << cellSz(1) << "," << cellSz(2) << "\n";
			//test if idx's to locations are working properly
			//for (unsigned int k = 0; k < numCellZ; k += (sz1i) / 2) {
			//	for (unsigned int j = 0; j < numCellY; j += (numCellY - 1) / 2) {
			//		for (unsigned int i = 0; i < numCellX; i += (sx1i) / 2) {
			//			Eigen::Vector3d res = getIdxLoc(IX(i, j, k));
			//			double res1 = getIdxSqDistFromCtr(IX(i, j, k));
			//			cout << "Res for i,j,k: " << i << "," << j << "," << k << " : " << res(0) << "," << res(1) << "," << res(2) << "\tSqDist from center : " << res1 << "\n";
			//		}
			//	}
			//}
		}//setCenter

		//set isMesh and init sphere boundaries if true
		inline void setIsMesh(bool val) {
			isMesh = val;
			if (isMesh) {
				initMeshBounds();
			}
		}
		//inline void set_bnd(int b, double* x) { if (isMesh) { if (b == 0) { set_bndDiffSphere(x); return; } else { set_bndSphere(b, x); return; } } set_bndCube(b, x); }
		//inline void set_bndSphere(int b, double* x);
		inline void set_bndCube(int b, double* x);
		inline int forceIDXBnd(int idx, int idxMax, int idxMin) { return (((idx > idxMin)&&(idx < idxMax) ? idx : (idx > idxMin) ? idxMax : idxMin)); }
		inline double forceIDXBndD(double val, double ubnd, double lbnd) { return (((val > lbnd) && (val < ubnd) ? val : (val  > lbnd) ? ubnd : lbnd)); }
		inline void addSource(double* cur, double* old) {			for (unsigned int i = 0; i < numCells; ++i) { if (isOOB[i]) { continue; } cur[i] += old[i]; }		}

		inline int IX(int x, int y, int z) {			return x + (y * numCellX) + (z * numCellXY);		}
		//given an index, return the x,y,z coords
		inline void invIX(int idx, int& x, int& y, int& z) {
			z = (int)(idx / numCellXY);
			int tmp = (idx % numCellXY);
			y = (int)(tmp/numCellX); 
			x = (int)(tmp % numCellX);
		}
		void resetOldVals();

		//vorticity confinement
		void vorticityConfinement(double* vortN);
		//vorticity particle method
		void vorticityParticles();

		//for mesh
		inline void set_bndDiffSphere(double* x);
		inline void set_bndSphere3(double* x, double* y, double* z);
		inline void set_bndDiffVec(eignVecTyp& vec);
		//inline void lin_solveSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double a, double c, int iter);

		//solve for all dirs simulataneously
		void lin_solveSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double a, double c);

		inline void diffuseSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double viscdiff){
			//void myFluidBox::diffuseSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double viscdiff, int iter) {
			double delVisc = (deltaT*viscdiff*numCells);
			lin_solveSphere(x, x0, y, y0, z, z0, delVisc, 1 + 6 * delVisc);
		}

		inline void advectSphere(double* _velx, double* _velx0, double* _vely, double* _vely0, double* _velz, double* _velz0);
		inline void projectSphere(double* velocX, double* velocY, double* velocZ, double* p, double* div);
		inline void diffSphDens(double* x, double* xOld, double viscdiff);
		inline void advSphDens(double* d, double* d0, double* velocX, double* velocY, double* velocZ);

		//for simple box
		//inline void lin_solve(int b, double* x, double* x0, double a, double c, int iter);
		inline void lin_solve(int b, double* x, double* x0, double a, double c) {
			//double cInv = 1.0 / c;
			for (unsigned int itr = 0; itr < slvIters; ++itr) {
				for (unsigned int k = 1; k < sz1i; ++k) {
					for (unsigned int j = 1; j < sy1i; ++j) {
						for (unsigned int i = 1; i < sx1i; ++i) {
							x[IX(i, j, k)] = (x0[IX(i, j, k)] + a *
								(x[IX(i + 1, j, k)]
									+ x[IX(i - 1, j, k)]
									+ x[IX(i, j + 1, k)]
									+ x[IX(i, j - 1, k)]
									+ x[IX(i, j, k + 1)]
									+ x[IX(i, j, k - 1)])) / c;
						}//for i
					}//for j
				}//for k
				set_bndCube(b, x);
			}//for itr
		}//lin solvers


//		inline void diffuse(int b, double* x, double* xOld, double viscdiff, int _numCells);
		inline void diffuse(int b, double* x, double* xOld, double viscdiff, int _numCells) {
			double delVisc = (deltaT*viscdiff*(_numCells));
			lin_solve(b, x, xOld, delVisc, 1 + 6 * delVisc);
		}
		inline void advect(int b, double* d, double* d0, double* velocX, double* velocY, double* velocZ);
		inline void project(double* velocX, double* velocY, double* velocZ, double* p, double* div);

		void myFluidBoxTimeStep();
		void myFluidSphereTimeStep();


		void myFluidBoxAddDensity(int x, int y, int z, double amount);
		void myFluidBoxAddForce(const Eigen::Ref<const Eigen::Vector3d>& cellLoc, const Eigen::Ref<const Eigen::Vector3d>& amount);
		Eigen::Vector3d getVelAtCell(const Eigen::Ref<const Eigen::Vector3d>&  cellLoc);

		std::string E3VecToStr(const Eigen::Ref<const Eigen::Vector3d>& vec) { stringstream ss; ss << vec(0) << "," << vec(1) << "," << vec(2); return ss.str(); }

		friend ostream& operator<<(ostream& out, const myFluidBox& fb) {
			out << "Fluidbox ID : " << fb.ID << " #cells in X:" << fb.numCellX << " #cells in Y:" << fb.numCellY << " #cells in Z:" << fb.numCellZ << " diff : " << fb.diff << " visc : " << fb.visc;
			out << "Center : " << fb.center(0) << ", " << fb.center(1) << ", " << fb.center(2) << "\n";
			//out << " size of s (array size): " << fb.s.size();
			out << endl;
			return out;
		}

	public:
		static unsigned int ID_gen;
		int ID;
		bool isMesh;								//fluid is inside a mesh
		double deltaT, diff, visc, radSq;			//radSq used for collision if sphere - is square of radius
		
		int halfNmCellX, halfNmCellY, halfNmCellZ;
		unsigned int numCellX, numCellY, numCellZ, numCellXY, numCells;	//corresponds to # of cells in each direction
		Eigen::Vector3d cellSz,						//size of cell in each dimension - defaults to 1
						startLoc,					//location to translate to for start of arrays, to render velocity vectors
						ctrSzHalfNC,				//precalculated (center - sz*halfNumCell)/sz for force query for particles
						center;

		//index in arrays of boundary object, mapped to ptr to fluid bnd obj
		typedef std::map<int, std::shared_ptr<myFluidBndObj>> sphereBndMap;		
		sphereBndMap sphereBnds;										//cells at boundary of sphere/shaped fluid box

		//use arrays for speed
		double *oldDensity, *density, 
			*Vx, *Vy, *Vz, 
			*Vx0, *Vy0, *Vz0;
		bool* isOOB;

	public : 
		//precomputed constants
		double
			vortEps,														//amount of vorticity allowed TODO make UI enterable
			hO2Sxd, hO2Syd, hO2Szd,									
			hSxd, hSyd, hSzd;
		unsigned int 
			memSetNumElems,													//holds size of arrays in bytes
			sz1i, sz2i,
			sx1i, sx2i,
			sy1i, sy2i,
			slvIters;															//# of reptitions of GS solver

		//sstd::vector<int> sphereExtBnds;									//cells where no fluid should go.
	};
}//partsys
#endif