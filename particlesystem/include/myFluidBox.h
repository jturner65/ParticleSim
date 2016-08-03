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
		myFluidBndObj(int _idx, const Eigen::Vector3d& _normToCtr, double _mag) : idx(_idx), norm(_normToCtr), mag(_mag) {}
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

		myFluidBox() :ID(++ID_gen) {}//std::cout<<"building fluid box ID : "<<ID<<std::endl;	}
		myFluidBox(int _numCellX, int _numCellY, int _numCellZ, double _diffusion, double _viscosity, double _deltaT, const Eigen::Vector3d& _cellSz) :
			ID(++ID_gen), numCellX(_numCellX), numCellY(_numCellY), numCellZ(_numCellZ), diff(_diffusion), visc(_viscosity), startLoc(0, 0, 0), sphereBnds(), ctrSzHalfNC(0, 0, 0), //sphereExtBnds(), 				//precalculated (center - sz*halfNumCell)/sz for force query for particles
			deltaT(_deltaT), center(0, 0, 0), numCellXY(numCellX * numCellY), numCells(numCellX * numCellY * numCellZ), cellSz(_cellSz), isMesh(false), radSq(0){
			initVecs();

			halfNmCellX = (numCellX / 2);
			halfNmCellY = (numCellY / 2);
			halfNmCellZ = (numCellZ / 2);

			//corresponds to x,y,z idx of "internal" array (inside single cube boundary layer
			sx1i = numCellX - 1;
			sy1i = numCellY - 1;
			sz1i = numCellZ - 1;
			//corresponds to x,y,z size of internal cube array of nodes used for sim
			sx2i = numCellX - 2;
			sy2i = numCellY - 2;
			sz2i = numCellZ - 2;
			hOSxd = 1.0 / (2.0* sx2i); hOSyd = 1.0 / (2.0* sy2i); hOSzd = 1.0 / (2.0* sz2i);
			hSxd = 0.5f* sx2i; hSyd = 0.5f* sy2i; hSzd = 0.5f* sz2i;
			memSetNumElems = sizeof(Vx0[0]) * numCells;
			//location to start rendering
			setStartLoc();
		}

		~myFluidBox() {
			delete[] oldDensity;
			delete[] density;
			delete[] isOOB;
			delete[] Vx;
			delete[] Vy;
			delete[] Vz;
			delete[] Vx0;
			delete[] Vy0;
			delete[] Vz0;
		}

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
			for (unsigned int i = 0; i < numCells; ++i) { vortVec[i].setZero(); }
			accelVecX.reserve(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { accelVecX[i].setZero(); }
			accelVecY.reserve(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { accelVecY[i].setZero(); }
			accelVecZ.reserve(numCells);
			for (unsigned int i = 0; i < numCells; ++i) { accelVecZ[i].setZero(); }
		}

		//if this fluidbox uses a mesh (like a sphere or some other mesh), set up boundary 
		//structure to hold idx of upper fwd left corner of cell containing bound and vector 
		//from center of cell to center of mesh, whose magnitude corresponds 
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
		inline bool checkIDXAtThreshold(double thresh, int _idx, int _x, int _y, int _z, double* sq_dists, double& amtInside) {
			amtInside = 0;
			double minVal = 9999999999, maxVal = -9999999999, d;
			int idx = 0, x, y, z;
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
		void setCenter(const Eigen::Vector3d& _ctr) { 
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
		inline void lin_solveSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double a, double c, unsigned int iter) {
			if (a == 0) {
				std::memcpy(&x, &x0, sizeof x0);
				std::memcpy(&y, &y0, sizeof y0);
				std::memcpy(&z, &z0, sizeof z0);
				set_bndSphere3(x, y, z);
			}
			else {
				unsigned int idx0, idx1, idx2, idx3, idx4, idx5, idx6;
				for (unsigned int itr = 0; itr < iter; ++itr) {
					for (unsigned int k = 1; k < sz1i; ++k) {
						for (unsigned int j = 1; j < sy1i; ++j) {
							for (unsigned int i = 1; i < sx1i; ++i) {
								idx0 = IX(i, j, k);
								if (isOOB[idx0]) { continue; }
								idx1 = IX(i + 1, j, k);
								idx2 = IX(i - 1, j, k);
								idx3 = IX(i, j + 1, k);
								idx4 = IX(i, j - 1, k);
								idx5 = IX(i, j, k + 1);
								idx6 = IX(i, j, k - 1);
								x[idx0] = (x0[idx0] + a * (x[idx1] + x[idx2] + x[idx3] + x[idx4] + x[idx5] + x[idx6])) / c;
								y[idx0] = (y0[idx0] + a * (y[idx1] + y[idx2] + y[idx3] + y[idx4] + y[idx5] + y[idx6])) / c;
								z[idx0] = (z0[idx0] + a * (z[idx1] + z[idx2] + z[idx3] + z[idx4] + z[idx5] + z[idx6])) / c;
							}//for i
						}//for j
					}//for k
					set_bndSphere3(x, y, z);
				}//for itr
			}//if a != 0
		}//lin_solveSphere


		inline void diffuseSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double viscdiff, int iter){
			//void myFluidBox::diffuseSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double viscdiff, int iter) {
			double delVisc = (deltaT*viscdiff*(numCells));
			lin_solveSphere(x, x0, y, y0, z, z0, delVisc, 1 + 6 * delVisc, iter);
		}

		inline void advectSphere(double* _velx, double* _velx0, double* _vely, double* _vely0, double* _velz, double* _velz0);
		inline void projectSphere(double* velocX, double* velocY, double* velocZ, double* p, double* div, int iter);
		inline void diffSphDens(double* x, double* xOld, double viscdiff, int iter, int _numCells);
		inline void advSphDens(double* d, double* d0, double* velocX, double* velocY, double* velocZ);

		//for simple box
		inline void lin_solve(int b, double* x, double* x0, double a, double c, int iter);
		inline void diffuse(int b, double* x, double* xOld, double viscdiff, int iter, int _numCells);
		inline void advect(int b, double* d, double* d0, double* velocX, double* velocY, double* velocZ);
		inline void project(double* velocX, double* velocY, double* velocZ, double* p, double* div, int iter);

		void myFluidBoxTimeStep();
		void myFluidSphereTimeStep();


		void myFluidBoxAddDensity(int x, int y, int z, double amount);
		void myFluidBoxAddForce(const Eigen::Vector3d& cellLoc, const Eigen::Vector3d& amount);
		Eigen::Vector3d getVelAtCell(const Eigen::Vector3d&  cellLoc);

		std::string E3VecToStr(const Eigen::Vector3d& vec) { stringstream ss; ss << vec(0) << "," << vec(1) << "," << vec(2); return ss.str(); }

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
		double hOSxd, hOSyd, hOSzd, hSxd, hSyd, hSzd;
		unsigned int 
			memSetNumElems,													//holds size of arrays in bytes
			sz1i, sz2i,
			sx1i, sx2i,
			sy1i, sy2i;

		//sstd::vector<int> sphereExtBnds;									//cells where no fluid should go.
	};
}//partsys
#endif