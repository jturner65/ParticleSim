#include <vector>
#include "..\include\myfluidBox.h"
#include <Eigen/Dense>

using namespace std;

namespace particleSystem{
	unsigned int myFluidBox::ID_gen = 0;

	myFluidBox::myFluidBox() :ID(++ID_gen) {}//std::cout<<"building fluid box ID : "<<ID<<std::endl;	}
	myFluidBox::myFluidBox(int _numCellX, int _numCellY, int _numCellZ, double _diffusion, double _viscosity, double _deltaT, int _sviters, const Eigen::Ref<const Eigen::Vector3d>& _cellSz) :
		ID(++ID_gen), numCellX(_numCellX), numCellY(_numCellY), numCellZ(_numCellZ), diff(_diffusion), visc(_viscosity), slvIters(_sviters), startLoc(0, 0, 0), sphereBnds(), ctrSzHalfNC(0, 0, 0), //sphereExtBnds(), 				//precalculated (center - sz*halfNumCell)/sz for force query for particles
		deltaT(_deltaT), center(0, 0, 0), numCellXY(numCellX * numCellY), numCells(numCellX * numCellY * numCellZ), cellSz(_cellSz), isMesh(false), radSq(0) {

		vortEps = deltaT * .01;			//TODO allow for UI input
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
		hO2Sxd = 1.0 / (2.0* sx2i); hO2Syd = 1.0 / (2.0* sy2i); hO2Szd = 1.0 / (2.0* sz2i);
		hSxd = 0.5f* sx2i; hSyd = 0.5f* sy2i; hSzd = 0.5f* sz2i;
		memSetNumElems = sizeof(Vx0[0]) * numCells;
		//location to start rendering
		setStartLoc();
	}

	myFluidBox::~myFluidBox() {
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
	void myFluidBox::set_bndCube(int b, double* x) {
		//int N = size;//for reading ease
		//int sz1 = numCellZ - 1, sz2 = numCellZ - 2,
		//	sx1 = numCellX - 1, sx2 = numCellX - 2,
		//	sy1 = sy1i, sy2 = numCellY - 2;
		if (b == 0) {//diffusion
			for (unsigned int j = 1; j < sy1i; ++j) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, j, 0)] = x[IX(i, j, 1)]; x[IX(i, j, sz1i)] = x[IX(i, j, sz2i)];}}
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, 0, k)] = x[IX(i, 1, k)]; x[IX(i, sy1i, k)] = x[IX(i, sy2i, k)];}}
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int j = 1; j < sy1i; ++j) {	x[IX(0, j, k)] = x[IX(1, j, k)]; x[IX(sx1i, j, k)] = x[IX(sx2i, j, k)];}}
		}
		else if (b == 1) {//reflecting in x
			for (unsigned int j = 1; j < sy1i; ++j) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, j, 0)] = x[IX(i, j, 1)];  x[IX(i, j, sz1i)] = x[IX(i, j, sz2i)]; } }
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, 0, k)] = x[IX(i, 1, k)];  x[IX(i, sy1i, k)] = x[IX(i, sy2i, k)]; } }
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int j = 1; j < sy1i; ++j) {	x[IX(0, j, k)] = -x[IX(1, j, k)]; x[IX(sx1i, j, k)] = -x[IX(sx2i, j, k)];} }
		}
		else if (b == 2) {//reflecting in y
			for (unsigned int j = 1; j < sy1i; ++j) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, j, 0)] = x[IX(i, j, 1)]; x[IX(i, j, sz1i)] = x[IX(i, j, sz2i)]; } }
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, 0, k)] = -x[IX(i, 1, k)];x[IX(i, sy1i, k)] = -x[IX(i, sy2i, k)];} }
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int j = 1; j < sy1i; ++j) {	x[IX(0, j, k)] = x[IX(1, j, k)]; x[IX(sx1i, j, k)] = x[IX(sx2i, j, k)]; } }
		}
		else if (b == 3) {//refelecting in z
			for (unsigned int j = 1; j < sy1i; ++j) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, j, 0)] = -x[IX(i, j, 1)]; x[IX(i, j, sz1i)] = -x[IX(i, j, sz2i)]; } }
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int i = 1; i < sx1i; ++i) {	x[IX(i, 0, k)] = x[IX(i, 1, k)];  x[IX(i, sy1i, k)] = x[IX(i, sy2i, k)]; } }
			for (unsigned int k = 1; k < sz1i; ++k) {for (unsigned int j = 1; j < sy1i; ++j) {	x[IX(0, j, k)] = x[IX(1, j, k)];  x[IX(sx1i, j, k)] = x[IX(sx2i, j, k)]; } }
		}

		// edges
		for (unsigned int i = 1; i < sx1i; ++i) {
			x[IX(i, 0, 0)] = 0.5f * (x[IX(i, 1, 0)] + x[IX(i, 0, 1)]);
			x[IX(i, sy1i, 0)] = 0.5f * (x[IX(i, sy2i, 0)] + x[IX(i, sy1i, 1)]);
			x[IX(i, 0, sz1i)] = 0.5f * (x[IX(i, 1, sz1i)] + x[IX(i, 0, sz2i)]);
			x[IX(i, sy1i, sz1i)] = 0.5f * (x[IX(i, sy2i, sz1i)] + x[IX(i, sy1i, sz2i)]);
		}

		for (unsigned int j = 1; j < sy1i; ++j) {
			x[IX(0, j, 0)] = 0.5f * (x[IX(1, j, 0)] + x[IX(0, j, 1)]);
			x[IX(sx1i, j, 0)] = 0.5f * (x[IX(sx2i, j, 0)] + x[IX(sx1i, j, 1)]);
			x[IX(0, j, sz1i)] = 0.5f * (x[IX(1, j, sz1i)] + x[IX(0, j, sz2i)]);
			x[IX(sx1i, j, sz1i)] = 0.5f * (x[IX(sx2i, j, sz1i)] + x[IX(sx1i, j, sz2i)]);
		}

		for (unsigned int k = 1; k < sz1i; ++k) {
			x[IX(0, 0, k)] = 0.5f * (x[IX(0, 1, k)] + x[IX(1, 0, k)]);
			x[IX(0, sy1i, k)] = 0.5f * (x[IX(0, sy2i, k)] + x[IX(1, sy1i, k)]);
			x[IX(sx1i, 0, k)] = 0.5f * (x[IX(sx1i, 1, k)] + x[IX(sx2i, 0, k)]);
			x[IX(sx1i, sy1i, k)] = 0.5f * (x[IX(sx1i, sy2i, k)] + x[IX(sx2i, sy1i, k)]);
		}

		// corners
		double calcVal = (1 / 3.0);
		x[IX(0, 0, 0)] = calcVal*(x[IX(0, 1, 0)] + x[IX(1, 0, 0)] + x[IX(0, 0, 1)]);
		x[IX(sx1i, 0, 0)] = calcVal*(x[IX(sx2i, 0, 0)] + x[IX(sx1i, 1, 0)] + x[IX(sx1i, 0, 1)]);
		x[IX(0, sy1i, 0)] = calcVal*(x[IX(0, sy2i, 0)] + x[IX(1, sy1i, 0)] + x[IX(0, sy1i, 1)]);
		x[IX(sx1i, sy1i, 0)] = calcVal*(x[IX(sx1i, sy2i, 0)] + x[IX(sx2i, sy1i, 0)] + x[IX(sx1i, sy1i, 1)]);

		x[IX(0, 0, sz1i)] = calcVal*(x[IX(0, 1, sz1i)] + x[IX(1, 0, sz1i)] + x[IX(0, 0, sz2i)]);
		x[IX(sx1i, 0, sz1i)] = calcVal*(x[IX(sx2i, 0, sz1i)] + x[IX(sx1i, 1, sz1i)] + x[IX(sx1i, 0, sz2i)]);
		x[IX(0, sy1i, sz1i)] = calcVal*(x[IX(0, sy2i, sz1i)] + x[IX(1, sy1i, sz1i)] + x[IX(0, sy1i, sz2i)]);
		x[IX(sx1i, sy1i, sz1i)] = calcVal*(x[IX(sx1i, sy2i, sz1i)] + x[IX(sx2i, sy1i, sz1i)] + x[IX(sx1i, sy1i, sz2i)]);
	}


	//handle advection for passed arrays of velocities
	void myFluidBox::advect(int b, double* d, double* d0, double* velocX, double* velocY, double* velocZ) {
		double s0, s1, t0, t1, u0, u1;
		//double tmp1, tmp2, tmp3;
		double x, y, z;

		//double idouble = 1, jdouble = 1, kdouble = 1;
		int i0, i1, j0, j1, k0, k1;

		double dtx = deltaT * (sx2i),
			dty = deltaT * (sy2i),
			dtz = deltaT * (sz2i);

		int IXidx;

		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					IXidx = IX(i, j, k);

					x = i - dtx * velocX[IXidx];
					y = j - dty * velocY[IXidx];
					z = k - dtz * velocZ[IXidx];
					x = forceIDXBndD(x, sx1i + 0.5, .5);
					i0 = floor(x);	i1 = i0 + 1.0;  s1 = x - i0;	s0 = 1.0 - s1;
					y = forceIDXBndD(y, sy1i + 0.5, .5);
					j0 = floor(y);	j1 = j0 + 1.0;	t1 = y - j0;	t0 = 1.0 - t1;
					z = forceIDXBndD(z, sz1i + 0.5, .5);
					k0 = floor(z);	k1 = k0 + 1.0;	u1 = z - k0;	u0 = 1.0 - u1;

					d[IXidx] =
						s0 * (t0 * (u0 * d0[IX(i0, j0, k0)]
						+ u1 * d0[IX(i0, j0, k1)])
						+ (t1 * (u0 * d0[IX(i0, j1, k0)]
						+ u1 * d0[IX(i0, j1, k1)])))
						+ s1 * (t0 * (u0 * d0[IX(i1, j0, k0)]
						+ u1 * d0[IX(i1, j0, k1)])
						+ (t1 * (u0 * d0[IX(i1, j1, k0)]
						+ u1 * d0[IX(i1, j1, k1)])));
				}//for i
			}//for j
		}//for k
		set_bndCube(b, d);
	}

	void myFluidBox::project(double* velocX, double* velocY, double* velocZ, double* p, double* div) {
		//double sx2 = sx2i, sy2 = numCellY - 2, sz2 = numCellZ - 2;
		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					div[IX(i, j, k)] = -(
						  (velocX[IX(i + 1, j, k)] - velocX[IX(i - 1, j, k)])*hO2Sxd
						+ (velocY[IX(i, j + 1, k)] - velocY[IX(i, j - 1, k)])*hO2Syd
						+ (velocZ[IX(i, j, k + 1)] - velocZ[IX(i, j, k - 1)])*hO2Szd
						);               
					p[IX(i, j, k)] = 0;
				}
			}
		}
		set_bndCube(0, div);
		//set_bnd(0, p);
		lin_solve(0, p, div, 1, 6);

		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					velocX[IX(i, j, k)] -= hSxd *(p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
					velocY[IX(i, j, k)] -= hSyd *(p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
					velocZ[IX(i, j, k)] -= hSzd *(p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
				}
			}
		}
		set_bndCube(1, velocX);
		set_bndCube(2, velocY);
		set_bndCube(3, velocZ);
	}

	//void myFluidBox::diffuse(int b, double* x, double* xOld, double viscdiff, int _numCells) {
	//	double delVisc = (deltaT*viscdiff*(numCells));
	//	lin_solve(b, x, xOld, delVisc, 1 + 6 * delVisc);
	//}

	void myFluidBox::myFluidBoxTimeStep() {
		addSource(Vx, Vx0);
		addSource(Vy, Vy0);
		addSource(Vz, Vz0);
		addSource(density, oldDensity);
		//if (isMesh) { myFluidSphereTimeStep(); return; }

		std::swap(Vx, Vx0); std::swap(Vy, Vy0); std::swap(Vz, Vz0);

		diffuse(1, Vx, Vx0, visc, numCellX);
		diffuse(2, Vy, Vy0, visc, numCellY);
		diffuse(3, Vz, Vz0, visc, numCellZ);

		project(Vx, Vy, Vz, Vx0, Vy0);

		std::swap(Vx, Vx0); std::swap(Vy, Vy0); std::swap(Vz, Vz0);

		advect(1, Vx, Vx0, Vx0, Vy0, Vz0);
		advect(2, Vy, Vy0, Vx0, Vy0, Vz0);
		advect(3, Vz, Vz0, Vx0, Vy0, Vz0);

		project(Vx, Vy, Vz, Vx0, Vy0);

		std::swap(density, oldDensity);
		diffuse(0, density, oldDensity, diff, numCellX);
		std::swap(density, oldDensity);
		advect(0, density, oldDensity, Vx0, Vy0, Vz0);
		resetOldVals();
	}//myFluidBoxTimeStep

	/////////////sphere stuff

	 //timestepping for sphere, to handle 3d bounds together instead of 1 dim at a time
	void myFluidBox::myFluidSphereTimeStep() {
		//addSource(Vx, Vx0);
		//addSource(Vy, Vy0);
		//addSource(Vz, Vz0);
		//addSource(density, oldDensity);

		std::swap(Vx, Vx0); std::swap(Vy, Vy0); std::swap(Vz, Vz0);
		double delVisc = (deltaT*visc*numCells);
		//diffusion of velocity
		lin_solveSphere(Vx, Vx0, Vy, Vy0, Vz, Vz0, delVisc, 1 + 6 * delVisc);

		projectSphere(Vx, Vy, Vz, Vx0, Vy0);

		std::swap(Vx, Vx0); std::swap(Vy, Vy0); std::swap(Vz, Vz0);
		advectSphere(Vx, Vx0, Vy, Vy0, Vz, Vz0);

		projectSphere(Vx, Vy, Vz, Vx0, Vy0);

		std::swap(density, oldDensity);
		diffSphDens(density, oldDensity, diff);
		std::swap(density, oldDensity);
		advSphDens(density, oldDensity, Vx, Vy, Vz);
		//vort confine here
		vorticityConfinement(oldDensity);
		//vorticity particle method
		//vorticityParticles();

		resetOldVals();
	}//myFluidSphereTimeStep


	//handle advection for passed arrays of velocities
	void myFluidBox::advectSphere(double* _velx, double* _velx0, double* _vely, double* _vely0, double* _velz, double* _velz0) {
		double s0, s1, t0, t1, u0, u1,x, y, z;

		int i0, i1, j0, j1, k0, k1;
		//precalced idx's
		int idx0, idx1, idx2, idx3, idx4, idx5, idx6, idx7;

		double dtx = deltaT * (sx2i),
			dty = deltaT * (sy2i),
			dtz = deltaT * (sz2i);
		int IXidx;

		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					IXidx = IX(i, j, k);
					if (isOOB[IXidx]) { continue; }

					x = i - dtx * _velx0[IXidx];
					x = forceIDXBndD(x, sx2i + 0.5, .5);
					i0 = floor(x);	i1 = i0 + 1.0;  s1 = x - i0;	s0 = 1.0 - s1;

					y = j - dty * _vely0[IXidx];
					y = forceIDXBndD(y, sy2i + 0.5, .5);
					j0 = floor(y);	j1 = j0 + 1.0;	t1 = y - j0;	t0 = 1.0 - t1;

					z = k - dtz * _velz0[IXidx];
					z = forceIDXBndD(z, sz2i + 0.5, .5);
					k0 = floor(z);	k1 = k0 + 1.0;	u1 = z - k0;	u0 = 1.0 - u1;

					idx0 = IX(i0, j0, k0);idx1 = IX(i0, j0, k1);idx2 = IX(i0, j1, k0);idx3 = IX(i0, j1, k1);
					idx4 = IX(i1, j0, k0);idx5 = IX(i1, j0, k1);idx6 = IX(i1, j1, k0);idx7 = IX(i1, j1, k1);

					_velx[IXidx] =	
						s0 * (t0 * (u0 * _velx0[idx0] + u1 * _velx0[idx1]) + (t1 * (u0 * _velx0[idx2] + u1 * _velx0[idx3]))) +
						s1 * (t0 * (u0 * _velx0[idx4] + u1 * _velx0[idx5]) + (t1 * (u0 * _velx0[idx6] + u1 * _velx0[idx7])));

					_vely[IXidx] =
						s0 * (t0 * (u0 * _vely0[idx0] + u1 * _vely0[idx1]) + (t1 * (u0 * _vely0[idx2] + u1 * _vely0[idx3]))) +
						s1 * (t0 * (u0 * _vely0[idx4] + u1 * _vely0[idx5]) + (t1 * (u0 * _vely0[idx6] + u1 * _vely0[idx7])));

					_velz[IXidx] =
						s0 * (t0 * (u0 * _velz0[idx0] + u1 * _velz0[idx1]) + (t1 * (u0 * _velz0[idx2] + u1 * _velz0[idx3]))) +
						s1 * (t0 * (u0 * _velz0[idx4] + u1 * _velz0[idx5]) + (t1 * (u0 * _velz0[idx6] + u1 * _velz0[idx7])));
				}//for i
			}//for j
		}//for k
		set_bndSphere3(_velx, _vely, _velz);
	}//advectSphere

	//diffuse density - 1 d but use sphere bnds
	void myFluidBox::diffSphDens(double* x, double* x0, double viscdiff) {
		double a = (deltaT*viscdiff*(numCells)), c = (1 + 6 * a);
		int idx;
		if(a==0){
			std::memcpy(&x, &x0, sizeof x0);
			set_bndDiffSphere(x);
		} else {
			for (unsigned int itr = 0; itr < slvIters; ++itr) {
				for (unsigned int k = 1; k < sz1i; ++k) {
					for (unsigned int j = 1; j < sy1i; ++j) {
						for (unsigned int i = 1; i < sx1i; ++i) {
							idx = IX(i, j, k);
							if (isOOB[idx]) { continue; }
							x[idx] = (x0[idx] + a *
								(x[IX(i + 1, j, k)]
									+ x[IX(i - 1, j, k)]
									+ x[IX(i, j + 1, k)]
									+ x[IX(i, j - 1, k)]
									+ x[IX(i, j, k + 1)]
									+ x[IX(i, j, k - 1)])) / c;
						}//for i
					}//for j
				}//for k
				set_bndDiffSphere(x);
			}
		}//for itr
	}//diffSphDens
	//advect density through sphere, using sphere bounds
	void myFluidBox::advSphDens(double* d, double* d0, double* velocX, double* velocY, double* velocZ) {
		double s0, s1, t0, t1, u0, u1, x, y, z, dtx = deltaT * (sx2i),dty = deltaT * (sy2i),dtz = deltaT * (sz2i);
		unsigned int i0, i1, j0, j1, k0, k1, i, j, k, IXidx;

		for (k = 1; k < sz1i; ++k) {
			for (j = 1; j < sy1i; ++j) {
				for (i = 1; i < sx1i; ++i) {
					IXidx = IX(i, j, k);
					if (isOOB[IXidx]) { continue; }
					x = i - dtx * velocX[IXidx];
					y = j - dty * velocY[IXidx];
					z = k - dtz * velocZ[IXidx];
					x = forceIDXBndD(x, sx2i + 0.5, .5);
					i0 = floor(x);	i1 = i0 + 1.0;  s1 = x - i0;	s0 = 1.0 - s1;

					y = forceIDXBndD(y, sy2i + 0.5, .5);
					j0 = floor(y);	j1 = j0 + 1.0;	t1 = y - j0;	t0 = 1.0 - t1;

					z = forceIDXBndD(z, sz2i + 0.5, .5);
					k0 = floor(z);	k1 = k0 + 1.0;	u1 = z - k0;	u0 = 1.0 - u1;

					d[IXidx] =
						s0 * (t0 * (u0 * d0[IX(i0, j0, k0)]	+ u1 * d0[IX(i0, j0, k1)]) + (t1 * (u0 * d0[IX(i0, j1, k0)] + u1 * d0[IX(i0, j1, k1)]))) + 
						s1 * (t0 * (u0 * d0[IX(i1, j0, k0)]	+ u1 * d0[IX(i1, j0, k1)]) + (t1 * (u0 * d0[IX(i1, j1, k0)]	+ u1 * d0[IX(i1, j1, k1)])));
				}//for i
			}//for j
		}//for k
		set_bndDiffSphere(d);
	}//advSphDens

	void myFluidBox::lin_solveSphere(double* x, double* x0, double* y, double* y0, double* z, double* z0, double a, double c) {
		if (a == 0) {
			std::memcpy(&x, &x0, sizeof x0);
			std::memcpy(&y, &y0, sizeof y0);
			std::memcpy(&z, &z0, sizeof z0);
			set_bndSphere3(x, y, z);
		}
		else {
			unsigned int idx0, idx1, idx2, idx3, idx4, idx5, idx6;
			for (unsigned int itr = 0; itr < slvIters; ++itr) {
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

	//vorticity confinement - add back vorticity details lost through numerical dissipation
	//vortN is unused array to hold calcs
	void myFluidBox::vorticityConfinement(double* vortN) {
		unsigned int idx, idx_ijp1k, idx_ijm1k, idx_ip1jk, idx_im1jk, idx_ijkm1, idx_ijkp1;
		//double vortEps = deltaT * .01;	//TODO change to allow for user input
		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					idx = IX(i, j, k);
					if (isOOB[idx]) { continue; }
					idx_ip1jk = IX(i + 1, j, k); idx_im1jk = IX(i - 1, j, k); idx_ijp1k = IX(i, j + 1, k); idx_ijm1k = IX(i, j - 1, k); idx_ijkp1 = IX(i, j, k + 1); idx_ijkm1 = IX(i, j, k - 1);
					//curl operation del cross u -> partial z w/respect to y is the finite diff of the z vels across the y coords
					vortVec[idx] << 
						((Vy[idx_ijkp1] - Vy[idx_ijkm1]) * hO2Syd) - ((Vz[idx_ijp1k] - Vz[idx_ijm1k]) * hO2Szd),
						((Vz[idx_ip1jk] - Vz[idx_im1jk]) * hO2Szd) - ((Vx[idx_ijkp1] - Vx[idx_ijkm1]) * hO2Sxd),
						((Vx[idx_ijp1k] - Vx[idx_ijm1k]) * hO2Sxd) - ((Vy[idx_ip1jk] - Vy[idx_im1jk]) * hO2Syd);
				}
			}
		}
		set_bndDiffVec(vortVec);
		for (idx = 0; idx < numCells; ++idx) { vortN[idx] = (isOOB[idx]) ?  0 : vortVec[idx].norm(); }

		Eigen::Vector3d eta, vf; eta.setZero();	vf.setZero();
		
		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {	
					idx = IX(i, j, k);
					if (vortN[idx] < .0000001) {	continue;}
					vortVec[idx].normalize();
					eta << ((vortN[IX(i + 1, j, k)] - vortN[IX(i - 1, j, k)]) * hO2Sxd), ((vortN[IX(i, j + 1, k)] - vortN[IX(i, j - 1, k)]) * hO2Syd), ( (vortN[IX(i, j, k + 1)] - vortN[IX(i, j, k - 1)]) * hO2Szd);
					eta.normalize();
					vf = vortEps * (eta.cross(vortVec[idx]));
					//cout << "Vx " << idx << " before :  " << Vx[idx] << " vortN : " << vortN[idx] << " invDivX : " << invDivX << " eta : " << eta(0) << "," << eta(1) << "," << eta(2) << " vf : " << vf(0) << "," << vf(1) << "," << vf(2) << " vort : " << vort(0) << "," << vort(1) << "," << vort(2);
					Vx[idx] += vf(0) * sx2i;
					Vy[idx] += vf(1) * sy2i;
					Vz[idx] += vf(2) * sz2i;	
					//cout << "Vx " << idx << " after :  " << Vx[idx]<<endl;
				}
			}
		}
		set_bndSphere3(Vx, Vy, Vz);
	}//vorticityConfinement

	//vorticity particle method TODO
	void myFluidBox::vorticityParticles() {
		int idx, idx_ijp1k, idx_ijm1k, idx_ip1jk, idx_im1jk, idx_ijkm1, idx_ijkp1;
		//find accelerations via finite diff
		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					idx = IX(i, j, k);
					if (isOOB[idx]) { continue; }
					idx_ip1jk = IX(i + 1, j, k);idx_im1jk = IX(i - 1, j, k);idx_ijp1k = IX(i, j + 1, k);idx_ijm1k = IX(i, j - 1, k);idx_ijkp1 = IX(i, j, k + 1);idx_ijkm1 = IX(i, j, k - 1);
												
					accelVecX[idx] << (Vx[idx_ip1jk] - Vx[idx_im1jk]) * hO2Sxd, (Vx[idx_ijp1k] - Vx[idx_ijm1k]) * hO2Syd, (Vx[idx_ijkp1] - Vx[idx_ijkm1]) * hO2Szd;//interpAccel(Vx, idx_ip1jk, idx_im1jk, idx_ijp1k, idx_ijm1k, idx_ijkp1, idx_ijkm1);
					accelVecY[idx] << (Vy[idx_ip1jk] - Vy[idx_im1jk]) * hO2Sxd, (Vy[idx_ijp1k] - Vy[idx_ijm1k]) * hO2Syd, (Vy[idx_ijkp1] - Vy[idx_ijkm1]) * hO2Szd;//interpAccel(Vy, idx_ip1jk, idx_im1jk, idx_ijp1k, idx_ijm1k, idx_ijkp1, idx_ijkm1);
					accelVecZ[idx] << (Vz[idx_ip1jk] - Vz[idx_im1jk]) * hO2Sxd, (Vz[idx_ijp1k] - Vz[idx_ijm1k]) * hO2Syd, (Vz[idx_ijkp1] - Vz[idx_ijkm1]) * hO2Szd;//interpAccel(Vz, idx_ip1jk, idx_im1jk, idx_ijp1k, idx_ijm1k, idx_ijkp1, idx_ijkm1);

				}
			}
		}
		set_bndDiffVec(accelVecX);
		set_bndDiffVec(accelVecY);
		set_bndDiffVec(accelVecZ);

		//vort particle code here TODO

	}//vorticityParticles

	void myFluidBox::projectSphere(double* velocX, double* velocY, double* velocZ, double* p, double* div) {
		unsigned int idx;
		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					idx = IX(i, j, k);
					p[idx] = 0;
					if (isOOB[idx]) {	div[idx] = 0;	continue; }
					//when using mac grid, change to * 1/sxi instead of 1/2sxi
					div[idx] = -(
						(velocX[IX(i + 1, j, k)] - velocX[IX(i - 1, j, k)]) * hO2Sxd
						+(velocY[IX(i, j + 1, k)] - velocY[IX(i, j - 1, k)]) * hO2Syd
						+(velocZ[IX(i, j, k + 1)] - velocZ[IX(i, j, k - 1)]) * hO2Szd
						);
				}
			}
		}
		set_bndDiffSphere(div);
		set_bndDiffSphere(p);
		for (unsigned int itr = 0; itr < slvIters; ++itr) {
			for (unsigned int k = 1; k < sz1i; ++k) {
				for (unsigned int j = 1; j < sy1i; ++j) {
					for (unsigned int i = 1; i < sx1i; ++i) {
						idx = IX(i, j, k);
						if (isOOB[idx]) { continue; }
						p[idx] = (div[idx] + (p[IX(i + 1, j, k)] + p[IX(i - 1, j, k)] + p[IX(i, j + 1, k)] + p[IX(i, j - 1, k)] + p[IX(i, j, k + 1)] + p[IX(i, j, k - 1)])) / 6.0;
					}//for i
				}//for j
			}//for k
			set_bndDiffSphere(p);
		}//for each iteration
		
		for (unsigned int k = 1; k < sz1i; ++k) {
			for (unsigned int j = 1; j < sy1i; ++j) {
				for (unsigned int i = 1; i < sx1i; ++i) {
					idx = IX(i, j, k);
					if (isOOB[idx]) { continue; }
					velocX[idx] -= hSxd *(p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
					velocY[idx] -= hSyd *(p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
					velocZ[idx] -= hSzd *(p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
				}
			}
		}
		set_bndSphere3(velocX, velocY, velocZ);
	}//projectSphere


	void myFluidBox::set_bndDiffSphere(double* x) {
		for (sphereBndMap::iterator it = sphereBnds.begin(); it != sphereBnds.end(); ++it) { x[it->first] *= it->second->mag; }//scale to amt of cube in bounds
	}
	void myFluidBox::set_bndDiffVec(eignVecTyp& egVec) {
		Eigen::Vector3d  velNorm(0, 0, 0);
		for (sphereBndMap::iterator it = sphereBnds.begin(); it != sphereBnds.end(); ++it) {
			//velNorm = (egVec[it->first].dot(it->second->norm)* it->second->mag)  *it->second->norm;			//velocity in the normal direction toward center of sphere
			velNorm = (egVec[it->first].dot(it->second->norm))  *it->second->norm;			//velocity in the normal direction toward center of sphere
			egVec[it->first] -= velNorm;													//remove velocity component in opposite direction of normal
			egVec[it->first] *= velNorm.norm();												//amplify remaining component by same amount - increase tangent velocity

		}
	}
	//address boundary layer values
	void myFluidBox::set_bndSphere3(double* x, double* y, double* z) {
		//int vIdx = b - 1, xIdx;
		Eigen::Vector3d velVec(0, 0, 0), velNorm(0, 0, 0);
		for (sphereBndMap::iterator it = sphereBnds.begin(); it != sphereBnds.end(); ++it) {
			velVec << x[it->first], y[it->first], z[it->first];
			//velNorm = (velVec.dot(it->second->norm)* it->second->mag)  *it->second->norm;			//velocity in the normal direction toward center of sphere
			velNorm = (velVec.dot(it->second->norm))  *it->second->norm;							//velocity in the normal direction toward center of sphere
																									//subtract velocity in direction of sphere wall
			x[it->first] -= velNorm(0);
			y[it->first] -= velNorm(1);
			z[it->first] -= velNorm(2);
			//scale result (tangent dir) by lost velocity magnitude - don't want to lose velocity, just redirect it
			x[it->first] *= velNorm.norm();
			y[it->first] *= velNorm.norm();
			z[it->first] *= velNorm.norm();
		}
	}//set_bndSphere

	void myFluidBox::myFluidBoxAddDensity(int x, int y, int z, double amount) { density[IX(x, y, z)] += amount; }

	void myFluidBox::resetOldVals() {
		//all elems same type
		//int numElems = sizeof(Vx0[0]) * numCells;
		memset(Vx0, 0, memSetNumElems);
		memset(Vy0, 0, memSetNumElems);
		memset(Vz0, 0, memSetNumElems);
		memset(oldDensity, 0, memSetNumElems);
	}

	void myFluidBox::myFluidBoxAddForce(const Eigen::Ref<const Eigen::Vector3d>& cellLoc, const Eigen::Ref<const Eigen::Vector3d>& amount) {
		//cout<<"force addition location in cube :("<< cellLoc (0)<<","<< cellLoc(1) <<","<< cellLoc(2) <<")"<<endl;

		int idx = IX(forceIDXBnd((int)(cellLoc(0)), sx1i,0),
					 forceIDXBnd((int)(cellLoc(1)), sy1i,0),
					 forceIDXBnd((int)(cellLoc(2)), sz1i,0));
		//static int iters = 0;
		//cout << "add force in fluidbox @idx : " << idx << " iter : " << iters++ << "\n";
		//Vx0[idx] = amount(0);
		//Vy0[idx] = amount(1);
		//Vz0[idx] = amount(2);
		Vx[idx] += amount(0);
		Vy[idx] += amount(1);
		Vz[idx] += amount(2);
	}

	//cellloc is a particle position - cell idx is going to be floor of each coord
	Eigen::Vector3d myFluidBox::getVelAtCell(const Eigen::Ref<const Eigen::Vector3d>& testLoc) {
		//ctrSzHalfNC == (ctr - (halfNumCell * cellSz))/cellSz
		//cout<<ctrSzHalfNC << "\n";
		Eigen::Vector3d tmpTestLoc = testLoc.cwiseQuotient(cellSz),
			testLocInFluid = tmpTestLoc - ctrSzHalfNC;
		
		//double locX = (testLoc(0) - center(0)) / cellSz(0) + halfNmCellX,
		//	locY = (testLoc(1) - center(1)) / cellSz(1) + halfNmCellY,
		//	locZ = (testLoc(2) - center(2)) / cellSz(2) + halfNmCellZ;
		int intLocX = (int)testLocInFluid(0),
			intLocY = (int)testLocInFluid(1),
			intLocZ = (int)testLocInFluid(2);

		double interpX = testLocInFluid(0) - intLocX, interpM1X = 1 - interpX,
			   interpY = testLocInFluid(1) - intLocY, interpM1Y = 1 - interpY,
			   interpZ = testLocInFluid(2) - intLocZ, interpM1Z = 1 - interpZ;
		//bound idx's
		//int tX[2], tY[2], tZ[2];
		//tX[0] = forceIDXBnd(intLocX, sx1i),
		//tY[0] = forceIDXBnd(intLocY, sy1i),
		//tZ[0] = forceIDXBnd(intLocZ, sz1i),
		//tX[1] = forceIDXBnd(intLocX + 1, sx1i),
		//tY[1] = forceIDXBnd(intLocY + 1, sy1i),
		//tZ[1] = forceIDXBnd(intLocZ + 1, sz1i);
	/*
		vector<int> idxs(8); int cnt = 0;
		for (unsigned int z = 0; z < 2; ++z) {for (unsigned int y = 0; y < 2; ++y) {for (unsigned int x = 0; x < 2; ++x) { idxs[cnt++] = IX(tX[x], tY[y], tZ[z]);}}}
	*/	
		int idx000 = IX(forceIDXBnd(intLocX, sx1i, 0), forceIDXBnd(intLocY, sy1i, 0), forceIDXBnd(intLocZ, sz1i, 0)),
			idx111 = IX(forceIDXBnd(intLocX + 1, sx1i, 0), forceIDXBnd(intLocY + 1, sy1i, 0), forceIDXBnd(intLocZ + 1, sz1i, 0));

		//int idx000 = IX(forceIDXBnd(intLocX, sx2i, 1), forceIDXBnd(intLocY, sy2i, 1), forceIDXBnd(intLocZ, sz2i, 1)),
		//	idx111 = IX(forceIDXBnd(intLocX + 1, sx2i, 1), forceIDXBnd(intLocY + 1, sy2i, 1), forceIDXBnd(intLocZ + 1, sz2i, 1));

		//double valx = interpM1X * Vx[idx000] + interpX*Vx[idx111],
		//	   valy = interpM1Y * Vy[idx000] + interpY*Vy[idx111],
		//	   valz = interpM1Z * Vz[idx000] + interpZ*Vz[idx111];

		//int idx = IX(forceIDXBnd((int)((testLoc(0) - center(0)) / cellSz(0) + halfNmCellX), sx1i),
		//			 forceIDXBnd((int)((testLoc(1) - center(1)) / cellSz(1) + halfNmCellY), sy1i),
		//			 forceIDXBnd((int)((testLoc(2) - center(2)) / cellSz(2) + halfNmCellZ), sz1i));

		//return Eigen::Vector3d(Vx[idx], Vy[idx], Vz[idx]);
		return Eigen::Vector3d(interpM1X * Vx[idx000] + interpX*Vx[idx111], interpM1Y * Vy[idx000] + interpY*Vy[idx111], interpM1Z * Vz[idx000] + interpZ*Vz[idx111]);
	}//getVelAtCell

}//namespace particleSystem

