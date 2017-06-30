/*
 * Copyright (c) 2014, 2017, Marcus Hirt, Miroslav Wengner
 *
 * Robo4J is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Robo4J is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Robo4J. If not, see <http://www.gnu.org/licenses/>.
 */

package com.robo4j.math.geometry;

import java.util.List;

/**
 *
 * Determines the center, radii, eigenvalues and eigenvectors of the ellipse
 * using an algorithm from Yury Petrov's Ellipsoid Fit MATLAB script.
 *
 * polynomial expression Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy +
 * 2Iz = 1
 *
 * @see http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
 *
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class FitElilipsoid {

	public Tuple3d center;
	public Tuple3d radii;
	public Tuple3d evecs0;
	public Tuple3d evecs1;
	public Tuple3d evecs2;

	private double[] evals;

	public FitElilipsoid() {
	}

	/**
	 * Fit points to the polynomial expression Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz
	 * + 2Fyz + 2Gx + 2Hy + 2Iz = 1 determine the center and radii of the fit
	 * ellipsoid.
	 * 
	 * @param points
	 *            input ellipsoid , then fit ellipsoid
	 */
	void fitEllipsoid(List<Tuple3d> ellipsoidInput) {

		VectorNd v = solveInputs(ellipsoidInput);

		// Form the algebraic form of the ellipsoid.
		Matrix4d a = formAlgebraicMatrix(v);

		// Find the center of the ellipsoid.
		center = findCenter(a);

		// Translate the algebraic form of the ellipsoid to the center.
		Matrix4d r = translateToCenter(center, a);

		// Generate a submatrix of r.
		Matrix3d subr = r.getSubMatrix3d0();

		// subr[i][j] = subr[i][j] / -r[3][3]).
		r.fitData(); // set internal data value storage
		double divr = -r.getData()[3][3];
		for (int i = 0; i < subr.getDimension(); i++) {
			for (int j = 0; j < subr.getDimension(); j++) {
				subr.setElement(i, j, subr.getData()[i][j] / divr);
			}
		}
		// adjust data to m_xy values
		subr.adjustValues();

		EigenDecomposition ed = new EigenDecomposition(subr, 0);
		evals = ed.getRealEigenvalues();
		evecs0 = ed.getEigenvector(0);
		evecs1 = ed.getEigenvector(1);
		evecs2 = ed.getEigenvector(2);

		// Find the radii of the ellipsoid.
		//TODO : bug is here
		radii = findRadii(evals);
	}

	// Private Methods

	/**
	 * Find the radii of the ellipsoid in ascending order.
	 * 
	 * @param evals
	 *            the eigenvalues of the ellipsoid.
	 * @return the radii of the ellipsoid.
	 */
	private Tuple3d findRadii(double[] evals) {
		double[] radii = new double[evals.length];

		// radii[i] = sqrt(1/eval[i]);
		for (int i = 0; i < evals.length; i++) {
			double number = 1d / Math.abs(evals[i]);
			double tmp = Math.sqrt(number);
			radii[i] = tmp;
		}

		return new Tuple3d(radii);
	}

	private Tuple3d findCenter(Matrix4d a) {
		// submatrix m(3,3)
		Matrix3d subA = a.getSubMatrix3d0();
		subA.multiplyByFactor(-1.0);

		Tuple3d subV = a.getRowVector3();
		// inv (dtd)
		//TODO: bug is here
		DecompositionSolver solver = new SingularDecompositionSolver(subA).getSolver();

		Matrix3d subAi = (Matrix3d) solver.getInverse();

		VectorNd resultVector = subAi.operate(subV.getRefData());
		return new Tuple3d(resultVector.gedDataRef());
	}

	private Matrix4d translateToCenter(Tuple3d center, Matrix4d a) {
		// Form the corresponding translation matrix.
		Matrix4d t = Matrix4d.createIdentity();

		Tuple3d centerMatrix = new Tuple3d();

		// centerMatrix.setRowVector(0, center);
		centerMatrix.set(center);

		t.setSubmatrixL3F0(centerMatrix);

		// Translate to the center.
		Matrix4d r = t.multiply(a).multiply(t.transpose());

		return r;
	}

	/**
	 * Solve the polynomial expression: Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz
	 * + 2Gx + 2Hy + 2Iz by input ellipsoid points.
	 *
	 * @param input
	 *            represent points in 3D
	 * @return
	 */
	private VectorNd solveInputs(List<Tuple3d> input) {

		final int pointNumbers = input.size();

		// init design matrix 9x9
		Matrix9d d = new Matrix9d();

		for (int i = 0; i < d.getDimension(); i++) {
			double xx = Math.pow(input.get(i).x, 2);
			double yy = Math.pow(input.get(i).y, 2);
			double zz = Math.pow(input.get(i).z, 2);
			double xy = 2 * (input.get(i).x * input.get(i).y);
			double xz = 2 * (input.get(i).x * input.get(i).z);
			double yz = 2 * (input.get(i).y * input.get(i).z);
			double x = 2 * input.get(i).x;
			double y = 2 * input.get(i).y;
			double z = 2 * input.get(i).z;

			d.setElement(i, 0, xx);
			d.setElement(i, 1, yy);
			d.setElement(i, 2, zz);
			d.setElement(i, 3, xy);
			d.setElement(i, 4, xz);
			d.setElement(i, 5, yz);
			d.setElement(i, 6, x);
			d.setElement(i, 7, y);
			d.setElement(i, 8, z);
		}

		// solve the normal system of equations
		// v = (( d' * d )^-1) * ( d' * ones.mapAddToSelf(1));

		Matrix9d dtdTranspose = d.transpose();
		Matrix9d dtd = (Matrix9d)dtdTranspose.multiply(d);
		VectorNd ones = new VectorNd(pointNumbers);
		ones.mapAddToSelf(1);

		// Multiply: d' * ones.mapAddToSelf(1)
		Matrix9d transposeMatrix9d = d.transpose();
		VectorNd dtOnes = transposeMatrix9d.operate(ones.gedDataRef(), ones.gedDataRef().length);


		// Find ( d' * d )^-1
		DecompositionSolver solver = new SingularDecompositionSolver(dtd).getSolver();

		 Matrix9d dtdi = (Matrix9d) solver.getInverse();

		// v = (( d' * d )^-1) * ( d' * ones.mapAddToSelf(1));

		return dtdi.operate(dtOnes.gedDataRef());
	}

	/**
	 * Create a matrix in the algebraic form of the polynomial Ax^2 + By^2 +
	 * Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz = 1.
	 *
	 * @param v
	 *            the vector polynomial.
	 * @return the matrix of the algebraic form of the polynomial.
	 */
	private Matrix4d formAlgebraicMatrix(VectorNd v) {
		// a =
		// [ Ax^2 2Dxy 2Exz 2Gx ]
		// [ 2Dxy By^2 2Fyz 2Hy ]
		// [ 2Exz 2Fyz Cz^2 2Iz ]
		// [ 2Gx 2Hy 2Iz -1 ] ]
		Matrix4d a = new Matrix4d();

//		a.setEntry(0, 0, v.getEntry(0));
//		a.setEntry(0, 1, v.getEntry(3));
//		a.setEntry(0, 2, v.getEntry(4));
//		a.setEntry(0, 3, v.getEntry(6));

//		a.setEntry(1, 0, v.getEntry(3));
//		a.setEntry(1, 1, v.getEntry(1));
//		a.setEntry(1, 2, v.getEntry(5));
//		a.setEntry(1, 3, v.getEntry(7));

//		a.setEntry(2, 0, v.getEntry(4));
//		a.setEntry(2, 1, v.getEntry(5));
//		a.setEntry(2, 2, v.getEntry(2));
//		a.setEntry(2, 3, v.getEntry(8));
//		a.setEntry(3, 0, v.getEntry(6));
//		a.setEntry(3, 1, v.getEntry(7));
//		a.setEntry(3, 2, v.getEntry(8));
//		a.setEntry(3, 3, -1);

		a.m11 = v.getEntry(0);
		a.m12 = v.getEntry(3);
		a.m13 = v.getEntry(4);
		a.m14 = v.getEntry(6);

		a.m21 = v.getEntry(3);
		a.m22 = v.getEntry(1);
		a.m23 = v.getEntry(5);
		a.m24 = v.getEntry(7);

		a.m31 = v.getEntry(4);
		a.m32 = v.getEntry(5);
		a.m33 = v.getEntry(2);
		a.m34 = v.getEntry(8);

		a.m41 = v.getEntry(6);
		a.m42 = v.getEntry(7);
		a.m43 = v.getEntry(8);
		a.m44 = -1;


		a.setElement(0,0, a.m11);
		a.setElement(0,1, a.m12);
		a.setElement(0,2, a.m13);
		a.setElement(0,3, a.m14);

		a.setElement(1,0, a.m21);
		a.setElement(1,1, a.m22);
		a.setElement(1,2, a.m23);
		a.setElement(1,3, a.m24);

		a.setElement(2,0, a.m31);
		a.setElement(2,1, a.m32);
		a.setElement(2,2, a.m33);
		a.setElement(2,3, a.m34);

		a.setElement(3,0, a.m41);
		a.setElement(3,1, a.m42);
		a.setElement(3,2, a.m43);
		a.setElement(3,3, a.m44);

		return a;
	}

}
