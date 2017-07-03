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

	public double[] evals;

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

		Tuple9d v = solveInputs(ellipsoidInput);

		// Form the algebraic form of the ellipsoid.
		Matrix4d a = formAlgebraicMatrix(v);

		// Find the center of the ellipsoid.
		center = findCenter(a);

		// Translate the algebraic form of the ellipsoid to the center.
		Matrix4d translatedCenter = translateToCenter(center, a);

		// Generate a submatrix of r.
		Matrix3d tmpSubMatrix = translatedCenter.getSubMatrix3d0();
		double divisor = -translatedCenter.m44;
		Matrix3d subMatrixDiv = (Matrix3d) tmpSubMatrix.divideByNumber(divisor);

		EigenDecomposition3 ed = new EigenDecomposition3(subMatrixDiv);
		evals = ed.getRealEigenvalues();
		evecs0 = ed.getEigenvector(0);
		evecs1 = ed.getEigenvector(1);
		evecs2 = ed.getEigenvector(2);

		// Find the radii of the ellipsoid.
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
	public Tuple3d findRadii(double[] evals) {
		double[] radii = new double[evals.length];

		// radii[i] = sqrt(1/eval[i]);
		for (int i = 0; i < evals.length; i++) {
			radii[i] =Math.sqrt(1 / evals[i]);
		}

		return new Tuple3d(radii);
	}

	private Tuple3d findCenter(Matrix4d a) {
		// submatrix m(3,3)
		Matrix3d subA = a.getSubMatrix3d0();
		subA.multiplyByFactor(-1.0);

		// create sub vector
		Tuple3d subV = a.getRowVector3();
		// inv (dtd)
		DecompositionSolver solver = new SingularDecompositionSolver(subA).getSolver();

		Matrix3d subAi = (Matrix3d) solver.getInverse();

		// subAi(3x3)*subVector(3)' = result
		return subAi.operateMultiplyByVector3(subV);
	}

	public Matrix4d translateToCenter(Tuple3d center, Matrix4d a) {
		a.fitData();
		// Form the corresponding translation matrix.
		Matrix4d t = Matrix4d.createIdentity();

		Tuple3d centerMatrix = new Tuple3d();

		// centerMatrix.setRowVector(0, center);
		centerMatrix.set(center);

		t.setSubmatrixL3F0(centerMatrix);
		t.fitData();

		// Translate to the center.
		Matrix4d tmp1 = t.multiply(a);
		Matrix4d transposeT = t.transpose();
		transposeT.fitData();

		Matrix4d r = tmp1.multiply(transposeT);

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
	private Tuple9d solveInputs(List<Tuple3d> input) {

		// init design matrix 9x9
		Matrix9d d = new Matrix9d();

		for (int i = 0; i < d.getRows(); i++) {
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
		//ones is vector of entry elements N~1000

		// Multiply: d' * ones.mapAddToSelf(1)
		Matrix9d transposeMatrix9d = d.transpose();

		// Matrix(9x9) * IdentVector(9)' = dtOnes
		Tuple9d ones9 = new Tuple9d(1,1,1,1,1,1,1,1,1);
		Tuple9d dtOnes = transposeMatrix9d.operateMultiplyByVector9(ones9);

		// Find ( d' * d )^-1
		DecompositionSolver solver = new SingularDecompositionSolver(dtd).getSolver();

		Matrix9d dtdi = (Matrix9d) solver.getInverse();

		// v = (( d' * d )^-1) * ( d' * IdentityVectory(9));
		return dtdi.operateMultiplyByVector9(dtOnes);
	}

	/**
	 * Create a matrix in the algebraic form of the polynomial Ax^2 + By^2 +
	 * Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz = 1.
	 *
	 * @param v
	 *            the vector polynomial.
	 * @return the matrix of the algebraic form of the polynomial.
	 */
	private Matrix4d formAlgebraicMatrix(Tuple9d v) {
		// a =
		// [ Ax^2 2Dxy 2Exz 2Gx ]
		// [ 2Dxy By^2 2Fyz 2Hy ]
		// [ 2Exz 2Fyz Cz^2 2Iz ]
		// [ 2Gx 2Hy 2Iz -1 ] ]
		Matrix4d a = new Matrix4d(v.x1, v.x4, v.x5, v.x7, v.x4, v.x2, v.x6, v.x8, v.x5, v.x6, v.x3, v.x9, v.x7, v.x8, v.x9, -1);

		//FIXME : remove
		a.fitData();

		return a;
	}

}
