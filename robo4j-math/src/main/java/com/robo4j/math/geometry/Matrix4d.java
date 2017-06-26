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

/**
 * A four dimensional matrix. (The things we do to save the array size field
 * bytes... ;))
 * 
 * @author Marcus Hirt (@hirt)
 * @author Miroslav Wengner (@miragemiko)
 */
public class Matrix4d {
	private static final int DIMENSION = 4;
	private double[][] data = new double[DIMENSION][DIMENSION];
	public double m11;
	public double m12;
	public double m13;
	public double m14;
	public double m21;
	public double m22;
	public double m23;
	public double m24;
	public double m31;
	public double m32;
	public double m33;
	public double m34;
	public double m41;
	public double m42;
	public double m43;
	public double m44;

	public Matrix4d() {
	}

	public Matrix4d(double m11, double m12, double m13, double m14, double m21, double m22, double m23, double m24, double m31, double m32, double m33,
					double m34, double m41, double m42, double m43, double m44) {
		this.m11 = m11;
		this.m12 = m12;
		this.m13 = m13;
		this.m14 = m14;
		this.m21 = m21;
		this.m22 = m22;
		this.m23 = m23;
		this.m24 = m24;
		this.m31 = m31;
		this.m32 = m32;
		this.m33 = m33;
		this.m34 = m34;
		this.m41 = m41;
		this.m42 = m42;
		this.m43 = m43;
		this.m44 = m44;
	}

	public Matrix4d(double[] matrix) {
		if (matrix.length != 16) {
			throw new IllegalArgumentException("Array argument for Matrix3f must be 9 elements long");
		}
		m11 = matrix[0];
		m12 = matrix[1];
		m13 = matrix[2];
		m14 = matrix[3];
		m21 = matrix[4];
		m22 = matrix[5];
		m23 = matrix[6];
		m24 = matrix[7];
		m31 = matrix[8];
		m32 = matrix[9];
		m33 = matrix[10];
		m34 = matrix[11];
		m41 = matrix[12];
		m42 = matrix[13];
		m43 = matrix[14];
		m44 = matrix[15];
	}

	public void fitData(){
		data[0][0] = m11;
		data[0][1] = m12;
		data[0][2] = m13;
		data[0][3] = m14;

		data[1][0] = m21;
		data[1][1] = m22;
		data[1][2] = m23;
		data[1][3] = m24;

		data[2][0] = m31;
		data[2][1] = m32;
		data[2][2] = m33;
		data[2][3] = m34;

		data[3][0] = m41;
		data[3][1] = m42;
		data[3][2] = m43;
		data[3][3] = m44;
	}

	public double[][] getData(){
		return data;
	}

	public void adjustValues(){
		m11 = data[0][0];
		m12 = data[0][1];
		m13 = data[0][2];
		m14 = data[0][3];

		m21 = data[1][0];
		m22 = data[1][1];
		m23 = data[1][2];
		m24 = data[1][3];

		m31 = data[2][0];
		m32 = data[2][1];
		m33 = data[2][2];
		m34 = data[2][3];

		m41 = data[3][0];
		m42 = data[3][1];
		m43 = data[3][2];
		m44 = data[3][3];
	}

	public void setElement(int row, int column, double val){
		data[row][column] = val;
	}

	/**
	 * Transforms the tuple by multiplying with this matrix.
	 * 
	 * @param tuple
	 *            the tuple to multiply with this matrix.
	 */
	public void transform(Tuple4d tuple) {
		tuple.set(m11 * tuple.x + m12 * tuple.y + m13 * tuple.z + m14 * tuple.t,
				m21 * tuple.x + m22 * tuple.y + m23 * tuple.z + m24 * tuple.t,
				m31 * tuple.x + m32 * tuple.y + m33 * tuple.z + m34 * tuple.t,
				m41 * tuple.x + m42 * tuple.y + m43 * tuple.z + m44 * tuple.t);
	}

	/**
	 * Transposes the matrix.
	 */
	public Matrix4d transpose() {
		double tmp = m12;
		m12 = m21;
		m21 = tmp;
		tmp = m13;
		m13 = m31;
		m31 = tmp;
		tmp = m14;
		m14 = m41;
		m41 = tmp;
		tmp = m23;
		m23 = m32;
		m32 = tmp;
		tmp = m34;
		m34 = m43;
		m43 = tmp;
		return this;
	}

	public void setSubmatrixL3F0(Tuple3d tuple3d){
		m31 = tuple3d.x;
		m32 = tuple3d.y;
		m33 = tuple3d.z;
	}

	/**
	 * Like transform, but creating a new tuple without changing the old one.
	 * 
	 * @param tuple
	 *            the tuple to multiply with.
	 * @return the result from multiplying this matrix with the tuple.
	 */
	public Tuple4d multiply(Tuple4d tuple) {
		double x = m11 * tuple.x + m12 * tuple.y + m13 * tuple.z + m14 * tuple.t;
		double y = m21 * tuple.x + m22 * tuple.y + m23 * tuple.z + m24 * tuple.t;
		double z = m31 * tuple.x + m32 * tuple.y + m33 * tuple.z + m34 * tuple.t;
		double t = m41 * tuple.x + m42 * tuple.y + m43 * tuple.z + m44 * tuple.t;
		return new Tuple4d(x, y, z, t);
	}

	public Matrix4d multiply(Matrix4d m){

		Matrix4d out = new Matrix4d();
		this.fitData();
		for(int row = 0; row < DIMENSION; ++row) {
			for(int col = 0; col < DIMENSION; ++col) {
				double sum = 0.0D;

				for(int i = 0; i < DIMENSION; ++i) {
					sum += data[row][i] * m.getData()[i][col];
				}

				out.setElement(row, col, sum);
			}
		}

		out.adjustValues();
		return out;
	}

	public Tuple3d getRowVector3(){
		return new Tuple3d(m41, 42, 43);
	}

	/**
	 * Creates an identity matrix.
	 */
	public static Matrix4d createIdentity() {
		return new Matrix4d(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
	}

	// submatrix starts from 0,0
	public Matrix3d getSubMatrix3d0() {
		return new Matrix3d(m11, m12, m13, m21, m22, m23, m31, m32, m33);
	}

	@Override
	public String toString() {
		return String.format(
				"m11:%f, m12:%f, m13:%f, m14:%f, m21:%f, m22:%f, m23:%f, m24:%f, m31:%f, m32:%f, m33:%f, m34:%f, m41:%f, 42:%f, m43:%f, m44:%f",
				m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44);
	}
}
