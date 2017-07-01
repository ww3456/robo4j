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

import com.robo4j.math.RoboOutOfRangeException;

/**
 * A three dimensional matrix.
 * 
 * @author Marcus Hirt (@hirt)
 * @author Miroslav Wengner (@miragemiko)
 */
public class Matrix3d implements Matrix {

	public static final int DIMENSION = 3;
	private double[][] data = new double[DIMENSION][DIMENSION];

	public double m11;
	public double m12;
	public double m13;
	public double m21;
	public double m22;
	public double m23;
	public double m31;
	public double m32;
	public double m33;

	public Matrix3d(double m11, double m12, double m13, double m21, double m22, double m23, double m31, double m32,
			double m33) {
		this.m11 = m11;
		this.m12 = m12;
		this.m13 = m13;
		this.m21 = m21;
		this.m22 = m22;
		this.m23 = m23;
		this.m31 = m31;
		this.m32 = m32;
		this.m33 = m33;
		data[0][0] = m11;
		data[0][1] = m12;
		data[0][2] = m13;

		data[1][0] = m21;
		data[1][1] = m22;
		data[1][2] = m23;

		data[2][0] = m31;
		data[2][1] = m32;
		data[2][2] = m33;
	}

	public Matrix3d(double[] matrix) {
		if (matrix.length != 9) {
			throw new IllegalArgumentException("Array argument for Matrix3f must be 9 elements long");
		}
		m11 = matrix[0];
		m12 = matrix[1];
		m13 = matrix[2];
		m21 = matrix[3];
		m22 = matrix[4];
		m23 = matrix[5];
		m31 = matrix[6];
		m32 = matrix[7];
		m33 = matrix[8];
	}

	public Matrix3d(double[][] matrix) {
		m11 = matrix[0][0];
		m12 = matrix[0][1];
		m13 = matrix[0][2];
		m21 = matrix[1][0];
		m22 = matrix[1][1];
		m23 = matrix[1][2];
		m31 = matrix[2][0];
		m32 = matrix[2][1];
		m33 = matrix[2][2];
		fitData();
	}

	@Override
	public double[][] getData() {
		data[0][0] = m11;
		data[0][1] = m12;
		data[0][2] = m13;
		data[1][0] = m21;
		data[1][1] = m22;
		data[1][2] = m23;
		data[2][0] = m31;
		data[2][1] = m32;
		data[2][2] = m33;

		return data;
	}

	@Override
	public int getDimension() {
		return DIMENSION;
	}

	/**
	 * add data to internal array structure
	 */
	@Override
	public void fitData() {
		data[0][0] = m11;
		data[0][1] = m12;
		data[0][2] = m13;

		data[1][0] = m21;
		data[1][1] = m22;
		data[1][2] = m23;

		data[2][0] = m31;
		data[2][1] = m32;
		data[2][2] = m33;

	}

	/**
	 * internal array structure  to values
	 */
	@Override
	public void adjustValues() {
		m11 = data[0][0];
		m12 = data[0][1];
		m13 = data[0][2];

		m21 = data[1][0];
		m22 = data[1][1];
		m23 = data[1][2];

		m31 = data[2][0];
		m32 = data[2][1];
		m33 = data[2][2];

	}

	public void setElement(int row, int column, double val) {
		data[row][column] = val;
	}

	/**
	 * Transforms the tuple by multiplying with this matrix.
	 * 
	 * @param tuple
	 *            the tuple to multiply with this matrix.
	 */
	public void transform(Tuple3d tuple) {
		tuple.set(m11 * tuple.x + m12 * tuple.y + m13 * tuple.z, m21 * tuple.x + m22 * tuple.y + m23 * tuple.z,
				m31 * tuple.x + m32 * tuple.y + m33 * tuple.z);
	}

	/**
	 * Like transform, but creating a new tuple without changing the old one.
	 * 
	 * @param tuple
	 *            the tuple to multiply with.
	 * @return the result from multiplying this matrix with the tuple.
	 */
	public Tuple3d multiply(Tuple3d tuple) {
		double x = m11 * tuple.x + m12 * tuple.y + m13 * tuple.z;
		double y = m21 * tuple.x + m22 * tuple.y + m23 * tuple.z;
		double z = m31 * tuple.x + m32 * tuple.y + m33 * tuple.z;
		return new Tuple3d(x, y, z);
	}

	@Override
	public Matrix multiply(Matrix m){
		double[][] tmp = new double[DIMENSION][DIMENSION];
		for(int i = 0; i < DIMENSION; i++) {         // rows from current matrix
			for(int j = 0; j < DIMENSION; j++) {     // columns matrix m
				for(int k = 0; k < DIMENSION; k++) { // columns from current matrix
					tmp[i][j] += data[i][k] * ((Matrix3d)m).getElement(k,j);
				}
			}
		}
		return new Matrix3d(tmp[0][0],tmp[0][1],tmp[0][2],
				tmp[1][0],tmp[1][1],tmp[1][2], tmp[2][0],tmp[2][1],tmp[2][2]);
	}

	public double getElement(int row, int column){
		return data[row][column];
	}

	/**
	 * Transposes the matrix.
	 */
	@Override
	public Matrix transpose() {
		double tmp = m12;
		m12 = m21;
		m21 = tmp;
		tmp = m13;
		m13 = m31;
		m31 = tmp;
		tmp = m23;
		m23 = m32;
		m32 = tmp;
		return new Matrix3d(m11, m12, m13, m21, m22, m23, m31, m32, m33);
	}

	/**
	 * Creates an identity matrix.
	 */
	public static Matrix3d createIdentity() {
		return new Matrix3d(1, 0, 0, 0, 1, 0, 0, 0, 1);
	}

	public void multiplyByFactor(double factor) {
		m11 = factor * m11;
		m12 = factor * m12;
		m13 = factor * m13;
		m21 = factor * m21;
		m22 = factor * m22;
		m23 = factor * m23;
		m31 = factor * m31;
		m32 = factor * m32;
		m33 = factor * m33;
	}

	public VectorNd operate(double[] v) {
		double[][] tmpData = getData();
		if (v.length != DIMENSION) {
			throw new RoboOutOfRangeException("dimension mismatch");
		} else {
			double[] result = new double[DIMENSION];

			for (int row = 0; row < DIMENSION; ++row) {
				double[] dataRow = tmpData[row];
				double sum = 0.0D;

				for (int i = 0; i < DIMENSION; ++i) {
					sum += dataRow[i] * v[i];
				}

				result[row] = sum;
			}

			return new VectorNd(result);
		}
	}

	@Override
	public String toString() {
		return String.format("m11:%f, m12:%f, m13:%f, m21:%f, m22:%f, m23:%f, m31:%f, m32:%f, m33:%f", m11, m12, m13,
				m21, m22, m23, m31, m32, m33);
	}
}
