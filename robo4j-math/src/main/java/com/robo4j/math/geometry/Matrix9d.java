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

import java.util.Arrays;

import com.robo4j.math.RoboOutOfRangeException;

/**
 *
 * @see FitElilipsoid
 *
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class Matrix9d implements Matrix {
	private static final int DIMENSION = 9;
	private double[][] data = new double[DIMENSION][DIMENSION];

	public Matrix9d() {
	}

	public Matrix9d(double[][] data) {
		this.data = data;
	}

	@Override
	public int getRows() {
		return DIMENSION;
	}

	@Override
	public int getColumns() {
		return DIMENSION;
	}

	@Override
	public double[][] getData() {
		return data;
	}

	public void setElement(int row, int column, double value) {
		if (checkRowColumn(row, column)) {
			data[row][column] = value;
		} else {
			throw new RoboOutOfRangeException(String.format("wrong numbers: row: %d, column: %d", row, column));
		}
	}

	public double getElement(int row, int column) {
		return data[row][column];
	}

	@Override
	public Matrix9d transpose() {
		double[][] temp = new double[DIMENSION][data[0].length];
		for (int i = 0; i < DIMENSION; i++)
			for (int j = 0; j < data[0].length; j++)
				temp[j][i] = data[i][j];
		return new Matrix9d(temp);
	}

	@Override
	public Matrix multiply(Matrix m) {

		double[][] tmp = new double[DIMENSION][DIMENSION];
		for (int i = 0; i < DIMENSION; i++) { // rows from current matrix
			for (int j = 0; j < DIMENSION; j++) { // columns matrix m
				for (int k = 0; k < DIMENSION; k++) { // columns from current
														// matrix
					tmp[i][j] += data[i][k] * ((Matrix9d) m).getElement(k, j);
				}
			}
		}
		return new Matrix9d(tmp);
	}

    /**
     * Result =Matrix(9x9) * OnesVector(9)'
     *
     * create vector
     */
    public Tuple9d operateMultiplyByVector9(Tuple9d v) {

            double r1 = data[0][0]*v.x1 + data[0][1]*v.x2 + data[0][2]*v.x3 + data[0][3]*v.x4 + data[0][4]*v.x5 + data[0][5]*v.x6 + data[0][6]*v.x7 + data[0][7]*v.x8 + data[0][8]*v.x9;
            double r2 = data[1][0]*v.x1 + data[1][1]*v.x2 + data[1][2]*v.x3 + data[1][3]*v.x4 + data[1][4]*v.x5 + data[1][5]*v.x6 + data[1][6]*v.x7 + data[1][7]*v.x8 + data[1][8]*v.x9;
            double r3 = data[2][0]*v.x1 + data[2][1]*v.x2 + data[2][2]*v.x3 + data[2][3]*v.x4 + data[2][4]*v.x5 + data[2][5]*v.x6 + data[2][6]*v.x7 + data[2][7]*v.x8 + data[2][8]*v.x9;
            double r4 = data[3][0]*v.x1 + data[3][1]*v.x2 + data[3][2]*v.x3 + data[3][3]*v.x4 + data[3][4]*v.x5 + data[3][5]*v.x6 + data[3][6]*v.x7 + data[3][7]*v.x8 + data[3][8]*v.x9;
            double r5 = data[4][0]*v.x1 + data[4][1]*v.x2 + data[4][2]*v.x3 + data[4][3]*v.x4 + data[4][4]*v.x5 + data[4][5]*v.x6 + data[4][6]*v.x7 + data[4][7]*v.x8 + data[4][8]*v.x9;
            double r6 = data[5][0]*v.x1 + data[5][1]*v.x2 + data[5][2]*v.x3 + data[5][3]*v.x4 + data[5][4]*v.x5 + data[5][5]*v.x6 + data[5][6]*v.x7 + data[5][7]*v.x8 + data[5][8]*v.x9;
            double r7 = data[6][0]*v.x1 + data[6][1]*v.x2 + data[6][2]*v.x3 + data[6][3]*v.x4 + data[6][4]*v.x5 + data[6][5]*v.x6 + data[6][6]*v.x7 + data[6][7]*v.x8 + data[6][8]*v.x9;
            double r8 = data[7][0]*v.x1 + data[7][1]*v.x2 + data[7][2]*v.x3 + data[7][3]*v.x4 + data[7][4]*v.x5 + data[7][5]*v.x6 + data[7][6]*v.x7 + data[7][7]*v.x8 + data[7][8]*v.x9;
            double r9 = data[8][0]*v.x1 + data[8][1]*v.x2 + data[8][2]*v.x3 + data[8][3]*v.x4 + data[8][4]*v.x5 + data[8][5]*v.x6 + data[8][6]*v.x7 + data[8][7]*v.x8 + data[8][8]*v.x9;

            return new Tuple9d(r1, r2, r3, r4, r5, r6, r7, r8, r9);
    }

	@Override
	public double getValue(int row, int column) {
        throw new RuntimeException("not implemented");
	}

	@Override
	public void fitData() {
		throw new RuntimeException("not implemented");
	}

	@Override
	public void adjustValues() {
		throw new RuntimeException("not implemented");
	}

	@Override
	public String toString() {
		return "Matrix9d{" + "data=" + Arrays.toString(data) + '}';
	}

	// Private Methods
	private boolean checkRowColumn(int row, int column) {
		return 0 <= row && row < DIMENSION && 0 <= column && column < DIMENSION;
	}
}
