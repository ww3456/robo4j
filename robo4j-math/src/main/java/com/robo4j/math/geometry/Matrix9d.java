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

import java.util.Arrays;

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

    public Matrix9d(double[][] data){
        this.data = data;
    }

    @Override
    public int getDimension(){
        return DIMENSION;
    }

    @Override
    public double[][] getData(){
        return data;
    }

    public void setElement(int row, int column, double value){
        if(checkRowColumn(row, column)){
            data[row][column] = value;
        } else {
            throw new RoboOutOfRangeException(String.format("wrong numbers: row: %d, column: %d", row, column));
        }
    }

    public double getElement(int row, int column){
        return data[row][column];
    }

    @Override
    public Matrix9d transpose(){
        double[][] temp = new double[DIMENSION][data[0].length];
        for (int i = 0; i < DIMENSION; i++)
            for (int j = 0; j < data[0].length; j++)
                temp[j][i] = data[i][j];
        return new Matrix9d(temp);
    }

    @Override
    public Matrix multiply(Matrix m){

        double[][] tmp = new double[DIMENSION][DIMENSION];
        for(int i = 0; i < DIMENSION; i++) {         // rows from current matrix
            for(int j = 0; j < DIMENSION; j++) {     // columns matrix m
                for(int k = 0; k < DIMENSION; k++) { // columns from current matrix
                    tmp[i][j] += data[i][k] * ((Matrix9d)m).getElement(k,j);
                }
            }
        }
        return new Matrix9d(tmp);
    }

    public VectorNd operate(double[] v){
        if(v.length != DIMENSION) {
            throw new RoboOutOfRangeException("dimension mismatch");
        } else {
            double[] result = new double[DIMENSION];

            for(int row = 0; row < DIMENSION; ++row) {
                double[] dataRow = data[row];
                double sum = 0.0D;

                for(int i = 0; i < DIMENSION; ++i) {
                    sum += dataRow[i] * v[i];
                }

                result[row] = sum;
            }

            return new VectorNd(result);
        }
    }

    public VectorNd operate(double[] v, int columns){
        if(v.length != columns) {
            throw new RoboOutOfRangeException("dimension mismatch");
        } else {
            double[] result = new double[DIMENSION];

            for(int row = 0; row < DIMENSION; ++row) {
                double[] dataRow = data[row];
                double sum = 0.0D;

                for(int i = 0; i < DIMENSION; ++i) {
                    sum += dataRow[i] * v[i];
                }

                result[row] = sum;
            }

            return new VectorNd(result);
        }
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
        return "Matrix9d{" +
                "data=" + Arrays.toString(data) +
                '}';
    }

    //Private Methods
    private boolean checkRowColumn(int row, int column){
        return 0 <= row && row < DIMENSION && 0 <= column && column < DIMENSION;
    }
}
