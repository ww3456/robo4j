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
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class EigenDecomposition {

    /**
     * Smallest positive number such that {@code 1 - EPSILON} is not
     * numerically equal to 1: {@value}.
     */
    public static final double EPSILON = 0x1.0p-53;

    private static final int DIMENSION = 3;
    private Matrix3d matrix3d;
    private double splitTolerance;
    private TriDiagonalTransformer transformer;
    private double[] main;
    private double[] secondary;
    private double[] realEigenvalues;
    private double[] imagEigenvalues;
    private double[][] eigenvectors = new double[DIMENSION][];

    public EigenDecomposition(Matrix3d matrix3d, double splitTolerance) {

        this.matrix3d = matrix3d;
        this.splitTolerance = splitTolerance;
        transformToTridiagonal(matrix3d);
        findEigenVectors(transformer.getQ().getData());
    }

    public double[] getRealEigenvalues() {
        return this.realEigenvalues.clone();
    }

    public Tuple3d getEigenvector(int i) {
        return new Tuple3d(eigenvectors[i]);
    }

    void transformToTridiagonal(Matrix3d matrix){
        transformer = new TriDiagonalTransformer(matrix);
        main = transformer.getMain();
        secondary = transformer.getSecondary();
    }

    public void findEigenVectors(double[][] householderMatrix) {

        double[][] z = householderMatrix.clone();
        int n = this.main.length;
        this.realEigenvalues = new double[n];
        this.imagEigenvalues = new double[n];
        double[] e = new double[n];

        for(int i = 0; i < n - 1; ++i) {
            this.realEigenvalues[i] = this.main[i];
            e[i] = this.secondary[i];
        }

        this.realEigenvalues[n - 1] = this.main[n - 1];
        e[n - 1] = 0;

        // Determine the largest main and secondary value in absolute term.
        double maxAbsoluteValue = 0;

        for (int i = 0; i < n; i++) {
            if (Math.abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = Math.abs(realEigenvalues[i]);
            }
            if (Math.abs(e[i]) > maxAbsoluteValue) {
                maxAbsoluteValue = Math.abs(e[i]);
            }
        }
        //fixed

        // Make null any main and secondary value too small to be significant
        if (maxAbsoluteValue != 0) {
            for (int i=0; i < n; i++) {
                if (Math.abs(realEigenvalues[i]) <= EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0;
                }
                if (Math.abs(e[i]) <= EPSILON * maxAbsoluteValue) {
                    e[i]=0;
                }
            }
        }

        //start
        for (int j = 0; j < n; j++) {
            int m;
            do {
                for (m = j; m < n - 1; m++) {
                    double delta = Math.abs(realEigenvalues[m]) +
                            Math.abs(realEigenvalues[m + 1]);
                    if (Math.abs(e[m]) + delta == delta) {
                        break;
                    }
                }
                if (m != j) {
                    double q = (realEigenvalues[j + 1] - realEigenvalues[j]) / (2 * e[j]);
                    double t = Math.sqrt(1 + q * q);
                    if (q < 0.0) {
                        q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q - t);
                    } else {
                        q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q + t);
                    }
                    double u = 0.0;
                    double s = 1.0;
                    double c = 1.0;
                    int i;
                    for (i = m - 1; i >= j; i--) {
                        double p = s * e[i];
                        double h = c * e[i];
                        if (Math.abs(p) >= Math.abs(q)) {
                            c = q / p;
                            t = Math.sqrt(c * c + 1.0);
                            e[i + 1] = p * t;
                            s = 1.0 / t;
                            c = c * s;
                        } else {
                            s = p / q;
                            t = Math.sqrt(s * s + 1.0);
                            e[i + 1] = q * t;
                            c = 1.0 / t;
                            s = s * c;
                        }
                        if (e[i + 1] == 0.0) {
                            realEigenvalues[i + 1] -= u;
                            e[m] = 0.0;
                            break;
                        }
                        q = realEigenvalues[i + 1] - u;
                        t = (realEigenvalues[i] - q) * s + 2.0 * c * h;
                        u = s * t;
                        realEigenvalues[i + 1] = q + u;
                        q = c * t - h;
                        for (int ia = 0; ia < n; ia++) {
                            p = z[ia][i + 1];
                            z[ia][i + 1] = s * z[ia][i] + c * p;
                            z[ia][i] = c * z[ia][i] - s * p;
                        }
                    }
                    if (t == 0.0 && i >= j) {
                        continue;
                    }
                    realEigenvalues[j] -= u;
                    e[j] = q;
                    e[m] = 0.0;
                }
            } while (m != j);
        }
        //stop


        //Sort the eigen values (and vectors) in increase order
        for (int i = 0; i < n; i++) {
            int k = i;
            double p = realEigenvalues[i];
            for (int j = i + 1; j < n; j++) {
                if (realEigenvalues[j] > p) {
                    k = j;
                    p = realEigenvalues[j];
                }
            }
            if (k != i) {
                realEigenvalues[k] = realEigenvalues[i];
                realEigenvalues[i] = p;
                for (int j = 0; j < n; j++) {
                    p = z[j][i];
                    z[j][i] = z[j][k];
                    z[j][k] = p;
                }
            }
        }

        // Determine the largest eigen value in absolute term.
        maxAbsoluteValue = 0;
        for (int i = 0; i < n; i++) {
            if (Math.abs(realEigenvalues[i]) > maxAbsoluteValue) {
                maxAbsoluteValue=Math.abs(realEigenvalues[i]);
            }
        }

        // Make null any eigen value too small to be significant
        if (maxAbsoluteValue!=0.0) {
            for (int i=0; i < n; i++) {
                if (Math.abs(realEigenvalues[i]) < EPSILON * maxAbsoluteValue) {
                    realEigenvalues[i] = 0;
                }
            }
        }


        final double[] tmp = new double[n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                tmp[j] = z[j][i];
            }
            eigenvectors[i] = tmp.clone();
        }

    }
}
