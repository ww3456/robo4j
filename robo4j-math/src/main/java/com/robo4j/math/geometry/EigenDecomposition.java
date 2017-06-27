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

    private void transformToTridiagonal(Matrix3d matrix){
        transformer = new TriDiagonalTransformer(matrix);
        main = transformer.getMain();
        secondary = transformer.getSecondary();
    }

    private void findEigenVectors(double[][] householderMatrix){
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
        e[n - 1] = 0.0D;
        double maxAbsoluteValue = 0.0D;

        int j;
        for(j = 0; j < n; ++j) {
            if(Math.abs(this.realEigenvalues[j]) > maxAbsoluteValue) {
                maxAbsoluteValue = Math.abs(this.realEigenvalues[j]);
            }

            if(Math.abs(e[j]) > maxAbsoluteValue) {
                maxAbsoluteValue = Math.abs(e[j]);
            }
        }

        if(maxAbsoluteValue != 0.0D) {
            for(j = 0; j < n; ++j) {
                if(Math.abs(this.realEigenvalues[j]) <= 1.1102230246251565E-16D * maxAbsoluteValue) {
                    this.realEigenvalues[j] = 0.0D;
                }

                if(Math.abs(e[j]) <= 1.1102230246251565E-16D * maxAbsoluteValue) {
                    e[j] = 0.0D;
                }
            }
        }

        int i;
        int m;
        for(j = 0; j < n; ++j) {
            i = 0;

            do {
                double q;
                for(m = j; m < n - 1; ++m) {
                    q = Math.abs(this.realEigenvalues[m]) + Math.abs(this.realEigenvalues[m + 1]);
                    if(Math.abs(e[m]) + q == q) {
                        break;
                    }
                }

                if(m != j) {

                    ++i;
                    q = (this.realEigenvalues[j + 1] - this.realEigenvalues[j]) / (2.0D * e[j]);
                    double t = Math.sqrt(1.0D + q * q);
                    if(q < 0.0D) {
                        q = this.realEigenvalues[m] - this.realEigenvalues[j] + e[j] / (q - t);
                    } else {
                        q = this.realEigenvalues[m] - this.realEigenvalues[j] + e[j] / (q + t);
                    }

                    double u = 0.0D;
                    double s = 1.0D;
                    double c = 1.0D;

                    int i1;
                    for(i1 = m - 1; i1 >= j; --i1) {
                        double p = s * e[i1];
                        double h = c * e[i1];
                        if(Math.abs(p) >= Math.abs(q)) {
                            c = q / p;
                            t = Math.sqrt(c * c + 1.0D);
                            e[i1 + 1] = p * t;
                            s = 1.0D / t;
                            c *= s;
                        } else {
                            s = p / q;
                            t = Math.sqrt(s * s + 1.0D);
                            e[i1 + 1] = q * t;
                            c = 1.0D / t;
                            s *= c;
                        }

                        if(e[i1 + 1] == 0.0D) {
                            this.realEigenvalues[i1 + 1] -= u;
                            e[m] = 0.0D;
                            break;
                        }

                        q = this.realEigenvalues[i1 + 1] - u;
                        t = (this.realEigenvalues[i1] - q) * s + 2.0D * c * h;
                        u = s * t;
                        this.realEigenvalues[i1 + 1] = q + u;
                        q = c * t - h;

                        for(int ia = 0; ia < n; ++ia) {
                            p = z[ia][i1 + 1];
                            z[ia][i1 + 1] = s * z[ia][i1] + c * p;
                            z[ia][i1] = c * z[ia][i1] - s * p;
                        }
                    }

                    if(t != 0.0D || i < j) {
                        this.realEigenvalues[j] -= u;
                        e[j] = q;
                        e[m] = 0.0D;
                    }
                }
            } while(m != j);
        }

        for(j = 0; j < n; ++j) {
            i = j;
            double p = this.realEigenvalues[j];

            int j1;
            for(j1 = j + 1; j1 < n; ++j1) {
                if(this.realEigenvalues[j1] > p) {
                    i = j1;
                    p = this.realEigenvalues[j1];
                }
            }

            if(i != j1) {
                this.realEigenvalues[i] = this.realEigenvalues[j1];
                this.realEigenvalues[j1] = p;

                for(j1 = 0; j1 < n; ++j) {
                    p = z[j1][j1];
                    z[j1][j1] = z[j1][i];
                    z[j1][i] = p;
                }
            }
        }

        maxAbsoluteValue = 0.0D;

        for(j = 0; j < n; ++j) {
            if(Math.abs(this.realEigenvalues[j]) > maxAbsoluteValue) {
                maxAbsoluteValue = Math.abs(this.realEigenvalues[j]);
            }
        }

        if(maxAbsoluteValue != 0.0D) {
            for(j = 0; j < n; ++j) {
                if(Math.abs(this.realEigenvalues[j]) < 1.1102230246251565E-16D * maxAbsoluteValue) {
                    this.realEigenvalues[j] = 0.0D;
                }
            }
        }

        double[] tmp = new double[n];

        for(i = 0; i < n; ++i) {
            for(m = 0; m < n; ++m) {
                tmp[m] = z[m][i];
            }

            this.eigenvectors[i] = tmp.clone();
        }
    }
}
