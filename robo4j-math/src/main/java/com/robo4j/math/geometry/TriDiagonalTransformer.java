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

/**
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class TriDiagonalTransformer {

    private int DIMENSION = 3;
    private final double[][] householderVectors;
    private final double[] main;
    private final double[] secondary;
    private Matrix3d matrix3d;
    private Matrix3d cachedQ;
    private Matrix3d cachedQt;
    private Matrix3d cachedT;

    public TriDiagonalTransformer(Matrix3d matrix3d) {
        this.matrix3d = matrix3d;
        householderVectors = matrix3d.getData();
        this.main = new double[DIMENSION];
        this.secondary = new double[DIMENSION - 1];

        this.cachedQ = null;
        this.cachedQt = null;
        this.cachedT = null;

        transform();

    }

    public double[] getMain(){
        return main;
    }

    public double[] getSecondary(){
        return secondary;
    }

    public Matrix3d getQ() {
        if(this.cachedQ == null) {
            this.cachedQ = (Matrix3d) this.getQT().transpose();
        }

        return this.cachedQ;
    }

    public Matrix3d getQT() {
        if(this.cachedQt == null) {
            int m = this.householderVectors.length;
            double[][] qta = new double[m][m];

            for(int k = m - 1; k >= 1; --k) {
                double[] hK = this.householderVectors[k - 1];
                qta[k][k] = 1.0D;
                if(hK[k] != 0.0D) {
                    double inv = 1.0D / (this.secondary[k - 1] * hK[k]);
                    double beta = 1.0D / this.secondary[k - 1];
                    qta[k][k] = 1.0D + beta * hK[k];

                    int j;
                    for(j = k + 1; j < m; ++j) {
                        qta[k][j] = beta * hK[j];
                    }

                    for(j = k + 1; j < m; ++j) {
                        beta = 0.0D;

                        int i;
                        for(i = k + 1; i < m; ++i) {
                            beta += qta[j][i] * hK[i];
                        }

                        beta *= inv;
                        qta[j][k] = beta * hK[k];

                        for(i = k + 1; i < m; ++i) {
                            qta[j][i] += beta * hK[i];
                        }
                    }
                }
            }

            qta[0][0] = 1.0D;
            this.cachedQt = new Matrix3d(qta[0][0], qta[0][1], qta[0][2],
                    qta[1][0], qta[1][1], qta[1][2], qta[2][0], qta[2][1], qta[2][2]);
        }

        return this.cachedQt;
    }

    /**
     * Transform original matrix to tridiagonal form.
     * <p>Transformation is done using Householder transforms.</p>
     */
    private void transform() {
        final int m = householderVectors.length;
        final double[] z = new double[m];
        for (int k = 0; k < m - 1; k++) {

            //zero-out a row and a column simultaneously
            final double[] hK = householderVectors[k];
            main[k] = hK[k];
            double xNormSqr = 0;
            for (int j = k + 1; j < m; ++j) {
                final double c = hK[j];
                xNormSqr += c * c;
            }
            final double a = (hK[k + 1] > 0) ? - Math.sqrt(xNormSqr) : Math.sqrt(xNormSqr);
            secondary[k] = a;
            if (a != 0.0) {
                // apply Householder transform from left and right simultaneously

                hK[k + 1] -= a;
                final double beta = -1 / (a * hK[k + 1]);

                // compute a = beta A v, where v is the Householder vector
                // this loop is written in such a way
                //   1) only the upper triangular part of the matrix is accessed
                //   2) access is cache-friendly for a matrix stored in rows
                Arrays.fill(z, k + 1, m, 0);
                for (int i = k + 1; i < m; ++i) {
                    final double[] hI = householderVectors[i];
                    final double hKI = hK[i];
                    double zI = hI[i] * hKI;
                    for (int j = i + 1; j < m; ++j) {
                        final double hIJ = hI[j];
                        zI   += hIJ * hK[j];
                        z[j] += hIJ * hKI;
                    }
                    z[i] = beta * (z[i] + zI);
                }

                // compute gamma = beta vT z / 2
                double gamma = 0;
                for (int i = k + 1; i < m; ++i) {
                    gamma += z[i] * hK[i];
                }
                gamma *= beta / 2;

                // compute z = z - gamma v
                for (int i = k + 1; i < m; ++i) {
                    z[i] -= gamma * hK[i];
                }

                // update matrix: A = A - v zT - z vT
                // only the upper triangular part of the matrix is updated
                for (int i = k + 1; i < m; ++i) {
                    final double[] hI = householderVectors[i];
                    for (int j = i; j < m; ++j) {
                        hI[j] -= hK[i] * z[j] + z[i] * hK[j];
                    }
                }
            }
        }
        main[m - 1] = householderVectors[m - 1][m - 1];
    }


}
