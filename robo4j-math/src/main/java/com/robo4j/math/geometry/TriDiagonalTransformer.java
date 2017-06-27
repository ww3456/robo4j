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

        this.transform();

    }

    public double[] getMain(){
        return main;
    }

    public double[] getSecondary(){
        return secondary;
    }

    public Matrix3d getQ() {
        if(this.cachedQ == null) {
            this.cachedQ = this.getQT().transpose();
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

    private void transform() {
        int m = this.householderVectors.length;
        double[] z = new double[m];

        for(int k = 0; k < m - 1; ++k) {
            double[] hK = this.householderVectors[k];
            this.main[k] = hK[k];
            double xNormSqr = 0.0D;

            for(int j = k + 1; j < m; ++j) {
                double c = hK[j];
                xNormSqr += c * c;
            }

            double a = hK[k + 1] > 0.0D?-Math.sqrt(xNormSqr):Math.sqrt(xNormSqr);
            this.secondary[k] = a;
            if(a != 0.0D) {
                hK[k + 1] -= a;
                double beta = -1.0D / (a * hK[k + 1]);
                Arrays.fill(z, k + 1, m, 0.0D);

                for(int i = k + 1; i < m; ++i) {
                    double[] hI = this.householderVectors[i];
                    double hKI = hK[i];
                    double zI = hI[i] * hKI;

                    for(int j = i + 1; j < m; ++j) {
                        double hIJ = hI[j];
                        zI += hIJ * hK[j];
                        z[j] += hIJ * hKI;
                    }

                    z[i] = beta * (z[i] + zI);
                }

                double gamma = 0.0D;

                int i;
                for(i = k + 1; i < m; ++i) {
                    gamma += z[i] * hK[i];
                }

                gamma *= beta / 2.0D;

                for(i = k + 1; i < m; ++i) {
                    z[i] -= gamma * hK[i];
                }

                for(i = k + 1; i < m; ++i) {
                    double[] hI = this.householderVectors[i];

                    for(int j = i; j < m; ++j) {
                        hI[j] -= hK[i] * z[j] + z[i] * hK[j];
                    }
                }
            }
        }

        this.main[m - 1] = this.householderVectors[m - 1][m - 1];
    }


}
