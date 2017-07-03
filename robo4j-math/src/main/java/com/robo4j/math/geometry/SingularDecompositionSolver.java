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
 * inspired by SingularValueDecomposition (Apache Math3)
 *
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class SingularDecompositionSolver {

    private final double[] singularValues;
    private final int n;
    private final double tol;
    private final Matrix cachedU;
    private Matrix cachedUt;
    private final Matrix cachedV;


    public SingularDecompositionSolver(Matrix matrix) {

        double[][] A = matrix.getData();
        this.n = matrix.getRows();

        this.singularValues = new double[this.n];
        double[][] U = new double[this.n][this.n];
        double[][] V = new double[this.n][this.n];
        double[] e = new double[this.n];
        double[] work = new double[this.n];
        int nct = min(this.n - 1, this.n);
        int nrt = max(0, this.n - 2);

        int p;
        int pp;
        int iter;
        for (p = 0; p < max(nct, nrt); ++p) {
            if (p < nct) {
                this.singularValues[p] = 0.0D;

                for (pp = p; pp < this.n; ++pp) {
                    this.singularValues[p] = hypot(this.singularValues[p], A[pp][p]);
                }

                if (this.singularValues[p] != 0.0D) {
                    if (A[p][p] < 0.0D) {
                        this.singularValues[p] = -this.singularValues[p];
                    }

                    for (pp = p; pp < this.n; ++pp) {
                        A[pp][p] /= this.singularValues[p];
                    }

                    ++A[p][p];
                }

                this.singularValues[p] = -this.singularValues[p];
            }

            double t;
            int i;
            for (pp = p + 1; pp < this.n; ++pp) {
                if (p < nct && this.singularValues[p] != 0.0D) {
                    t = 0.0D;

                    for (i = p; i < this.n; ++i) {
                        t += A[i][p] * A[i][pp];
                    }

                    t = -t / A[p][p];

                    for (i = p; i < this.n; ++i) {
                        A[i][pp] += t * A[i][p];
                    }
                }

                e[pp] = A[p][pp];
            }

            if (p < nct) {
                for (pp = p; pp < this.n; ++pp) {
                    U[pp][p] = A[pp][p];
                }
            }

            if (p < nrt) {
                e[p] = 0.0D;

                for (pp = p + 1; pp < this.n; ++pp) {
                    e[p] = hypot(e[p], e[pp]);
                }

                if (e[p] != 0.0D) {
                    if (e[p + 1] < 0.0D) {
                        e[p] = -e[p];
                    }

                    for (pp = p + 1; pp < this.n; ++pp) {
                        e[pp] /= e[p];
                    }

                    ++e[p + 1];
                }

                e[p] = -e[p];
                if (p + 1 < this.n && e[p] != 0.0D) {
                    for (pp = p + 1; pp < this.n; ++pp) {
                        work[pp] = 0.0D;
                    }

                    for (pp = p + 1; pp < this.n; ++pp) {
                        for (iter = p + 1; iter < this.n; ++iter) {
                            work[iter] += e[pp] * A[iter][pp];
                        }
                    }

                    for (pp = p + 1; pp < this.n; ++pp) {
                        t = -e[pp] / e[p + 1];

                        for (i = p + 1; i < this.n; ++i) {
                            A[i][pp] += t * work[i];
                        }
                    }
                }

                for (pp = p + 1; pp < this.n; ++pp) {
                    V[pp][p] = e[pp];
                }
            }
        }

        p = this.n;
        if (nct < this.n) {
            this.singularValues[nct] = A[nct][nct];
        }

        if (nrt + 1 < p) {
            e[nrt] = A[nrt][p - 1];
        }

        e[p - 1] = 0.0D;

        for (pp = nct; pp < this.n; ++pp) {
            for (iter = 0; iter < this.n; ++iter) {
                U[iter][pp] = 0.0D;
            }

            U[pp][pp] = 1.0D;
        }

        double t;
        int i;
        for (pp = nct - 1; pp >= 0; --pp) {
            if (this.singularValues[pp] != 0.0D) {
                for (iter = pp + 1; iter < this.n; ++iter) {
                    t = 0.0D;

                    for (i = pp; i < this.n; ++i) {
                        t += U[i][pp] * U[i][iter];
                    }

                    t = -t / U[pp][pp];

                    for (i = pp; i < this.n; ++i) {
                        U[i][iter] += t * U[i][pp];
                    }
                }

                for (iter = pp; iter < this.n; ++iter) {
                    U[iter][pp] = -U[iter][pp];
                }

                ++U[pp][pp];

                for (iter = 0; iter < pp - 1; ++iter) {
                    U[iter][pp] = 0.0D;
                }
            } else {
                for (iter = 0; iter < this.n; ++iter) {
                    U[iter][pp] = 0.0D;
                }

                U[pp][pp] = 1.0D;
            }
        }

        for (pp = this.n - 1; pp >= 0; --pp) {
            if (pp < nrt && e[pp] != 0.0D) {
                for (iter = pp + 1; iter < this.n; ++iter) {
                    t = 0.0D;

                    for (i = pp + 1; i < this.n; ++i) {
                        t += V[i][pp] * V[i][iter];
                    }

                    t = -t / V[pp + 1][pp];

                    for (i = pp + 1; i < this.n; ++i) {
                        V[i][iter] += t * V[i][pp];
                    }
                }
            }

            for (iter = 0; iter < this.n; ++iter) {
                V[iter][pp] = 0.0D;
            }

            V[pp][pp] = 1.0D;
        }

        pp = p - 1;
        iter = 0;

        while (true) {
            label315:
            while (p > 0) {
                int k;
                double t1;
                for (k = p - 2; k >= 0; --k) {
                    t1 = 1.6033346880071782E-291D + 2.220446049250313E-16D * (abs(this.singularValues[k]) + abs(this.singularValues[k + 1]));
                    if (abs(e[k]) <= t1) {
                        e[k] = 0.0D;
                        break;
                    }
                }

                byte kase;
                if (k == p - 2) {
                    kase = 4;
                } else {
                    for (i = p - 1; i >= k && i != k; --i) {
                        double t2 = (i != p ? abs(e[i]) : 0.0D) + (i != k + 1 ? abs(e[i - 1]) : 0.0D);
                        if (abs(this.singularValues[i]) <= 1.6033346880071782E-291D + 2.220446049250313E-16D * t2) {
                            this.singularValues[i] = 0.0D;
                            break;
                        }
                    }

                    if (i == k) {
                        kase = 3;
                    } else if (i == p - 1) {
                        kase = 1;
                    } else {
                        kase = 2;
                        k = i;
                    }
                }

                ++k;
                double t3;
                double cs;
                double sn;
                int i1;
                switch (kase) {
                    case 1:
                        t = e[p - 2];
                        e[p - 2] = 0.0D;
                        i1 = p - 2;

                        while (true) {
                            if (i1 < k) {
                                continue label315;
                            }

                            t3 = hypot(this.singularValues[i1], t);
                            cs = this.singularValues[i1] / t3;
                            sn = t3 / t3;
                            this.singularValues[i1] = t3;
                            if (i1 != k) {
                                t3 = -sn * e[i1 - 1];
                                e[i1 - 1] = cs * e[i1 - 1];
                            }

                            for (i1 = 0; i1 < this.n; ++i1) {
                                t3 = cs * V[i1][i1] + sn * V[i1][p - 1];
                                V[i1][p - 1] = -sn * V[i1][i1] + cs * V[i1][p - 1];
                                V[i1][i1] = t3;
                            }

                            --i1;
                        }
                    case 2:
                        t = e[k - 1];
                        e[k - 1] = 0.0D;
                        i = k;

                        while (true) {
                            if (i >= p) {
                                continue label315;
                            }

                            t3 = hypot(this.singularValues[i], t);
                            cs = this.singularValues[i] / t3;
                            sn = t3 / t3;
                            this.singularValues[i] = t3;
                            t = -sn * e[i];
                            e[i] = cs * e[i];

                            for (i = 0; i < this.n; ++i) {
                                t = cs * U[i][i] + sn * U[i][k - 1];
                                U[i][k - 1] = -sn * U[i][i] + cs * U[i][k - 1];
                                U[i][i] = t;
                            }

                            ++i;
                        }
                    case 3:
                        t3 = max(abs(this.singularValues[p - 1]), abs(this.singularValues[p - 2]));
                        double scale = max(max(max(t3, abs(e[p - 2])), abs(this.singularValues[k])), abs(e[k]));
                        double sp = this.singularValues[p - 1] / scale;
                        double spm1 = this.singularValues[p - 2] / scale;
                        double epm1 = e[p - 2] / scale;
                        double sk = this.singularValues[k] / scale;
                        double ek = e[k] / scale;
                        double b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) / 2.0D;
                        double c = sp * epm1 * sp * epm1;
                        double shift = 0.0D;
                        if (b != 0.0D || c != 0.0D) {
                            shift = Math.sqrt(b * b + c);
                            if (b < 0.0D) {
                                shift = -shift;
                            }

                            shift = c / (b + shift);
                        }

                        double f = (sk + sp) * (sk - sp) + shift;
                        double g = sk * ek;

                        for (int j = k; j < p - 1; ++j) {
                            double t4 = hypot(f, g);
                            double cs1 = f / t4;
                            double sn1 = g / t4;
                            if (j != k) {
                                e[j - 1] = t4;
                            }

                            f = cs1 * this.singularValues[j] + sn1 * e[j];
                            e[j] = cs1 * e[j] - sn1 * this.singularValues[j];
                            g = sn1 * this.singularValues[j + 1];
                            this.singularValues[j + 1] = cs1 * this.singularValues[j + 1];

                            int i2;
                            for (i2 = 0; i2 < this.n; ++i2) {
                                t4 = cs1 * V[i2][j] + sn1 * V[i2][j + 1];
                                V[i2][j + 1] = -sn1 * V[i2][j] + cs1 * V[i2][j + 1];
                                V[i2][j] = t4;
                            }

                            t4 = hypot(f, g);
                            cs = f / t4;
                            sn = g / t4;
                            this.singularValues[j] = t4;
                            f = cs * e[j] + sn * this.singularValues[j + 1];
                            this.singularValues[j + 1] = -sn * e[j] + cs * this.singularValues[j + 1];
                            g = sn * e[j + 1];
                            e[j + 1] = cs * e[j + 1];
                            if (j < this.n - 1) {
                                for (i = 0; i < this.n; ++i) {
                                    t4 = cs * U[i][j] + sn * U[i][j + 1];
                                    U[i][j + 1] = -sn * U[i][j] + cs * U[i][j + 1];
                                    U[i][j] = t4;
                                }
                            }
                        }

                        e[p - 2] = f;
                        ++iter;
                        continue;
                }

                if (this.singularValues[k] <= 0.0D) {
                    this.singularValues[k] = this.singularValues[k] < 0.0D ? -this.singularValues[k] : 0.0D;

                    for (i = 0; i <= pp; ++i) {
                        V[i][k] = -V[i][k];
                    }
                }

                for (; k < pp && this.singularValues[k] < this.singularValues[k + 1]; ++k) {
                    t = this.singularValues[k];
                    this.singularValues[k] = this.singularValues[k + 1];
                    this.singularValues[k + 1] = t;
                    if (k < this.n - 1) {
                        for (i = 0; i < this.n; ++i) {
                            t = V[i][k + 1];
                            V[i][k + 1] = V[i][k];
                            V[i][k] = t;
                        }
                    }

                    if (k < this.n - 1) {
                        for (i = 0; i < this.n; ++i) {
                            t = U[i][k + 1];
                            U[i][k + 1] = U[i][k];
                            U[i][k] = t;
                        }
                    }
                }

                iter = 0;
                --p;
            }

            this.tol = max((double) this.n * this.singularValues[0] * 2.220446049250313E-16D, Math.sqrt(2.2250738585072014E-308D));

            switch (matrix.getRows()) {
                case 3:
                    this.cachedU = new Matrix3d(U);
                    this.cachedV = new Matrix3d(V);
                    return;
                case 9:
                    this.cachedU = new Matrix9d(U);
                    this.cachedV = new Matrix9d(V);
                    return;
                default:
                    throw new RuntimeException("not supported matrix");
            }

        }

    }

    public Matrix getUT() {
        if (this.cachedUt == null) {
            this.cachedUt = cachedU.transpose();
        }

        return this.cachedUt;
    }


    public int getRank() {
        int r = 0;

        for (int i = 0; i < this.singularValues.length; ++i) {
            if (this.singularValues[i] > this.tol) {
                ++r;
            }
        }

        return r;
    }

    public DecompositionSolver getSolver() {

        return new SingularDecompositionSolver.Solver(this.singularValues, this.getUT(), this.cachedV,
                this.getRank() == this.n, this.tol);
    }

    //FastMath
    private double hypot(double x, double y) {
        if (!Double.isInfinite(x) && !Double.isInfinite(y)) {
            if (!Double.isNaN(x) && !Double.isNaN(y)) {
                int expX = getExponent(x);
                int expY = getExponent(y);
                if (expX > expY + 27) {
                    return abs(x);
                } else if (expY > expX + 27) {
                    return abs(y);
                } else {
                    int middleExp = (expX + expY) / 2;
                    double scaledX = scalb(x, -middleExp);
                    double scaledY = scalb(y, -middleExp);
                    double scaledH = Math.sqrt(scaledX * scaledX + scaledY * scaledY);
                    return scalb(scaledH, middleExp);
                }
            } else {
                return 0.0D / 0.0;
            }
        } else {
            return 1.0D / 0.0;
        }
    }

    private static int getExponent(double d) {
        return (int) (Double.doubleToLongBits(d) >>> 52 & 2047L) - 1023;
    }

    private static double abs(double x) {
        return x < 0.0D ? -x : (x == 0.0D ? 0.0D : x);
    }

    private static double scalb(double d, int n) {
        if (n > -1023 && n < 1024) {
            return d * Double.longBitsToDouble((long) (n + 1023) << 52);
        } else if (!Double.isNaN(d) && !Double.isInfinite(d) && d != 0.0D) {
            if (n < -2098) {
                return d > 0.0D ? 0.0D : -0.0D;
            } else if (n > 2097) {
                return d > 0.0D ? 1.0D / 0.0 : -1.0D / 0.0;
            } else {
                long bits = Double.doubleToLongBits(d);
                long sign = bits & -9223372036854775808L;
                int exponent = (int) (bits >>> 52) & 2047;
                long mantissa = bits & 4503599627370495L;
                int scaledExponent = exponent + n;
                if (n < 0) {
                    if (scaledExponent > 0) {
                        return Double.longBitsToDouble(sign | (long) scaledExponent << 52 | mantissa);
                    } else if (scaledExponent > -53) {
                        mantissa |= 4503599627370496L;
                        long mostSignificantLostBit = mantissa & 1L << -scaledExponent;
                        mantissa >>>= 1 - scaledExponent;
                        if (mostSignificantLostBit != 0L) {
                            ++mantissa;
                        }

                        return Double.longBitsToDouble(sign | mantissa);
                    } else {
                        return sign == 0L ? 0.0D : -0.0D;
                    }
                } else if (exponent != 0) {
                    return scaledExponent < 2047 ? Double.longBitsToDouble(sign | (long) scaledExponent << 52 | mantissa) : (sign == 0L ? 1.0D / 0.0 : -1.0D / 0.0);
                } else {
                    while (mantissa >>> 52 != 1L) {
                        mantissa <<= 1;
                        --scaledExponent;
                    }

                    ++scaledExponent;
                    mantissa &= 4503599627370495L;
                    return scaledExponent < 2047 ? Double.longBitsToDouble(sign | (long) scaledExponent << 52 | mantissa) : (sign == 0L ? 1.0D / 0.0 : -1.0D / 0.0);
                }
            }
        } else {
            return d;
        }
    }

    public static double max(double a, double b) {
        if (a > b) {
            return a;
        } else if (a < b) {
            return b;
        } else if (a != b) {
            return 0.0D / 0.0;
        } else {
            long bits = Double.doubleToRawLongBits(a);
            return bits == -9223372036854775808L ? b : a;
        }
    }

    private static int min(int a, int b) {
        return a <= b ? a : b;
    }

    private static int max(int a, int b) {
        return a <= b ? b : a;
    }


    private static class Solver implements DecompositionSolver {
        private final Matrix pseudoInverse;
        private boolean nonSingular;

        private Solver(double[] singularValues, Matrix uT, Matrix v, boolean nonSingular, double tol) {
            double[][] suT = uT.getData();

            for (int i = 0; i < singularValues.length; ++i) {
                double a;
                if (singularValues[i] > tol) {
                    a = 1.0D / singularValues[i];
                } else {
                    a = 0.0D;
                }

                double[] suTi = suT[i];

                for (int j = 0; j < suTi.length; ++j) {
                    suTi[j] *= a;
                }
            }


            switch (v.getRows()) {
                case 3:
                    this.pseudoInverse = v.multiply(new Matrix3d(suT));
                    this.nonSingular = nonSingular;
                    break;
                case 9:
                    this.pseudoInverse = v.multiply(new Matrix9d(suT));
                    this.nonSingular = nonSingular;
                    break;
                default:
                    throw new RuntimeException("not supported matrix operation");
            }
        }


        public Matrix getInverse() {
            return this.pseudoInverse;
        }
    }

}
