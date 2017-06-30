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

import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public class SampleDataGenerator {

    public List<Tuple3d> generatePoints(double a, double b,
                                        double c, double shiftx, double shifty, double shiftz,
                                        double noiseIntensity) {
        List<Tuple3d> points = new ArrayList<>();
        double[] x;
        double[] y;
        double[] z;

        int numPoints = 1000;

        x = new double[numPoints];
        y = new double[numPoints];
        z = new double[numPoints];
        Random r = new Random();

        for (int i = 0; i < numPoints; i++) {
            double s = Math.toRadians(r.nextInt(360));
            double t = Math.toRadians(r.nextInt(360));

            x[i] = a * Math.cos(s) * Math.cos(t);
            y[i] = b * Math.cos(s) * Math.sin(t);
            z[i] = c * Math.sin(s);
        }

        double angle = Math.toRadians((Math.PI / 6));

        double[] xt = new double[numPoints];
        double[] yt = new double[numPoints];

        for (int i = 0; i < numPoints; i++) {
            xt[i] = x[i] * Math.cos(angle) - y[i] * Math.sin(angle);
            yt[i] = x[i] * Math.sin(angle) + y[i] * Math.cos(angle);
        }

        for (int i = 0; i < numPoints; i++) {
            x[i] = xt[i] + shiftx;
            y[i] = yt[i] + shifty;
            z[i] = z[i] + shiftz;
        }

        for (int i = 0; i < numPoints; i++) {
            x[i] = x[i] + r.nextDouble() * noiseIntensity;
            y[i] = y[i] + r.nextDouble() * noiseIntensity;
            z[i] = z[i] + r.nextDouble() * noiseIntensity;
        }

        Tuple3d tsp;

        for (int i = 0; i < numPoints; i++) {
            tsp = new Tuple3d(x[i], y[i], z[i]);
            points.add(tsp);
        }

        return points;
    }
}