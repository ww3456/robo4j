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

/**
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class EllipsoidFitTest {

    private static final List<Tuple3d> CONTROL_ELLIPSOID_POINTS = new ArrayList<>();
    static {
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(1.218234969106912, 1.5622932396950777, 2.908127141636472));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(3.406914196934893, 2.0905555158425777, 1.9887962320370343));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.2647839110418686, 2.2647839110418686, 3.153654761400918));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));
        CONTROL_ELLIPSOID_POINTS.add(new Tuple3d(2.0624005905759786, 1.96337625716086, 3.200989974576277));

    }



    public static void main(String[] args) {
        System.out.println("ellipsoid fit test");



        // Fit the ellipsoid points to a polynomial
        FitElilipsoid ellipsoidFit = new FitElilipsoid();
        ellipsoidFit.fitEllipsoid(CONTROL_ELLIPSOID_POINTS);

        System.out.println("CENTER: " + ellipsoidFit.center);
        System.out.println("RADII: " + ellipsoidFit.radii);
    }

}
