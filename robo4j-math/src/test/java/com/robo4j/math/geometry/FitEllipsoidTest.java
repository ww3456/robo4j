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

import org.junit.Test;

/**
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class FitEllipsoidTest {

    @Test
    public void translateToCenter(){
        Matrix4d matrix4d1 = new Matrix4d(
                1,2,0,0,
                0,1,0,0,
                0,0,1,0,
                0,0,2,1);
        Tuple3d tuple3d = new Tuple3d(2,2,2);

        FitElilipsoid fitElilipsoid = new FitElilipsoid();
        Matrix4d result = fitElilipsoid.translateToCenter(tuple3d, matrix4d1);

        System.out.println("Translate: " + result);
    }

    @Test
    public void findRadiiTest(){
        FitElilipsoid fitElilipsoid = new FitElilipsoid();

        double[] array = {1D,2D,3D};
        Tuple3d result = fitElilipsoid.findRadii(array);
        System.out.println("Radii: " + result);
    }


}
