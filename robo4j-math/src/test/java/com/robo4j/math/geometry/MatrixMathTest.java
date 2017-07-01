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
public class MatrixMathTest {


    @Test
    public void multiply1(){

        Matrix4d matrix4d1 = new Matrix4d(1,2,0,0,0,1,0,0,0,0,1,0,0,0,2,1);
        Matrix4d matrix4d2 = new Matrix4d(1,2,0,0,0,1,0,0,0,0,1,0,0,0,2,1);


        Matrix4d result = matrix4d1.multiply(matrix4d2).multiply(matrix4d2);

        System.out.println("Result : " + result);
    }

    @Test
    public void transpose4d1(){
        Matrix4d matrix4d1 = new Matrix4d(
                //@formatter:off
                1,2,0,0,
                0,1,0,0,
                0,0,1,0,
                0,3,2,1
                //@formatter:on
        );

        Matrix4d transposeMatrix = matrix4d1.transpose();

        System.out.println("Transpose: " + transposeMatrix);

    }

    @Test
    public void transpose9d1(){

        double[][] values = new double[][]{
                //@formatter:off
                {1, 0, 0, 0, 0, 0, 0,0, 9},
                {0, 1, 0, 0, 0, 0, 0,0, 0},
                {0, 0, 1, 0, 0, 0, 0,0, 0},
                {0, 0, 0, 1, 0, 0, 0,0, 0},
                {0, 0, 0, 0, 1, 0, 0,0, 0},
                {0, 0, 0, 0, 0, 1, 0,0, 0},
                {0, 0, 0, 0, 0, 0, 1,0, 0},
                {0, 0, 0, 0, 0, 0, 0,1, 3},
                {0, 0, 0, 0, 0, 0, 0,0, 1},
                //@formatter:on
        };
        Matrix9d matrix9d = new Matrix9d(values);
        Matrix9d transpose =  matrix9d.transpose();
        System.out.println("Transpose:" + transpose);
    }
}
