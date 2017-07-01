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

import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class EigenDecompositionTest {

    @Test
    public void eigenVectorsTest1(){

        Matrix3d matrix3d = new Matrix3d(1,0,0,0,1,0,0,0,1);
        EigenDecomposition ed = new EigenDecomposition(matrix3d, 0);

        System.out.println("matrix: " + matrix3d);
        System.out.println("eigenValues: " + Stream.of(ed.getRealEigenvalues()).collect(Collectors.toList()));
        System.out.println("eigen0: " + ed.getEigenvector(0));
        System.out.println("eigen1: " + ed.getEigenvector(1));
        System.out.println("eigen2: " + ed.getEigenvector(2));


    }

    @Test
    public void eigenVectorsTest2(){

        Matrix3d matrix3d = new Matrix3d(
                1,7,3,
                7,4,-5,
                3,-5,6);
        EigenDecomposition ed = new EigenDecomposition(matrix3d, 0);

        System.out.println("matrix: " + matrix3d);
        System.out.println("eigenValues: " + Stream.of(ed.getRealEigenvalues()).collect(Collectors.toList()));
        System.out.println("eigen0: " + ed.getEigenvector(0));
        System.out.println("eigen1: " + ed.getEigenvector(1));
        System.out.println("eigen2: " + ed.getEigenvector(2));


    }

    @Test
    public void findEigenVectorsTest(){


        Matrix3d matrix3d = new Matrix3d(
                1,0,0,
                0,1,0,
                0,0,1);
        EigenDecomposition ed = new EigenDecomposition(matrix3d, 0);




    }

}
