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
public class VectorNd {

	private double[] data;

	public VectorNd(int size) {
		this.data = new double[size];
	}

	public VectorNd(double[] data){
	    this.data = data;
    }

	private void setValue(double val) {
		for (int i = 0; i < data.length; i++) {
			data[i] = val;
		}
	}

	public double[] gedDataRef(){
	    return data;
    }

    public double getEntry(int i){
		return data[i];
	}

}
