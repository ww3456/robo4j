package com.robo4j.math.geometry;

/**
 * Created by miroslavwengner on 28.06.17.
 */
public class Matrix9rNrd {

    private static final int DIMENSION_ROWS = 9;
    private double[][] data = new double[DIMENSION_ROWS][];

    public Matrix9rNrd() {
    }

    public Matrix9rNrd(double[][] data){
        if(checkRowColumn(data.length, data[0].length)){
            this.data = data;
        }
    }

    public Matrix9d transpose(){
        double[][] temp = new double[DIMENSION_ROWS][];
        for (int i = 0; i < DIMENSION_ROWS; i++)
            for (int j = 0; j < data[0].length; j++)
                temp[j][i] = data[i][j];
        return new Matrix9d(temp);
    }


    //Private Methods
    private boolean checkRowColumn(int row, int column){
        return 0 <= row && row < DIMENSION_ROWS && 0 <= column;
    }
}

