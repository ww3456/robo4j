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

package com.robo4j.socket.http.util;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.SocketChannel;

/**
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public final class SocketUtil {

    /**
     * Used for normal message
     *
     * @param channel
     * @param buffer
     * @return
     * @throws IOException
     */
    public static  int readBuffer(SocketChannel channel, ByteBuffer buffer) throws IOException {
        buffer.rewind();
        int numberRead = channel.read(buffer);
        int position = 0;
        int totalRead = numberRead;
        while (numberRead != 0 && position <= totalRead) {
            numberRead = channel.read(buffer);
            if(numberRead > 0){
                totalRead += numberRead;
            } else if(numberRead == -1){
                numberRead = 0;
            }
            position++;
        }
        return totalRead;
    }

    /**
     * used for messages contains Base64 messages
     *
     * @param channel
     * @param buffer
     * @param stopper
     * @return
     * @throws IOException
     */
    public static int readBuffer(SocketChannel channel, ByteBuffer buffer, int stopper) throws IOException {
        int numberRead = channel.read(buffer);
        int position = 0;
        int totalRead = numberRead;
        while (numberRead != stopper && position <= totalRead) {
            numberRead = channel.read(buffer);
            if(numberRead > 0){
                totalRead += numberRead;
            }
            position++;
        }
        return totalRead;
    }

    public static int writeBuffer(SocketChannel channel, ByteBuffer buffer) throws IOException {
        int numberWriten = channel.write(buffer);
        int totalWritten = numberWriten;

        while(numberWriten > 0 && buffer.hasRemaining()){
            numberWriten = channel.write(buffer);
            totalWritten += numberWriten;
        }
        return totalWritten;
    }
}
