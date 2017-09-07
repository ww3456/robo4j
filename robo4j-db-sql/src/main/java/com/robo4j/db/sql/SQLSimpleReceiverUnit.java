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

package com.robo4j.db.sql;

import com.robo4j.core.RoboContext;
import com.robo4j.core.RoboUnit;
import com.robo4j.db.sql.dto.ERoboResponse;

/**
 * Simple receiver for ERobo Entities
 *
 * @author Marcus Hirt (@hirt)
 * @author Miro Wengner (@miragemiko)
 */
public class SQLSimpleReceiverUnit extends RoboUnit<ERoboResponse> {

    public SQLSimpleReceiverUnit(RoboContext context, String id) {
        super(ERoboResponse.class, context, id);
    }

    @Override
    public void onMessage(ERoboResponse message) {
        System.out.println(getClass().getSimpleName() + " response message: " + message);
    }
}
