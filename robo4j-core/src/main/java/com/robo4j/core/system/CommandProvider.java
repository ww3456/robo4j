/*
 * Copyright (C)  2016. Miroslav Wengner, Marcus Hirt
 * This CommandProvider.java  is part of robo4j.
 *
 *  robo4j is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  robo4j is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with robo4j .  If not, see <http://www.gnu.org/licenses/>.
 */

package com.robo4j.core.system;

import com.robo4j.commons.command.CommandTargetEnum;
import com.robo4j.commons.command.GenericCommand;
import com.robo4j.commons.command.PlatformCommandEnum;
import com.robo4j.commons.enums.RoboHardwareEnum;

/**
 * RoboUnit doesn't included the definition of the process method. Reason
 * because process can be done by Agent without human input
 *
 * @author Miroslav Wengner (@miragemiko)
 * @since 10.06.2016
 */
public interface CommandProvider {

	boolean process(GenericCommand<? extends RoboHardwareEnum<?, CommandTargetEnum>> command);
}
