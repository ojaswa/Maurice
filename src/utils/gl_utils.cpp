/***************************************************************************
**                                                                        **
**  Maurice - a 2D deformable mapping/animation program                   **
**  Copyright (C) 2014-2016 Graphics Research Group, IIIT Delhi           **
**                                                                        **
**  This program is free software: you can redistribute it and/or modify  **
**  it under the terms of the GNU General Public License as published by  **
**  the Free Software Foundation, either version 3 of the License, or     **
**  (at your option) any later version.                                   **
**                                                                        **
**  This program is distributed in the hope that it will be useful,       **
**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
**  GNU General Public License for more details.                          **
**                                                                        **
**  You should have received a copy of the GNU General Public License     **
**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
**                                                                        **
****************************************************************************
**           Author: Ojaswa Sharma                                        **
**           E-mail: ojaswa@iiitd.ac.in                                   **
**           Date  : 03.09.2014                                           **
****************************************************************************/
#include "gl_utils.h"

#ifdef __APPLE__
#include <OpenGL/glu.h>
#else
#include <GL/glu.h>
#endif
#include <stdio.h>

int printOglError(const char *file, int line)
{
	GLenum glErr;
	int retCode = 0;
	glErr = glGetError();
	if(glErr != GL_NO_ERROR)
	{
		fprintf(stderr, "glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
		retCode = 1;
	}
	return retCode;
}
