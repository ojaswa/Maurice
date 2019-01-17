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
**           Date  : 14.09.2015                                           **
****************************************************************************/
attribute vec2 vVertex;
attribute vec2 vTexCoord;
uniform mat4 vProjection;
uniform vec2 vImageSize;
varying vec2 fTexCoord;

void main(void) {
        vec2 ndcs_vertex = 2.0*vVertex/vImageSize - vec2(1.0, 1.0);
        gl_Position = vProjection * vec4(ndcs_vertex, 0.0, 1.0);

        fTexCoord = vTexCoord;
        gl_FrontColor = gl_Color;
}
