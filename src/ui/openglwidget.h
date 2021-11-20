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
**           Date  : 28.10.2015                                           **
****************************************************************************/
#ifndef OPENGLWIDGET_H
#define OPENGLWIDGET_H

#include<QOpenGLWidget>
#include<QOpenGLFunctions>
#include<QElapsedTimer>

#include <shape.h>

#define GLM_FORCE_RADIANS
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#define HANDLE_SELECT_THRESHOLD 10.0 // In pixel units
#define NODE_SELECT_SIZE 5.0

#define M_PI 3.14159265f

#define printOpenGLError() printOglError(__FILE__, __LINE__)
#define degreeToRadians(X) ((X)*M_PI/180.0f)

class OpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT

public:
    OpenGLWidget(QWidget *parent);
    ~OpenGLWidget();
    Shape* getShape();
    void deformVertices();
protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    GLuint createProgram(const char *vshader_filename, const char* fshader_filename);
    char* getShaderCode(const char* filename);
    void printLog(GLuint object);
    GLuint createShader(const char* filename, GLenum type);
    int printOglError(const char *file, int line);
    void setupProjectionTransformation();
    void createShape(); //Gets called (using signal-slot mechanism) when an image is read
    void setupLineHandleDisplay();
    glm::vec2 screenToImage(int x, int y);
    void createCage();


    void mouseMoveEvent(QMouseEvent* event);
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent( QKeyEvent *k );
    void timerEvent(QTimerEvent *event);

private:
    int m_screenWidth, m_screenHeight;

    // Shape rendering variables
    GLuint m_shapeProgram;
    GLint m_shapeVertexAttrib, m_shapeTextureAttrib;
    GLint m_shapeProjectionUniform, m_shapeUniformTexture, m_shapeImageSizeUniform;
    GLint m_shapeUniformWireframe;
    GLuint m_shapeVBO, m_shapeIBO, m_texVBO;
    GLuint m_shapeTexture;
    Shape *m_shape;
    bool m_isShapeTriangulated;
    bool m_drawWireframe;
    bool m_modeAddLineHandle;
    bool m_modeAddPolylineHandle;
    glm::mat4 m_projection;
    GLfloat *m_shapeVertices;

    //Line drawing variables
    GLuint m_lineProgram;
    GLint m_lineVertexAttrib;
    GLint m_lineProjectionUniform, m_lineImageSizeUniform;
    std::vector<float> m_linePoints; //last two values point to currently moving mouse point
    bool m_isLineDrawing;
    int m_selectedNodeIndex;
    int m_selectedHandleIndex;
    bool m_hideHandles;
    bool m_showCage;
    GLuint m_cageVBO;
    bool m_recordVideo;
    long m_videoFrame;
    QString m_videoPath;
    int m_timerId;
    QElapsedTimer m_elapsedTimer;

public slots:
    void imageLoaded();
    void toggleDisplayMesh(bool);
    void toggleDisplayCage(bool);
    void toggleAddLineHandle(bool);
    void toggleAddPolylineHandle(bool);
    void toggleUseEuclideanDistanceWeights(bool);
    void resetDeformation();
    void deleteAllHandles();
    void takeSnapshot();
    void takeSnapshot(QString filename);
    void hideHandles(bool);
    void dilateCage(int);
    void simplifyCage(int);
};

#endif // OPENGLWIDGET_H
