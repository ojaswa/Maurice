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
**           Date  : 10.10.2015                                           **
****************************************************************************/
#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "cagedialog.h"
#include "openglwidget.h"
#include "shape.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    OpenGLWidget *openGLWidget;


private slots:
    void on_action_Open_Image_triggered();

    void on_actionAdd_Line_handles_triggered(bool checked);

    void on_actionAdd_Polyline_handle_triggered(bool checked);

    void on_actionModify_cage_triggered();

private:
    Ui::MainWindow *ui;
    CageDialog *cageDialog;
    Shape *m_shape;

signals:
    void imageLoaded();
};

#endif // MAINWINDOW_H
