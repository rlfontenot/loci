#include <QApplication>

 #include "mainwindow.h"

 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
     MainWindow mainWin;
     mainWin.show();
     qApp->setStyleSheet("QPushButton,QAction,QToolBar,QMenuBar,QToolButton{color:darkgreen; font:bold,large}");
     
     return app.exec();
 }
