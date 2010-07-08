#ifndef CUTDIALOG_H
#define CUTDIALOG_H

#include <QMainWindow>
#include "grid.h"
#include "pages.h"
#include "glviewer.h"

class QLabel;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QCheckBox;
class QSlider;

class CutDialog: public QMainWindow
{
  Q_OBJECT

public:
  CutDialog(LoadInfo, QWidget *parent = 0);
  
  signals:

  
 
  // void extremePressed();
 

private slots:
void setInfo();
  void planeSelected(int);
  void reset();
  void cut();
  void load();
  void showExtremeNodes(int);
  void updateVars(QString iter);
  
  
private:
  QSlider* xslider1;
  QSlider* yslider1;
  QSlider* zslider1;
  QSlider* xslider2;
  QSlider* yslider2;
  QSlider* zslider2;
  DoubleEdit* xEditor1;
  DoubleEdit* yEditor1;
  DoubleEdit* zEditor1;
  DoubleEdit* xEditor2;
  DoubleEdit* yEditor2;
  DoubleEdit* zEditor2;
  QGroupBox* rotateBox;

  QLabel *caseLabel;
  QComboBox *comboIter;
  QComboBox *comboVar;
  QSpinBox* extrSpinBox;
  
  QGroupBox* toolbox;

  
  GLViewer *viewer; 
  cutplane_info info;  // Data structure passed from main window
  LoadInfo ld_info;
  double size;
  void createToolBar();
 
};

#endif
