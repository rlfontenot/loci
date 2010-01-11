#ifndef CUTDIALOG_H
#define CUTDIALOG_H

#include <QDialog>
#include "grid.h"
#include "pages.h"


class QLabel;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QCheckBox;
class QSlider;

class CutDialog: public QWidget
{
  Q_OBJECT

public:
  CutDialog(LoadInfo ldInfo, float size, QWidget *parent = 0);
  signals:
  void cutInfoChanged(cutplane_info&);
   void loadInfoChanged(const LoadInfo&);
  void cutPressed();
  
private slots:
void setInfo();
  void planeSelected(int);
  void reset();
  void cut();
  void updateVars(QString iter);
  
  
private:
  QSlider* xslider1;
  QSlider* yslider1;
  QSlider* zslider1;
  QSlider* xslider2;
  QSlider* yslider2;
  QSlider* zslider2;
  FloatEdit* xEditor1;
  FloatEdit* yEditor1;
  FloatEdit* zEditor1;
  FloatEdit* xEditor2;
  FloatEdit* yEditor2;
  FloatEdit* zEditor2;
  QGroupBox* rotateBox;

  QLabel *caseLabel;
  QComboBox *comboIter;
  QComboBox *comboVar;





  
  cutplane_info info;  // Data structure passed from main window
  LoadInfo ld_info;
 
};

#endif
