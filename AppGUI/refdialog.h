#ifndef REFDIALOG_H
#define REFDIALOG_H

#include <QDialog>



class QLabel;
class QComboBox;
class QDoubleSpinBox;
class QGroupBox;
class QCheckBox;
class QSlider;
class QTextEdit;
class QButtonGroup;
class RefDialog: public QWidget
{
  Q_OBJECT

public:
  RefDialog(QString basefile,  QWidget *parent = 0);
  signals:
  void textChanged();
  
private slots:
void save();
  void runScript();
  void getfile(int);
  void balanceChanged(int i);
  void optChanged(int i);
  void modeChanged(int i);
  void setFold(double);
  void setLevels(int);
  void setText();
  void setTol(double);
private:
  QString basefile;
  QString tagfile;
  QString xmlfile;
  QString planfile;
  QString rplanfile;
  QString outfile;
  QString scriptfile;
  
  double tol;
  double fold;
  int levels;
  
  int mode;
  int balanceOpt;
  bool xmlOpt;
  
  QButtonGroup* fileButtonGroup;
  
  QTextEdit* display;
  
 
  
 
};

#endif
