#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QTextEdit>
#include <QButtonGroup>
#include <QProcess>
#include <QPushButton>
class MainWindow : public QWidget
{
  Q_OBJECT
  public:
  MainWindow(QWidget* parent=0);
 
public slots:
  void postProcess();
  void check(const QString&);
  void vcheck();
  // void help(const QString&);
  void vmerge(); //vogmerge button clicked
  void vcut(); //vogcut button clicked
  void fvmAdapt(); //FVMAdapt
  void import();
  void showQuality(QString, QProcess::ExitStatus, QString);
  void processItem(int);
  void generateVar();
  void pb();
  void showText(const QString&);
signals:
 
private:
  bool waitForQualityFile;
  // QStringList helpInfo;
   QTextEdit* display;
  QButtonGroup* buttonGroup;
  
};

class MyPushButton: public QPushButton{
  Q_OBJECT
  public:
  MyPushButton(const QString& title, QWidget* parent=0);
signals:
  void showText(const QString&);
protected:
  void mouseMoveEvent(QMouseEvent* event);
};
#endif

