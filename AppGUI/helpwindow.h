#ifndef HELPWINDOW_H
#define HELPWINDOW_H
#include <QWidget>
#include <QString>
class HelpWindow: public QWidget
{
  Q_OBJECT
  
  public:
  HelpWindow(QString filename, QWidget *parent = 0);
 
signals:

public slots:

private:
};
#endif  
