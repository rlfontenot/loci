

#include <QHBoxLayout>
#include <QTextEdit>
#include <QDir>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include "helpwindow.h"

HelpWindow::HelpWindow(QString filename, QWidget* parent):QWidget(parent){
  setWindowTitle(tr("Help Window"));
  char* resourcepath = getenv("CHEMDEMOPATH");
  if(resourcepath){
    filename = QString(resourcepath) + "/doc/html/" + filename;
    QWidget::setAttribute(Qt::WA_DeleteOnClose);
    setWindowTitle("help window");
    QHBoxLayout* mainLayout = new QHBoxLayout;
    QTextEdit* edit = new QTextEdit;
    QFile file(filename);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
      QMessageBox::information(window(), "help window",
                               tr("Cannot open ") + filename + tr(" for reading!"));
    }
  
    QTextStream in(&file);
    QString text = in.readAll();
    edit->setReadOnly(true);
    edit->setHtml(text);
    mainLayout->addWidget(edit);
    setLayout(mainLayout);
    setMinimumSize(700, 700);
  }
}




