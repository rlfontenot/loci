
#ifndef GETFILE_H
#define GETFILE_H
#include <QDomElement>
#include <QtGui>

class QComboBox;
class QDir;
class QLabel;
class QPushButton;
class QTableWidget;
class QTableWidgetItem;


class GetFileWindow  : public QDialog
{
  Q_OBJECT

 public:
  GetFileWindow(QString exp, QString& fileSelected, QWidget *parent = 0  );
 
  signals:
  void fileNameSelected(QString);
                                
 private slots:
  void browse();
  void find();
  void openFileOfItem(int row, int column);
public slots:
  void addDirectory(QString);
  void setFileName(QString);
private:
  QStringList findFiles(const QDir &directory, const QStringList &files,
                        const QString &text);
  void showFiles(const QDir &directory, const QStringList &files);
  QPushButton *createButton(const QString &text, const char *member);
  QComboBox *createComboBox(const QString &text = QString());
  void createFilesTable();
  
  QComboBox *fileComboBox;
  QComboBox *textComboBox;
  QComboBox *directoryComboBox;
  QLabel *fileLabel;
  QLabel *textLabel;
  QLabel *directoryLabel;
  QLabel *filesFoundLabel;
  QPushButton *browseButton;
  QTableWidget *filesTable;
  QString selectedFileName;
  QLabel* fileNameLabel;
  QStringList nameFilter;
};

class FindFileWindow  : public GetFileWindow
{
  Q_OBJECT
  
public:
  FindFileWindow(const QDomElement& my_elem, QString& fileSelected, QWidget *parent = 0  );
public slots:
void updateElem(QString);
private:
  QDomElement myelem;
}; 

#endif
