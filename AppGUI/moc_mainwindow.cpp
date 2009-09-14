/****************************************************************************
** Meta object code from reading C++ file 'mainwindow.h'
**
** Created: Mon Sep 14 09:59:16 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "mainwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'mainwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_showDelegate[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_showDelegate[] = {
    "showDelegate\0"
};

const QMetaObject showDelegate::staticMetaObject = {
    { &QItemDelegate::staticMetaObject, qt_meta_stringdata_showDelegate,
      qt_meta_data_showDelegate, 0 }
};

const QMetaObject *showDelegate::metaObject() const
{
    return &staticMetaObject;
}

void *showDelegate::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_showDelegate))
        return static_cast<void*>(const_cast< showDelegate*>(this));
    return QItemDelegate::qt_metacast(_clname);
}

int showDelegate::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QItemDelegate::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_colorDelegate[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_colorDelegate[] = {
    "colorDelegate\0"
};

const QMetaObject colorDelegate::staticMetaObject = {
    { &QItemDelegate::staticMetaObject, qt_meta_stringdata_colorDelegate,
      qt_meta_data_colorDelegate, 0 }
};

const QMetaObject *colorDelegate::metaObject() const
{
    return &staticMetaObject;
}

void *colorDelegate::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_colorDelegate))
        return static_cast<void*>(const_cast< colorDelegate*>(this));
    return QItemDelegate::qt_metacast(_clname);
}

int colorDelegate::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QItemDelegate::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_MainWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
      32,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      12,   11,   11,   11, 0x05,
      36,   11,   11,   11, 0x05,
      51,   11,   11,   11, 0x05,

 // slots: signature, parameters, type, tag, flags
      71,   11,   11,   11, 0x0a,
      81,   11,   11,   11, 0x0a,
     100,   92,   11,   11, 0x0a,
     127,  122,   11,   11, 0x0a,
     158,   11,  153,   11, 0x0a,
     175,   11,   11,   11, 0x0a,
     181,   11,   11,   11, 0x0a,
     195,   11,   11,   11, 0x0a,
     205,   11,   11,   11, 0x0a,
     215,   11,  153,   11, 0x0a,
     225,   11,  153,   11, 0x0a,
     235,   11,  153,   11, 0x0a,
     247,   11,   11,   11, 0x0a,
     265,   11,   11,   11, 0x0a,
     284,   11,   11,   11, 0x0a,
     301,   11,   11,   11, 0x0a,
     318,   11,   11,   11, 0x0a,
     331,   11,   11,   11, 0x0a,
     346,  344,   11,   11, 0x0a,
     384,   11,   11,   11, 0x0a,
     399,   11,   11,   11, 0x0a,
     426,   11,   11,   11, 0x0a,
     445,   11,   11,   11, 0x0a,
     462,   11,   11,   11, 0x0a,
     473,   11,   11,   11, 0x0a,
     495,   11,   11,   11, 0x0a,
     516,   11,   11,   11, 0x0a,
     533,   11,   11,   11, 0x0a,
     551,   11,   11,   11, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_MainWindow[] = {
    "MainWindow\0\0setCurrent(QModelIndex)\0"
    "stateChanged()\0componentsChanged()\0"
    "newCase()\0openCase()\0theelem\0"
    "setGrid(QDomElement&)\0elem\0"
    "setBoundary(QDomElement&)\0bool\0"
    "selectBoundary()\0cut()\0resetSlider()\0"
    "openSca()\0openVog()\0saveVar()\0saveXml()\0"
    "saveImage()\0aboutPreprocess()\0"
    "aboutPostprocess()\0showDisplayBar()\0"
    "hideDisplayBar()\0showVisBar()\0"
    "hideVisBar()\0,\0showBoundary(QModelIndex,QModelIndex)\0"
    "toggleViewer()\0setCurrentObj(QModelIndex)\0"
    "selectCurrent(int)\0bdWindowClosed()\0"
    "snapshot()\0updateStatus(QString)\0"
    "updateStatusTip(int)\0clearAllStatus()\0"
    "clearLastStatus()\0changePage(int)\0"
};

const QMetaObject MainWindow::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_MainWindow,
      qt_meta_data_MainWindow, 0 }
};

const QMetaObject *MainWindow::metaObject() const
{
    return &staticMetaObject;
}

void *MainWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_MainWindow))
        return static_cast<void*>(const_cast< MainWindow*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int MainWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: setCurrent((*reinterpret_cast< QModelIndex(*)>(_a[1]))); break;
        case 1: stateChanged(); break;
        case 2: componentsChanged(); break;
        case 3: newCase(); break;
        case 4: openCase(); break;
        case 5: setGrid((*reinterpret_cast< QDomElement(*)>(_a[1]))); break;
        case 6: setBoundary((*reinterpret_cast< QDomElement(*)>(_a[1]))); break;
        case 7: { bool _r = selectBoundary();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 8: cut(); break;
        case 9: resetSlider(); break;
        case 10: openSca(); break;
        case 11: openVog(); break;
        case 12: { bool _r = saveVar();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 13: { bool _r = saveXml();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 14: { bool _r = saveImage();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 15: aboutPreprocess(); break;
        case 16: aboutPostprocess(); break;
        case 17: showDisplayBar(); break;
        case 18: hideDisplayBar(); break;
        case 19: showVisBar(); break;
        case 20: hideVisBar(); break;
        case 21: showBoundary((*reinterpret_cast< QModelIndex(*)>(_a[1])),(*reinterpret_cast< QModelIndex(*)>(_a[2]))); break;
        case 22: toggleViewer(); break;
        case 23: setCurrentObj((*reinterpret_cast< QModelIndex(*)>(_a[1]))); break;
        case 24: selectCurrent((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 25: bdWindowClosed(); break;
        case 26: snapshot(); break;
        case 27: updateStatus((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 28: updateStatusTip((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 29: clearAllStatus(); break;
        case 30: clearLastStatus(); break;
        case 31: changePage((*reinterpret_cast< int(*)>(_a[1]))); break;
        }
        _id -= 32;
    }
    return _id;
}

// SIGNAL 0
void MainWindow::setCurrent(QModelIndex _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void MainWindow::stateChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 1, 0);
}

// SIGNAL 2
void MainWindow::componentsChanged()
{
    QMetaObject::activate(this, &staticMetaObject, 2, 0);
}
QT_END_MOC_NAMESPACE
