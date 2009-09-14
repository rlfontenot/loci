/****************************************************************************
** Meta object code from reading C++ file 'physicswindow.h'
**
** Created: Mon Sep 14 09:59:17 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "physicswindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'physicswindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_FloatEditDelegate[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_FloatEditDelegate[] = {
    "FloatEditDelegate\0"
};

const QMetaObject FloatEditDelegate::staticMetaObject = {
    { &QItemDelegate::staticMetaObject, qt_meta_stringdata_FloatEditDelegate,
      qt_meta_data_FloatEditDelegate, 0 }
};

const QMetaObject *FloatEditDelegate::metaObject() const
{
    return &staticMetaObject;
}

void *FloatEditDelegate::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_FloatEditDelegate))
        return static_cast<void*>(const_cast< FloatEditDelegate*>(this));
    return QItemDelegate::qt_metacast(_clname);
}

int FloatEditDelegate::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QItemDelegate::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
static const uint qt_meta_data_CpWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      10,    9,    9,    9, 0x0a,
      41,    9,   36,    9, 0x0a,
      48,    9,    9,    9, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_CpWindow[] = {
    "CpWindow\0\0setNumberOfIntervals(int)\0"
    "bool\0save()\0replaceButtonToggled(bool)\0"
};

const QMetaObject CpWindow::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_CpWindow,
      qt_meta_data_CpWindow, 0 }
};

const QMetaObject *CpWindow::metaObject() const
{
    return &staticMetaObject;
}

void *CpWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_CpWindow))
        return static_cast<void*>(const_cast< CpWindow*>(this));
    return QWidget::qt_metacast(_clname);
}

int CpWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: setNumberOfIntervals((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 1: { bool _r = save();
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = _r; }  break;
        case 2: replaceButtonToggled((*reinterpret_cast< bool(*)>(_a[1]))); break;
        }
        _id -= 3;
    }
    return _id;
}
static const uint qt_meta_data_PhysicsWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x05,
      43,   14,   14,   14, 0x05,

 // slots: signature, parameters, type, tag, flags
      69,   14,   14,   14, 0x08,
      90,   14,   14,   14, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_PhysicsWindow[] = {
    "PhysicsWindow\0\0parentStateChanged(QString)\0"
    "stateUpdated(QStringList)\0"
    "updateState(QString)\0checkStatus()\0"
};

const QMetaObject PhysicsWindow::staticMetaObject = {
    { &GeneralWindow::staticMetaObject, qt_meta_stringdata_PhysicsWindow,
      qt_meta_data_PhysicsWindow, 0 }
};

const QMetaObject *PhysicsWindow::metaObject() const
{
    return &staticMetaObject;
}

void *PhysicsWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_PhysicsWindow))
        return static_cast<void*>(const_cast< PhysicsWindow*>(this));
    return GeneralWindow::qt_metacast(_clname);
}

int PhysicsWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: parentStateChanged((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: stateUpdated((*reinterpret_cast< const QStringList(*)>(_a[1]))); break;
        case 2: updateState((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 3: checkStatus(); break;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void PhysicsWindow::parentStateChanged(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void PhysicsWindow::stateUpdated(const QStringList & _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_END_MOC_NAMESPACE
