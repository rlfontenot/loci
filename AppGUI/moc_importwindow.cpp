/****************************************************************************
** Meta object code from reading C++ file 'importwindow.h'
**
** Created: Mon Sep 14 09:59:17 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "importwindow.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'importwindow.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_VogOption[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_VogOption[] = {
    "VogOption\0\0unitButtonClicked(int)\0"
};

const QMetaObject VogOption::staticMetaObject = {
    { &QGroupBox::staticMetaObject, qt_meta_stringdata_VogOption,
      qt_meta_data_VogOption, 0 }
};

const QMetaObject *VogOption::metaObject() const
{
    return &staticMetaObject;
}

void *VogOption::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_VogOption))
        return static_cast<void*>(const_cast< VogOption*>(this));
    return QGroupBox::qt_metacast(_clname);
}

int VogOption::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGroupBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: unitButtonClicked((*reinterpret_cast< int(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_XdrOption[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      11,   10,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_XdrOption[] = {
    "XdrOption\0\0update(int)\0"
};

const QMetaObject XdrOption::staticMetaObject = {
    { &QGroupBox::staticMetaObject, qt_meta_stringdata_XdrOption,
      qt_meta_data_XdrOption, 0 }
};

const QMetaObject *XdrOption::metaObject() const
{
    return &staticMetaObject;
}

void *XdrOption::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_XdrOption))
        return static_cast<void*>(const_cast< XdrOption*>(this));
    return QGroupBox::qt_metacast(_clname);
}

int XdrOption::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QGroupBox::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: update((*reinterpret_cast< int(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}
static const uint qt_meta_data_ImportWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      31,   14,   13,   13, 0x0a,
      77,   13,   13,   13, 0x0a,
      87,   13,   13,   13, 0x0a,
      95,   13,   13,   13, 0x0a,
     117,   13,   13,   13, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_ImportWindow[] = {
    "ImportWindow\0\0current,previous\0"
    "changePage(QListWidgetItem*,QListWidgetItem*)\0"
    "convert()\0check()\0usuageButtonClicked()\0"
    "updateFileName(QString)\0"
};

const QMetaObject ImportWindow::staticMetaObject = {
    { &GeneralWindow::staticMetaObject, qt_meta_stringdata_ImportWindow,
      qt_meta_data_ImportWindow, 0 }
};

const QMetaObject *ImportWindow::metaObject() const
{
    return &staticMetaObject;
}

void *ImportWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_ImportWindow))
        return static_cast<void*>(const_cast< ImportWindow*>(this));
    return GeneralWindow::qt_metacast(_clname);
}

int ImportWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GeneralWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: changePage((*reinterpret_cast< QListWidgetItem*(*)>(_a[1])),(*reinterpret_cast< QListWidgetItem*(*)>(_a[2]))); break;
        case 1: convert(); break;
        case 2: check(); break;
        case 3: usuageButtonClicked(); break;
        case 4: updateFileName((*reinterpret_cast< QString(*)>(_a[1]))); break;
        }
        _id -= 5;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
