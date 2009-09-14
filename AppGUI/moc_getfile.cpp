/****************************************************************************
** Meta object code from reading C++ file 'getfile.h'
**
** Created: Mon Sep 14 09:59:16 2009
**      by: The Qt Meta Object Compiler version 59 (Qt 4.4.3)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "getfile.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'getfile.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.4.3. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_GetFileWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       4,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // signals: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x05,

 // slots: signature, parameters, type, tag, flags
      41,   14,   14,   14, 0x08,
      50,   14,   14,   14, 0x08,
      68,   57,   14,   14, 0x08,

       0        // eod
};

static const char qt_meta_stringdata_GetFileWindow[] = {
    "GetFileWindow\0\0fileNameSelected(QString)\0"
    "browse()\0find()\0row,column\0"
    "openFileOfItem(int,int)\0"
};

const QMetaObject GetFileWindow::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_GetFileWindow,
      qt_meta_data_GetFileWindow, 0 }
};

const QMetaObject *GetFileWindow::metaObject() const
{
    return &staticMetaObject;
}

void *GetFileWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_GetFileWindow))
        return static_cast<void*>(const_cast< GetFileWindow*>(this));
    return QDialog::qt_metacast(_clname);
}

int GetFileWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: fileNameSelected((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 1: browse(); break;
        case 2: find(); break;
        case 3: openFileOfItem((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        }
        _id -= 4;
    }
    return _id;
}

// SIGNAL 0
void GetFileWindow::fileNameSelected(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
static const uint qt_meta_data_FindFileWindow[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      16,   15,   15,   15, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_FindFileWindow[] = {
    "FindFileWindow\0\0updateElem(QString)\0"
};

const QMetaObject FindFileWindow::staticMetaObject = {
    { &GetFileWindow::staticMetaObject, qt_meta_stringdata_FindFileWindow,
      qt_meta_data_FindFileWindow, 0 }
};

const QMetaObject *FindFileWindow::metaObject() const
{
    return &staticMetaObject;
}

void *FindFileWindow::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_FindFileWindow))
        return static_cast<void*>(const_cast< FindFileWindow*>(this));
    return GetFileWindow::qt_metacast(_clname);
}

int FindFileWindow::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = GetFileWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: updateElem((*reinterpret_cast< QString(*)>(_a[1]))); break;
        }
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
