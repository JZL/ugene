/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2015 UniPro <ugene@unipro.ru>
 * http://ugene.unipro.ru
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include <QtCore/QEventLoop>
#include <QtCore/QTimer>
#include <QtTest/QSpontaneKeyEvent>
#include <QtTest>
#if (QT_VERSION < 0x050000) //Qt 5
#include <QtGui/QApplication>
#include <QtGui/QDesktopWidget>
#include <QtGui/QPixmap>
#else
#include <QtWidgets/QApplication>
#include <QtWidgets/QDesktopWidget>
#include <QtGui/QScreen>
#endif

#include "api/GTGlobals.h"
#include <U2Core/AppContext.h>

#ifdef Q_OS_WIN
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace U2 {
namespace {
    void sysSleep(int sec) {
#ifdef Q_OS_WIN
        Sleep(1000*sec);
#else
        sleep(sec);
#endif
    }
}

void GTGlobals::sleep(int msec) {
    QTest::qWait((msec));
}

void GTGlobals::systemSleep(int sec) {
    sysSleep(sec);
}


void GTGlobals::sendEvent(QObject *obj, QEvent *e) {
    QSpontaneKeyEvent::setSpontaneous(e);
    qApp->notify(obj, e);
}

void GTGlobals::takeScreenShot(QString path) {
#if (QT_VERSION < 0x050000) // deprecated method
    QPixmap originalPixmap = QPixmap::grabWindow(QApplication::desktop()->winId());
#else
    QPixmap originalPixmap = QGuiApplication::primaryScreen()->grabWindow(QApplication::desktop()->winId());
#endif
    originalPixmap.save(path);
}

GTGlobals::FindOptions::FindOptions(bool fail) :
    depth(INFINITE_DEPTH),
    failIfNull(fail)
{
}

void GTGlobals::GUITestFail(){
    uiLog.trace("\nGT_DEBUG_MESSAGE !!!FIRST FAIL");
}

} //namespace

