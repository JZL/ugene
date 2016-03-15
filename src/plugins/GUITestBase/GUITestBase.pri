#include (GUITestBase.pri)

PLUGIN_ID=GUITestBase
PLUGIN_NAME=GUI Test Base
PLUGIN_VENDOR=Unipro
include( ../../ugene_plugin_common.pri )

QT += testlib
greaterThan(QT_MAJOR_VERSION, 4): QT += webkitwidgets

INCLUDEPATH += ../../corelibs/U2View/_tmp/ ../../libs_3rdparty/QSpec/src
LIBS +=-L../../_release -lhumimit

!debug_and_release|build_pass {

    CONFIG(debug, debug|release) {
        LIBS -= -L../../_release -lhumimit
        LIBS += -L../../_debug -lhumimitd
    }
}

unix {
    !macx {
    LIBS += -lXtst
    }
    macx {
    QMAKE_LFLAGS += -framework ApplicationServices
    }
}

win32 {
    LIBS += User32.lib Gdi32.lib
}

macx {
    LIBS += -framework AppKit
}



