/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2017 UniPro <ugene@unipro.ru>
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

#ifndef _U2_MCA_EDITOR_WGT_H_
#define _U2_MCA_EDITOR_WGT_H_

#include "view_rendering/MaEditorWgt.h"

namespace U2 {

class McaEditor;
class McaEditorReferenceArea;
class McaEditorSequenceArea;
class McaReferenceCharController;

class McaEditorWgt : public MaEditorWgt {
    Q_OBJECT
public:
    McaEditorWgt(McaEditor* editor);

    McaEditor* getEditor() const;
    McaEditorSequenceArea* getSequenceArea() const;
    McaReferenceCharController* getRefCharController() const;

protected:
    void initSeqArea(GScrollBar* shBar, GScrollBar* cvBar);
    void initOverviewArea();
    void initNameList(QScrollBar* nhBar);
    void initConsensusArea();
    void initStatusBar();

private:
    McaEditorReferenceArea*     refArea;
    McaReferenceCharController* refCharController;
};

}   // namespace U2

#endif // _U2_MCA_EDITOR_WGT_H_
