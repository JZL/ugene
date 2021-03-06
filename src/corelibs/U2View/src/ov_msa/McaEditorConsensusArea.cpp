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

#include <QToolBar>

#include <U2Algorithm/MSAConsensusAlgorithmRegistry.h>
#include <U2Algorithm/BuiltInConsensusAlgorithms.h>

#include <U2Core/AppContext.h>

#include <U2Gui/GUIUtils.h>

#include "McaEditorConsensusArea.h"
#include "McaEditor.h"
#include "MSAEditor.h" // for menu names consts
#include "view_rendering/McaConsensusAreaRenderer.h"

#include "ov_msa/MaConsensusMismatchController.h"

namespace U2 {

/************************************************************************/
/* McaEditorConsensusArea */
/************************************************************************/
McaEditorConsensusArea::McaEditorConsensusArea(McaEditorWgt *ui)
    : MaEditorConsensusArea(ui) {
    MSAConsensusAlgorithmFactory* algoFactory = AppContext::getMSAConsensusAlgorithmRegistry()->getAlgorithmFactory(BuiltInConsensusAlgorithms::LEVITSKY_ALGO);
    setConsensusAlgorithm(algoFactory);

    mismatchController = new MaConsensusMismatchController(this, consensusCache, editor);

    initRenderer();
    setupFontAndHeight();
}

void McaEditorConsensusArea::sl_buildStaticToolbar(GObjectView *, QToolBar *t) {
    t->addAction(mismatchController->getPrevAction());
    t->addAction(mismatchController->getNextAction());
}

void McaEditorConsensusArea::initRenderer() {
    renderer = new McaConsensusAreaRenderer(this);
}

void McaEditorConsensusArea::buildMenu(QMenu* m) {
    QMenu* copyMenu = GUIUtils::findSubMenu(m, MSAE_MENU_COPY);
    SAFE_POINT(copyMenu != NULL, "copyMenu", );
    copyMenu->addAction(copyConsensusAction);
    copyMenu->addAction(copyConsensusWithGapsAction);

    m->addAction(mismatchController->getNextAction());
    m->addAction(mismatchController->getPrevAction());
    m->addSeparator();
}

bool McaEditorConsensusArea::highlightConsensusChar(int pos) {
    return consensusSettings.highlightMismatches && mismatchController->isMismatch(pos);
}

} // namespace
