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

#include <QFontDialog>

#include <U2Algorithm/MsaHighlightingScheme.h>

#include <U2Core/AppContext.h>
#include <U2Core/Counter.h>
#include <U2Core/DocumentModel.h>
#include <U2Core/MultipleSequenceAlignmentObject.h>
#include <U2Core/SaveDocumentTask.h>
#include <U2Core/Settings.h>
#include <U2Core/U2DbiUtils.h>
#include <U2Core/U2OpStatusUtils.h>

#include <U2Gui/ExportDocumentDialogController.h>
#include <U2Gui/ExportObjectUtils.h>
#include <U2Gui/GUIUtils.h>

#include <U2View/MSAEditorSequenceArea.h>
#include <U2View/MSAEditorOffsetsView.h>
#include <U2View/MSAEditorOverviewArea.h>

#include "MaEditor.h"
#include "helpers/ScrollController.h"
#include "view_rendering/MaEditorWgt.h"

namespace U2 {

SNPSettings::SNPSettings()
    : seqId(U2MsaRow::INVALID_ROW_ID) {
}

const float MaEditor::zoomMult = 1.25;

MaEditor::MaEditor(GObjectViewFactoryId factoryId, const QString &viewName, GObject *obj)
    : GObjectView(factoryId, viewName),
      ui(NULL),
      cachedColumnWidth(0),
      exportHighlightedAction(NULL)
{
    maObject = qobject_cast<MultipleAlignmentObject*>(obj);
    objects.append(maObject);

    onObjectAdded(maObject);

    requiredObjects.append(maObject);
    GCOUNTER(cvar,tvar,factoryId);

    if (!U2DbiUtils::isDbiReadOnly(maObject->getEntityRef().dbiRef)) {
        U2OpStatus2Log os;
        maObject->setTrackMod(os, TrackOnUpdate);
    }

    // SANGER_TODO: move to separate method
    // do that in createWidget along with initActions?
    saveAlignmentAction = new QAction(QIcon(":core/images/msa_save.png"), tr("Save alignment"), this);
    saveAlignmentAction->setObjectName("Save alignment");
    connect(saveAlignmentAction, SIGNAL(triggered()), SLOT(sl_saveAlignment()));

    saveAlignmentAsAction = new QAction(QIcon(":core/images/msa_save_as.png"), tr("Save alignment as"), this);
    saveAlignmentAsAction->setObjectName("Save alignment as");
    connect(saveAlignmentAsAction, SIGNAL(triggered()), SLOT(sl_saveAlignmentAs()));

    zoomInAction = new QAction(QIcon(":core/images/zoom_in.png"), tr("Zoom In"), this);
    zoomInAction->setObjectName("Zoom In");
    connect(zoomInAction, SIGNAL(triggered()), SLOT(sl_zoomIn()));

    zoomOutAction = new QAction(QIcon(":core/images/zoom_out.png"), tr("Zoom Out"), this);
    zoomOutAction->setObjectName("Zoom Out");
    connect(zoomOutAction, SIGNAL(triggered()), SLOT(sl_zoomOut()));

    zoomToSelectionAction = new QAction(QIcon(":core/images/zoom_reg.png"), tr("Zoom To Selection"), this);
    zoomToSelectionAction->setObjectName("Zoom To Selection");
    connect(zoomToSelectionAction, SIGNAL(triggered()), SLOT(sl_zoomToSelection()));

    resetZoomAction = new QAction(QIcon(":core/images/zoom_whole.png"), tr("Reset Zoom"), this);
    resetZoomAction->setObjectName("Reset Zoom");
    connect(resetZoomAction, SIGNAL(triggered()), SLOT(sl_resetZoom()));

    changeFontAction = new QAction(QIcon(":core/images/font.png"), tr("Change Font"), this);
    changeFontAction->setObjectName("Change Font");
    connect(changeFontAction, SIGNAL(triggered()), SLOT(sl_changeFont()));

    exportHighlightedAction = new QAction(tr("Export highlighted"), this);
    exportHighlightedAction->setObjectName("Export highlighted");
    connect(exportHighlightedAction, SIGNAL(triggered()), this, SLOT(sl_exportHighlighted()));
    exportHighlightedAction->setDisabled(true);

    connect(maObject, SIGNAL(si_lockedStateChanged()), SLOT(sl_lockedStateChanged()));
    connect(this, SIGNAL(si_zoomOperationPerformed(bool)), SLOT(sl_resetColumnWidthCache()));
    connect(this, SIGNAL(si_fontChanged(QFont)), SLOT(sl_resetColumnWidthCache()));

    Settings* s = AppContext::getSettings();
    zoomFactor = MOBJECT_DEFAULT_ZOOM_FACTOR;
    font.setFamily(s->getValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_FAMILY, MOBJECT_DEFAULT_FONT_FAMILY).toString());
    font.setPointSize(s->getValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_SIZE, MOBJECT_DEFAULT_FONT_SIZE).toInt());
    font.setItalic(s->getValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_ITALIC, false).toBool());
    font.setBold(s->getValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_BOLD, false).toBool());
    calcFontPixelToPointSizeCoef();

    if ( (font.pointSize() == MOBJECT_MIN_FONT_SIZE) && (zoomFactor < 1.0f) ) {
        resizeMode = ResizeMode_OnlyContent;
    } else {
        resizeMode = ResizeMode_FontAndContent;
    }
}

void MaEditor::buildStaticToolbar(QToolBar* tb) {
    // SANGER_TODO: the first three actions are not yet valid for MCA!!
    tb->addAction(ui->getCopyFormattedSelectionAction());

    tb->addAction(saveAlignmentAction);
    tb->addAction(saveAlignmentAsAction);

    tb->addAction(zoomInAction);
    tb->addAction(zoomOutAction);
    tb->addAction(zoomToSelectionAction);
    tb->addAction(resetZoomAction);

    tb->addAction(showOverviewAction);
    tb->addAction(changeFontAction);
}

void MaEditor::buildStaticMenu(QMenu* m) {
    addCopyMenu(m);
    addEditMenu(m);
    addExportMenu(m);
    addViewMenu(m);
    addLoadMenu(m);
    addAlignMenu(m);

    GObjectView::buildStaticMenu(m);

    GUIUtils::disableEmptySubmenus(m);
}

int MaEditor::getAlignmentLen() const {
    return maObject->getLength();
}

int MaEditor::getNumSequences() const {
    return maObject->getNumRows();
}

bool MaEditor::isAlignmentEmpty() const {
    return getAlignmentLen() == 0 || getNumSequences() == 0;
}

const QRect& MaEditor::getCurrentSelection() const {
    return ui->getSequenceArea()->getSelection().getRect();
}

int MaEditor::getRowContentIndent(int) const {
    return 0;
}

int MaEditor::getSequenceRowHeight() const {
    QFontMetrics fm(font, ui);
    return fm.height() * zoomMult;
}

int MaEditor::getColumnWidth() const {
    if (0 == cachedColumnWidth) {
        QFontMetrics fm(font, ui);
        cachedColumnWidth = fm.width('W') * zoomMult;

        cachedColumnWidth = (int)(cachedColumnWidth * zoomFactor);
        cachedColumnWidth = qMax(cachedColumnWidth, MOBJECT_MIN_COLUMN_WIDTH);

    }
    return cachedColumnWidth;
}

QVariantMap MaEditor::getHighlightingSettings(const QString &highlightingFactoryId) const {
    const QVariant v = snp.highlightSchemeSettings.value(highlightingFactoryId);
    if (v.isNull()) {
        return QVariantMap();
    } else {
        CHECK(v.type() == QVariant::Map, QVariantMap());
        return v.toMap();
    }
}

void MaEditor::saveHighlightingSettings( const QString &highlightingFactoryId, const QVariantMap &settingsMap /* = QVariant()*/ ) {
    snp.highlightSchemeSettings.insert(highlightingFactoryId, QVariant(settingsMap));
}

void MaEditor::setReference(qint64 sequenceId) {
    if(sequenceId == U2MsaRow::INVALID_ROW_ID){
        exportHighlightedAction->setDisabled(true);
    }else{
        exportHighlightedAction->setEnabled(true);
    }
    if(snp.seqId != sequenceId) {
        snp.seqId = sequenceId;
        emit si_referenceSeqChanged(sequenceId);
    }
    //REDRAW OTHER WIDGETS
}

void MaEditor::updateReference(){
    if(maObject->getRowPosById(snp.seqId) == -1){
        setReference(U2MsaRow::INVALID_ROW_ID);
    }
}

void MaEditor::resetCollapsibleModel() {
    MSACollapsibleItemModel *collapsibleModel = ui->getCollapseModel();
    SAFE_POINT(NULL != collapsibleModel, "NULL collapsible model!", );
    collapsibleModel->reset();
}

void MaEditor::sl_zoomIn() {
    int pSize = font.pointSize();

    if (resizeMode == ResizeMode_OnlyContent) {
        zoomFactor *= zoomMult;
    } else if ( (pSize < MOBJECT_MAX_FONT_SIZE) && (resizeMode == ResizeMode_FontAndContent) ) {
        font.setPointSize(pSize+1);
        setFont(font);
    }

    bool resizeModeChanged = false;

    if (zoomFactor >= 1) {
        ResizeMode oldMode = resizeMode;
        resizeMode = ResizeMode_FontAndContent;
        resizeModeChanged = resizeMode != oldMode;
        zoomFactor = 1;
    }
    updateActions();

    emit si_zoomOperationPerformed(resizeModeChanged);
}

void MaEditor::sl_zoomOut() {
    int pSize = font.pointSize();

    bool resizeModeChanged = false;

    if (pSize > MOBJECT_MIN_FONT_SIZE) {
        font.setPointSize(pSize-1);
        setFont(font);
    } else {
        SAFE_POINT(zoomMult > 0, QString("Incorrect value of MSAEditor::zoomMult"),);
        zoomFactor /= zoomMult;
        ResizeMode oldMode = resizeMode;
        resizeMode = ResizeMode_OnlyContent;
        resizeModeChanged = resizeMode != oldMode;
    }
    updateActions();

    emit si_zoomOperationPerformed(resizeModeChanged);
}

void MaEditor::sl_zoomToSelection()
{
    ResizeMode oldMode = resizeMode;
    int seqAreaWidth =  ui->getSequenceArea()->width();
    MaEditorSelection selection = ui->getSequenceArea()->getSelection();
    if (selection.isNull()) {
        return;
    }
    int selectionWidth = selection.width();
    float pixelsPerBase = (seqAreaWidth / float(selectionWidth)) * zoomMult;
    int fontPointSize = int(pixelsPerBase / fontPixelToPointSize);
    if (fontPointSize >= MOBJECT_MIN_FONT_SIZE) {
        if (fontPointSize > MOBJECT_MAX_FONT_SIZE) {
            fontPointSize = MOBJECT_MAX_FONT_SIZE;
        }
        font.setPointSize(fontPointSize);
        setFont(font);
        resizeMode = ResizeMode_FontAndContent;
        zoomFactor = 1;
    } else {
        if (font.pointSize() != MOBJECT_MIN_FONT_SIZE) {
            font.setPointSize(MOBJECT_MIN_FONT_SIZE);
            setFont(font);
        }
        zoomFactor = pixelsPerBase / (MOBJECT_MIN_FONT_SIZE * fontPixelToPointSize);
        resizeMode = ResizeMode_OnlyContent;
    }
    ui->getScrollController()->setFirstVisibleBase(selection.x());
    ui->getScrollController()->setFirstVisibleRowByNumber(selection.y());

    updateActions();

    emit si_zoomOperationPerformed(resizeMode != oldMode);
}

void MaEditor::sl_resetZoom() {
    QFont f = getFont();
    f.setPointSize(MOBJECT_DEFAULT_FONT_SIZE);
    setFont(f);
    zoomFactor = MOBJECT_DEFAULT_ZOOM_FACTOR;
    ResizeMode oldMode = resizeMode;
    resizeMode = ResizeMode_FontAndContent;
    emit si_zoomOperationPerformed(resizeMode != oldMode);

    updateActions();
}

void MaEditor::sl_saveAlignment(){
    AppContext::getTaskScheduler()->registerTopLevelTask(new SaveDocumentTask(maObject->getDocument()));
}

void MaEditor::sl_saveAlignmentAs(){

    Document* srcDoc = maObject->getDocument();
    if (srcDoc == NULL) {
        return;
    }
    if (!srcDoc->isLoaded()) {
        return;
    }

    QObjectScopedPointer<ExportDocumentDialogController> dialog = new ExportDocumentDialogController(srcDoc, ui);
    dialog->setAddToProjectFlag(true);
    dialog->setWindowTitle(tr("Save Alignment"));
    ExportObjectUtils::export2Document(dialog);
}

void MaEditor::sl_changeFont() {
    bool ok = false;
    // QFontDialog::DontUseNativeDialog - no color selector, affects only Mac OS
    QFont f = QFontDialog::getFont(&ok, font, widget, tr("Characters Font"), QFontDialog::DontUseNativeDialog);
    if (!ok) {
        return;
    }
    setFont(f);
    updateActions();
    emit si_completeUpdate();
}

void MaEditor::sl_lockedStateChanged() {
    updateActions();
}

void MaEditor::sl_exportHighlighted(){
    QObjectScopedPointer<ExportHighligtingDialogController> d = new ExportHighligtingDialogController(ui, (QWidget*)AppContext::getMainWindow()->getQMainWindow());
    d->exec();
    CHECK(!d.isNull(), );

    if (d->result() == QDialog::Accepted){
        AppContext::getTaskScheduler()->registerTopLevelTask(new ExportHighligtningTask(d.data(), ui->getSequenceArea()));
    }
}

void MaEditor::sl_resetColumnWidthCache() {
    cachedColumnWidth = 0;
}

void MaEditor::initActions() {
    saveScreenshotAction = new QAction(QIcon(":/core/images/cam2.png"), tr("Export as image"), this);
    saveScreenshotAction->setObjectName("Export as image");
    connect(saveScreenshotAction, SIGNAL(triggered()), ui, SLOT(sl_saveScreenshot()));

    showOverviewAction = new QAction(QIcon(":/core/images/msa_show_overview.png"), tr("Overview"), this);
    showOverviewAction->setObjectName("Show overview");
    showOverviewAction->setCheckable(true);
    showOverviewAction->setChecked(true);
    connect(showOverviewAction, SIGNAL(triggered()), ui->getOverviewArea(), SLOT(sl_show()));
}

void MaEditor::addCopyMenu(QMenu* m) {
    QMenu* cm = m->addMenu(tr("Copy/Paste"));
    cm->menuAction()->setObjectName(MSAE_MENU_COPY);
}

void MaEditor::addEditMenu(QMenu* m) {
    QMenu* em = m->addMenu(tr("Edit"));
    em->menuAction()->setObjectName(MSAE_MENU_EDIT);
}

void MaEditor::addExportMenu(QMenu* m) {
    QMenu* em = m->addMenu(tr("Export"));
    em->menuAction()->setObjectName(MSAE_MENU_EXPORT);
    em->addAction(exportHighlightedAction);
    if(!ui->getSequenceArea()->getCurrentHighlightingScheme()->getFactory()->isRefFree() &&
                getReferenceRowId() != U2MsaRow::INVALID_ROW_ID){
        exportHighlightedAction->setEnabled(true);
    }else{
        exportHighlightedAction->setDisabled(true);
    }
}

void MaEditor::addViewMenu(QMenu* m) {
    QMenu* em = m->addMenu(tr("View"));
    em->menuAction()->setObjectName(MSAE_MENU_VIEW);
    if (ui->getOffsetsViewController() != NULL) {
        em->addAction(ui->getOffsetsViewController()->getToggleColumnsViewAction());
    }
}

void MaEditor::addLoadMenu( QMenu* m ) {
    QMenu* lsm = m->addMenu(tr("Add"));
    lsm->menuAction()->setObjectName(MSAE_MENU_LOAD);
}

void MaEditor::addAlignMenu(QMenu* m) {
    QMenu* em = m->addMenu(tr("Align"));
    em->setIcon(QIcon(":core/images/align.png"));
    em->menuAction()->setObjectName(MSAE_MENU_ALIGN);
}

static void saveFont(const QFont& f) {
    Settings* s = AppContext::getSettings();
    s->setValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_FAMILY, f.family());
    s->setValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_SIZE, f.pointSize());
    s->setValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_ITALIC, f.italic());
    s->setValue(MOBJECT_SETTINGS_ROOT + MOBJECT_SETTINGS_FONT_BOLD, f.bold());
}

void MaEditor::setFont(const QFont& f) {
    int pSize = f.pointSize();
    font = f;
    calcFontPixelToPointSizeCoef();
    font.setPointSize(qBound(MOBJECT_MIN_FONT_SIZE, pSize, MOBJECT_MAX_FONT_SIZE));
    emit si_fontChanged(font);
    saveFont(font);
}

void MaEditor::calcFontPixelToPointSizeCoef() {
    QFontInfo info(font);
    fontPixelToPointSize = (float) info.pixelSize() / (float) info.pointSize();
}

void MaEditor::setFirstVisibleBase(int firstPos) {
    if (ui->getSequenceArea()->isPosInRange(firstPos)) {
        ui->getScrollController()->setFirstVisibleBase(firstPos);
    }
}

void MaEditor::setZoomFactor(float newZoomFactor) {
    zoomFactor = newZoomFactor;
    sl_resetColumnWidthCache();
}

void MaEditor::updateActions() {
    zoomInAction->setEnabled(font.pointSize() < MOBJECT_MAX_FONT_SIZE);
    zoomOutAction->setEnabled( getColumnWidth() > MOBJECT_MIN_COLUMN_WIDTH );
    zoomToSelectionAction->setEnabled( font.pointSize() < MOBJECT_MAX_FONT_SIZE);
    changeFontAction->setEnabled( resizeMode == ResizeMode_FontAndContent);
}

} // namespace
