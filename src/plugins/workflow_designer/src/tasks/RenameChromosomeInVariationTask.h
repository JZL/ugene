/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2016 UniPro <ugene@unipro.ru>
 * http://ugene.net
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
#ifndef _U2_RENAME_CHROMOSOME_IN_VARIATION_TASK_H_
#define _U2_RENAME_CHROMOSOME_IN_VARIATION_TASK_H_

#include <U2Core/Task.h>

namespace U2 {

class DocumentFormat;
class GObject;
class LoadDocumentTask;
class SaveDocumentTask;
class U2VariantTrack;

class RenameChromosomeInVariationTask : public Task {
    Q_OBJECT
public:
    RenameChromosomeInVariationTask(const QList<GObject *> &objects, const QStringList &prefixesToReplace, const QString &prefixReplaceWith);

private:
    void run();

    bool replaceSequenceName(U2VariantTrack &variantTrack) const;

    const QList<GObject *> objects;
    const QStringList prefixesToReplace;
    const QString prefixReplaceWith;
};

class RenameChromosomeInVariationFileTask : public Task {
    Q_OBJECT
public:
    RenameChromosomeInVariationFileTask(const QString &srcFileUrl, const QString &dstFileUrl, const QStringList &prefixesToReplace, const QString &prefixReplaceWith);

    QString getDstFileUrl() const;

private:
    void prepare();
    QList<Task *> onSubTaskFinished(Task *subTask);

    Task * initRenameTask();
    Task * initSaveTask();

    QList<GObject *> getVariantTrackObjects();
    DocumentFormat * getFormat();

    const QString srcFileUrl;
    const QString dstFileUrl;
    const QStringList prefixesToReplace;
    const QString prefixReplaceWith;

    LoadDocumentTask *loadTask;
    RenameChromosomeInVariationTask *renameTask;
    SaveDocumentTask *saveTask;
};

}   // namespace U2

#endif // _U2_RENAME_CHROMOSOME_IN_VARIATION_TASK_H_
