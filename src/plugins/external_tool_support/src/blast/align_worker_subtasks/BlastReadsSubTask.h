/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2016 UniPro <ugene@unipro.ru>
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

#ifndef _U2_BLAST_READS_SUBTASK_H_
#define _U2_BLAST_READS_SUBTASK_H_

#include <U2Core/Task.h>

#include <U2Lang/DbiDataHandler.h>
#include <U2Lang/DbiDataStorage.h>


namespace U2 {

class BlastAllSupportTask;

namespace Workflow {

/************************************************************************/
/* BlastReadsSubTask */
/************************************************************************/
class BlastReadsSubTask : public Task {
    Q_OBJECT
public:
    BlastReadsSubTask(const QString& dbPath,
                      const QList<SharedDbiDataHandler> &reads,
                      DbiDataStorage *storage);

    void prepare();
private:
    const QString dbPath;
    const QList<SharedDbiDataHandler> reads;

    DbiDataStorage *storage;
};

/************************************************************************/
/* BlastAndSwReadTask */
/************************************************************************/
class BlastAndSwReadTask : public Task {
    Q_OBJECT
public:
    BlastAndSwReadTask(const QString& dbPath,
                       const SharedDbiDataHandler& read,
                       DbiDataStorage *storage);

    void prepare();
    QList<Task*> onSubTaskFinished(Task *subTask);
    void run();
private:
    QByteArray getRead();

    const QString dbPath;
    const SharedDbiDataHandler read;
    DbiDataStorage *storage;

    BlastAllSupportTask*    blastTask;

    QString blastResultDir;
};

} // namespace Workflow
} // namespace U2

#endif // _U2_BLAST_READS_SUBTASK_H_
