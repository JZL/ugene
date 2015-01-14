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

#ifndef _U2_WORKFLOW_TASKS_REGISTRY_
#define _U2_WORKFLOW_TASKS_REGISTRY_

#include <U2Lang/ReadDocumentTaskFactory.h>

namespace U2 {
namespace Workflow {

class U2LANG_EXPORT WorkflowTasksRegistry {
public:
    WorkflowTasksRegistry();
    ~WorkflowTasksRegistry();

    /**
     * Returns true if the factory was successfully registered.
     * Returns false if the factory with this id already exists
     */
    bool registerReadDocumentTaskFactory(ReadDocumentTaskFactory *factory);
    ReadDocumentTaskFactory *getReadDocumentTaskFactory(const QString &id);

private:
    QMap<QString, ReadDocumentTaskFactory*> readTasks;
};

} // Workflow
} // U2

#endif // _U2_WORKFLOW_TASKS_REGISTRY_
