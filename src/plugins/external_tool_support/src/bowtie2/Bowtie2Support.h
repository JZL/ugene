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

#ifndef _U2_BOWTIE2_SUPPORT_H_
#define _U2_BOWTIE2_SUPPORT_H_

#include <U2Core/ExternalToolRegistry.h>

#define ET_BOWTIE2_ALIGN "Bowtie 2 aligner"
#define ET_BOWTIE2_BUILD "Bowtie 2 build indexer"
#define ET_BOWTIE2_INSPECT "Bowtie 2 index inspector"

#define BOWTIE2_TMP_DIR "bowtie2"

namespace U2 {

class Bowtie2Support : public ExternalTool
{
    Q_OBJECT

public:
    Bowtie2Support(const QString& name, const QString& path = "");
};

} // namespace U2

#endif
