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

#ifndef _U2_MCA_ROW_INNER_DATA_H_
#define _U2_MCA_ROW_INNER_DATA_H_

#include <U2Core/DNAChromatogramObject.h>
#include <U2Core/DNASequence.h>
#include <U2Core/U2Mca.h>
#include <U2Core/U2Sequence.h>

namespace U2 {

class McaRowMemoryData {
public:
    McaRowMemoryData();

    DNAChromatogram chromatogram;
    U2MsaRowGapModel gapModel;
    DNASequence sequence;
    U2Region workingArea;
    qint64 rowLength;
    QVariantMap additionalInfo;
};

class McaRowDatabaseData {
public:
    McaRowDatabaseData();

    U2Chromatogram chromatogram;
    U2Sequence sequence;
    U2MsaRowGapModel gapModel;
    U2Region workingArea;
    qint64 rowLength;
    QVariantMap additionalInfo;
};

}   // namespace U2

#endif // _U2_MCA_ROW_INNER_DATA_H_