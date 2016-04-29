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

#ifndef _U2_MSA_HIGHLIGHTING_SCHEME_DISAGREEMENTS_H_
#define _U2_MSA_HIGHLIGHTING_SCHEME_DISAGREEMENTS_H_

#include "MsaHighlightingScheme.h"

namespace U2 {

class MsaHighlightingSchemeDisagreements : public MsaHighlightingScheme{
    Q_OBJECT
public:
    MsaHighlightingSchemeDisagreements(QObject *parent, const MsaHighlightingSchemeFactory *factory, MAlignmentObject *maObj);

    void process(const char refChar, char &seqChar, QColor &color, bool &hightlight, int refCharColumn, int refCharRow) const;
};

class U2ALGORITHM_EXPORT MsaHighlightingSchemeDisagreementsFactory : public MsaHighlightingSchemeFactory {
public:
    MsaHighlightingSchemeDisagreementsFactory(QObject *parent, const QString &id, const QString &name, DNAAlphabetType alphabetType);

    MsaHighlightingScheme * create(QObject *parent, MAlignmentObject *maObj) const;
};

}   // nmaespace U2

#endif // _U2_MSA_HIGHLIGHTING_SCHEME_DISAGREEMENTS_H_
