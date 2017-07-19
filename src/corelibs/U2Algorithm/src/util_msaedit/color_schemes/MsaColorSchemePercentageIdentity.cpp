/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2017 UniPro <ugene@unipro.ru>
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

#include <QColor>

#include <U2Algorithm/MSAConsensusUtils.h>

#include <U2Core/MultipleAlignment.h>
#include <U2Core/MultipleAlignmentObject.h>

#include "MsaColorSchemePercentageIdentity.h"

namespace U2 {

MsaColorSchemePercentageIdentity::MsaColorSchemePercentageIdentity(QObject *parent, const MsaColorSchemeFactory *factory, MultipleAlignmentObject *maObj)
    : MsaColorScheme(parent, factory, maObj),
      cacheVersion(0),
      objVersion(1)
{
    mask4[0]=81;
    mask4[1]=61;
    mask4[2]=41;
    mask4[3]=25;

    colorsByRange[0] = QColor("#6464FF");
    colorsByRange[1] = QColor("#9999FF");
    colorsByRange[2] = QColor("#CCCCFF");
    colorsByRange[3] = QColor();

    connect(maObj, SIGNAL(si_alignmentChanged(const MultipleAlignment &, const MaModificationInfo &)), SLOT(sl_alignmentChanged()));
}

QColor MsaColorSchemePercentageIdentity::getColor(int /*seq*/, int pos, char c) const {
    if (c == U2Msa::GAP_CHAR) {
        //I think grey is clearer because it shows it's not nothing
        return QColor(250, 250, 250);
    }

    const MultipleAlignment ma = maObj->getMultipleAlignment();
    int nSeq = ma->getNumRows();
    for(int seq = 0; seq<nSeq;seq++){
            char seqC = ma->charAt(seq, pos);
            if(seqC == c){
                    /*
                            //If no Hsl, can use a bitmasked rgb value
                            tmpColor = colorFull/(seq+1);
                            return QColor(tmpColor>>(4*4), tmpColor>>(4*2)&0xff, tmpColor&0xff);
                    */
                    //-30 bc I don't want the wrap around of red on both sides
                    return QColor().fromHsl((((358-30)*(seq+1))/(nSeq)), 255, 128);
            }
    }
    //Should never return because, at least, the sequence should match itself
    return QColor();
}

void MsaColorSchemePercentageIdentity::sl_alignmentChanged() {
    objVersion++;
}

void MsaColorSchemePercentageIdentity::updateCache() const {
    if (cacheVersion == objVersion) {
        return;
    }
    const MultipleAlignment msa = maObj->getMultipleAlignment();
    int aliLen = msa->getLength();
    indentCache.resize(aliLen);
    for (int i = 0; i < aliLen; i++) {
        indentCache[i] = MSAConsensusUtils::packConsensusCharsToInt(msa, i, mask4, true);
    }
    cacheVersion = objVersion;
}

MsaColorSchemePercentageIdentityFactory::MsaColorSchemePercentageIdentityFactory(QObject *parent, const QString &id, const QString &name, const AlphabetFlags &supportedAlphabets) : MsaColorSchemeFactory(parent, id, name, supportedAlphabets)
{

}

MsaColorScheme * MsaColorSchemePercentageIdentityFactory::create(QObject *parent, MultipleAlignmentObject *maObj) const {
    return new MsaColorSchemePercentageIdentity(parent, this, maObj);
}

}   // namespace U2
