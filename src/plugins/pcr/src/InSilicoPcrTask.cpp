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

#include <U2Core/AppContext.h>
#include <U2Core/Counter.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNATranslation.h>
#include <U2Core/L10n.h>

#include "PrimerStatistics.h"

#include "InSilicoPcrTask.h"

namespace U2 {

InSilicoPcrTaskSettings::InSilicoPcrTaskSettings()
: isCircular(false), forwardMismatches(0), reverseMismatches(0), maxProductSize(0)
{

}

InSilicoPcrProduct::InSilicoPcrProduct()
: ta(0.0), forwardPrimerMatchLength(0), reversePrimerMatchLength(0)
{

}

InSilicoPcrTask::InSilicoPcrTask(const InSilicoPcrTaskSettings &settings)
: Task(tr("In Silico PCR"), TaskFlags(TaskFlag_ReportingIsSupported) | TaskFlag_ReportingIsEnabled | TaskFlags_FOSE_COSC),
settings(settings), forwardSearch(NULL), reverseSearch(NULL)
{
    GCOUNTER(cvar, tvar, "InSilicoPcrTask");
}

namespace {
    int getMaxError(const QByteArray &primer, uint mismatches) {
        return qMin(int(mismatches), primer.length() / 2);
    }
}

FindAlgorithmTaskSettings InSilicoPcrTask::getFindPatternSettings(U2Strand::Direction direction) {
    FindAlgorithmTaskSettings result;
    const DNAAlphabet *alphabet = AppContext::getDNAAlphabetRegistry()->findById(BaseDNAAlphabetIds::NUCL_DNA_DEFAULT());
    SAFE_POINT_EXT(NULL != alphabet, setError(L10N::nullPointerError("DNA Alphabet")), result);
    DNATranslation *translator = AppContext::getDNATranslationRegistry()->lookupComplementTranslation(alphabet);
    SAFE_POINT_EXT(NULL != translator, setError(L10N::nullPointerError("DNA Translator")), result);

    result.sequence = settings.sequence;
    result.searchIsCircular = settings.isCircular;
    result.searchRegion.length = settings.sequence.length();
    result.patternSettings = FindAlgorithmPatternSettings_Subst;
    result.strand = FindAlgorithmStrand_Both;

    if (U2Strand::Direct == direction) {
        result.pattern = settings.forwardPrimer;
        result.maxErr = getMaxError(settings.forwardPrimer, settings.forwardMismatches);
    } else {
        result.pattern = settings.reversePrimer;
        result.maxErr = getMaxError(settings.reversePrimer, settings.reverseMismatches);
    }

    result.complementTT = translator;

    return result;
}

void InSilicoPcrTask::prepare() {
    FindAlgorithmTaskSettings forwardSettings = getFindPatternSettings(U2Strand::Direct);
    CHECK_OP(stateInfo, );
    FindAlgorithmTaskSettings reverseSettings = getFindPatternSettings(U2Strand::Complementary);
    CHECK_OP(stateInfo, );
    forwardSearch = new FindAlgorithmTask(forwardSettings);
    reverseSearch = new FindAlgorithmTask(reverseSettings);
    addSubTask(forwardSearch);
    addSubTask(reverseSearch);
}

void InSilicoPcrTask::run() {
    QList<FindAlgorithmResult> forwardResults = forwardSearch->popResults();
    QList<FindAlgorithmResult> reverseResults = reverseSearch->popResults();
    algoLog.details(tr("Forward primers found: %1").arg(forwardResults.size()));
    algoLog.details(tr("Reverse primers found: %1").arg(reverseResults.size()));

    int minProductSize = qMax(settings.forwardPrimer.length(), settings.reversePrimer.length());
    foreach (const FindAlgorithmResult &forward, forwardResults) {
        foreach (const FindAlgorithmResult &reverse, reverseResults) {
            CHECK(!isCanceled(), );
            if (forward.strand == reverse.strand) {
                continue;
            }
            U2Region left = forward.region;
            U2Region right = reverse.region;
            if (forward.strand.isCompementary()) {
                qSwap(left, right);
            }
            qint64 productSize = getProductSize(left, right);
            if (productSize >= minProductSize && productSize <= qint64(settings.maxProductSize)) {
                InSilicoPcrProduct product = createResult(left, U2Region(left.startPos, productSize), right, forward.strand.getDirection());
                results << product;
            }
        }
    }
}

QString InSilicoPcrTask::generateReport() const {
    PrimersPairStatistics calc(settings.forwardPrimer, settings.reversePrimer);
    QString spaces;
    for (int i=0; i<150; i++) {
        spaces += "&nbsp;";
    }
    return tr("Products found: %1").arg(results.size()) +
           "<br>" +
           spaces +
           "<br>" +
           tr("Primers details:") +
           calc.generateReport();
}

InSilicoPcrProduct InSilicoPcrTask::createResult(const U2Region &leftPrimer, const U2Region &product, const U2Region &rightPrimer, U2Strand::Direction direction) const {
    QByteArray productSequence = settings.sequence.mid(product.startPos, product.length);
    if (productSequence.length() < product.length) {
        assert(settings.isCircular);
        productSequence += settings.sequence.left(product.endPos() - settings.sequence.length());
    }

    InSilicoPcrProduct result;
    result.region = product;
    result.ta = PrimerStatistics::getAnnealingTemperature(productSequence, settings.forwardPrimer, settings.reversePrimer);
    result.forwardPrimerMatchLength = leftPrimer.length;
    result.reversePrimerMatchLength = rightPrimer.length;
    result.forwardPrimer = settings.forwardPrimer;
    result.reversePrimer = settings.reversePrimer;
    if (U2Strand::Complementary == direction) {
        qSwap(result.forwardPrimer, result.reversePrimer);
    }
    return result;
}

qint64 InSilicoPcrTask::getProductSize(const U2Region &left, const U2Region &right) const {
    qint64 result = right.endPos() - left.startPos;
    if (result < 0 && settings.isCircular) {
        return right.endPos() + (settings.sequence.length() - left.startPos);
    }
    return result;
}

const QList<InSilicoPcrProduct> & InSilicoPcrTask::getResults() const {
    return results;
}

const InSilicoPcrTaskSettings & InSilicoPcrTask::getSettings() const {
    return settings;
}

} // U2
