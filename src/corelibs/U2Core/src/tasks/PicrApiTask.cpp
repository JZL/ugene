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

#include <QtCore/QEventLoop>
#include <QtCore/QTimer>

#include <QtNetwork/QAuthenticator>

#include <QtXml/QDomDocument>

#include <U2Core/DocumentModel.h>
#include <U2Core/LoadDocumentTask.h>
#include <U2Core/AppContext.h>
#include <U2Core/AppSettings.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/U2AlphabetUtils.h>
#include <U2Core/U2SequenceUtils.h>
#include <U2Core/GObjectRelationRoles.h>

#include "LoadDASDocumentTask.h"
#include "PicrApiTask.h"

namespace U2 {

//////////////////////////////////////////////////////////////////////////
//PicrElement
PicrElement::PicrElement(const QString& _accessionNumber,
            const int _accessionVersion,
            const QString& _databaseName,
            const QString& _databaseDescription,
            const int _ncbiId,
            const int _taxonId,
            int _identity):
    accessionNumber(_accessionNumber),
    accessionVersion(_accessionVersion),
    databaseName(_databaseName),
    databaseDescription(_databaseDescription),
    ncbiId(_ncbiId),
    taxonId(_taxonId),
    identity(_identity) {
}

bool PicrElement::operator == (const PicrElement& other) const {
    return accessionNumber == other.accessionNumber &&
            databaseName == other.databaseName;
}

bool PicrElement::operator > (const PicrElement &other) const {
    return CommonDasSettings::getPriority(databaseName) > CommonDasSettings::getPriority(other.databaseName);
}

//////////////////////////////////////////////////////////////////////////
//XmlPicrIdsParser
const QString XmlPicrParser::PICR_ACCESSION_RETURN = "ns2:getUPIForAccessionReturn";
const QString XmlPicrParser::PICR_SEQUENCE_RETURN = "ns2:getUPIForSequenceReturn";
const QString XmlPicrParser::PICR_BLAST_RETURN = "ns2:getUPIForBlastReturn";
const QString XmlPicrParser::PICR_SEQUENCE = "sequence";
const QString XmlPicrParser::PICR_IDENTICAL_CROSS_REFERENCE = "identicalCrossReferences";
const QString XmlPicrParser::PICR_ACCESSION_NUMBER = "accession";
const QString XmlPicrParser::PICR_ACCESSION_VERSION = "accessionVersion";
const QString XmlPicrParser::PICR_DATABASE_DESCRIPTION = "databaseDescription";
const QString XmlPicrParser::PICR_DATABASE_NAME = "databaseName";
const QString XmlPicrParser::PICR_NCBI_GI = "gi";
const QString XmlPicrParser::PICR_TAXON_ID = "taxonId";
const QString XmlPicrParser::SERVER_ERROR_STRING = "500 internal server error";

XmlPicrParser::XmlPicrParser(ResponseType _responseType) :
    responseType(_responseType) {
}

void XmlPicrParser::parse(const QByteArray& data) {
    if (QString(data).contains(SERVER_ERROR_STRING)) {
        setError(QString("Something gone wrong on the server (%1). Please, try later.").arg(SERVER_ERROR_STRING));
        return;
    }

    QDomDocument pDoc;
    pDoc.setContent(data);

    QDomElement docElement = pDoc.documentElement();
    if (docElement.childNodes().isEmpty()) {
        setError("No IDs found.");
        return;
    }

    QDomNodeList returnList;
    QDomElement returnElement;
    if (responseType == Accession) {
        // Find the longest matching sequence
        returnList = docElement.elementsByTagName(PICR_ACCESSION_RETURN);
        int maxLength = 0;
        int maxIndex = 0;
        for (int i = 0; i < returnList.count(); ++i) {
            QDomNodeList sequenceNodeList = returnList.item(i).toElement().elementsByTagName(PICR_SEQUENCE);
            if (1 != sequenceNodeList.count()) {
                setError(QString("Invalid XML: %1 sequences found in the \"%2\" part").
                         arg(sequenceNodeList.count()).
                         arg(PICR_ACCESSION_RETURN));
                return;
            }
            QDomElement sequenceElement = sequenceNodeList.item(0).toElement();

            if (maxLength < sequenceElement.text().length()) {
                maxLength = sequenceElement.text().length();
                maxIndex = i;
            }
        }

        returnElement = returnList.item(maxIndex).toElement();
        if (returnElement.isNull()) {
            setError(QString("Invalid XML: \"%1\" element is not found").arg(PICR_ACCESSION_RETURN));
            return;
        }
    }
    else if (responseType == ExactSequence) {
        // Find the longest matching sequence
        returnList = docElement.elementsByTagName(PICR_SEQUENCE_RETURN);
        int maxLength = 0;
        int maxIndex = 0;
        for (int i = 0; i < returnList.count(); ++i) {
            QDomNodeList sequenceNodeList = returnList.item(i).toElement().elementsByTagName(PICR_SEQUENCE);
            if (1 != sequenceNodeList.count()) {
                setError(QString("Invalid XML: %1 sequences found in the \"%2\" part").
                         arg(sequenceNodeList.count()).
                         arg(PICR_SEQUENCE_RETURN));
                return;
            }
            QDomElement sequenceElement = sequenceNodeList.item(0).toElement();

            if (maxLength < sequenceElement.text().length()) {
                maxLength = sequenceElement.text().length();
                maxIndex = i;
            }
        }

        returnElement = returnList.item(maxIndex).toElement();
        returnElement = returnList.item(maxIndex).toElement();
        if (returnElement.isNull()) {
            setError(QString("Invalid XML: \"%1\" element is not found").arg(PICR_SEQUENCE_RETURN));
            return;
        }
    }
    else if (responseType == Blast) {
        returnElement = docElement;
    }

    QDomNodeList icrList = returnElement.elementsByTagName(PICR_IDENTICAL_CROSS_REFERENCE);
    if (icrList.isEmpty()) {
        setError(QString("No %1 tag").arg(PICR_IDENTICAL_CROSS_REFERENCE));
        return;
    }

    for (int i = 0; i < icrList.count(); ++i) {
        QDomNode icr = icrList.item(i);
        PicrElement element;

        for (int j = 0; j < icr.childNodes().count(); j++) {
            QDomElement childElement = icr.childNodes().item(j).toElement();
            if (childElement.tagName() == PICR_ACCESSION_NUMBER) {
                element.accessionNumber = childElement.text();
            }
            if (childElement.tagName() == PICR_ACCESSION_VERSION) {
                bool ok = false;
                int value = childElement.text().toInt(&ok);
                if (ok) {
                    element.accessionVersion = value;
                }
            }
            if (childElement.tagName() == PICR_DATABASE_DESCRIPTION) {
                element.databaseDescription = childElement.text();
            }
            if (childElement.tagName() == PICR_DATABASE_NAME) {
                element.databaseName = childElement.text();
            }
            if (childElement.tagName() == PICR_NCBI_GI) {
                bool ok = false;
                int value = childElement.text().toInt(&ok);
                if (ok) {
                    element.ncbiId = value;
                }
            }
            if (childElement.tagName() == PICR_TAXON_ID) {
                bool ok = false;
                int value = childElement.text().toInt(&ok);
                if (ok) {
                    element.taxonId = value;
                }
            }
        }

        // Sometimes the "accession" element inculdes accession version. It should be splitted.
        if (element.accessionVersion < 0 &&
                (element.accessionNumber.contains(QRegExp("\\.|\\-")))) {
            QStringList accession = element.accessionNumber.split(QRegExp("\\.|\\-"));
            SAFE_POINT(accession.count() > 1, "Invalid accession field in the XML response", );

            element.accessionNumber = accession.takeFirst();
            bool ok = false;
            int value = accession.first().toInt(&ok);
            if (ok) {
                element.accessionVersion = value;
            }
        }

        if (!results.contains(element)) {
            results << element;
        }
    }

    return;
}

QString XmlPicrParser::getError() const {
    return error;
}

void XmlPicrParser::setError(const QString& val) {
    error = val;
}

QList<PicrElement> XmlPicrParser::getResults() {
    return results;
}

bool XmlPicrParser::isResultsEmpty() {
    return results.isEmpty();
}

//////////////////////////////////////////////////////////////////////////
//CommonDasSettings
const QString CommonDasSettings::dbSwissprot = "SWISSPROT";
const QString CommonDasSettings::dbTrembl = "TREMBL";
const QString CommonDasSettings::dbSwissprotVarsplic = "SWISSPROT_VARSPLIC";
const QString CommonDasSettings::dbTremblVarsplic = "TREMBL_VARSPLIC";
const QString CommonDasSettings::databaseStr = "&database=" + CommonDasSettings::dbSwissprot +
                                               "&database=" + CommonDasSettings::dbTrembl +
                                               "&database=" + CommonDasSettings::dbSwissprotVarsplic +
                                               "&database=" + CommonDasSettings::dbTremblVarsplic;

int CommonDasSettings::getPriority(const QString &databaseName) {
    return getPriorityTable().value(databaseName, -1);
}

QHash<QString, int> CommonDasSettings::getPriorityTable() {
    static QHash<QString, int> priorityTable;
    if (Q_UNLIKELY(priorityTable.isEmpty())) {
        priorityTable.insert(dbTremblVarsplic, 0);
        priorityTable.insert(dbSwissprotVarsplic, 1);
        priorityTable.insert(dbTrembl, 2);
        priorityTable.insert(dbSwissprot, 3);
    }
    return priorityTable;
}

//////////////////////////////////////////////////////////////////////////
//ConvertDasIdTask
const QString ConvertDasIdTask::baseUrl = QString("http://www.ebi.ac.uk/Tools/picr/rest/getUPIForAccession");
const QString ConvertDasIdTask::databasePart = CommonDasSettings::databaseStr;
const QString ConvertDasIdTask::parametersPart = QString("&includeattributes=true");

ConvertDasIdTask::ConvertDasIdTask(const QString& resId) :
    Task(tr("Convert resource id: %1").arg(resId), TaskFlags_FOSCOE | TaskFlag_MinimizeSubtaskErrorText),
    resourceId(resId),
    loop(NULL),
    downloadReply(NULL),
    networkManager(NULL) {
}

ConvertDasIdTask::~ConvertDasIdTask() {
    delete loop;
    delete networkManager;
}

void ConvertDasIdTask::run() {
    if (stateInfo.isCanceled()) {
        return;
    }
    stateInfo.progress = 0;
    ioLog.trace("Sending request to PICR...");

    loop = new QEventLoop;

    networkManager = new QNetworkAccessManager();
    connect(networkManager, SIGNAL(finished(QNetworkReply*)), this, SLOT(sl_replyFinished(QNetworkReply*)));
    connect(networkManager, SIGNAL(proxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)), this, SLOT(onProxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)));
    NetworkConfiguration* nc = AppContext::getAppSettings()->getNetworkConfiguration();

    ioLog.trace("Downloading xml file...");

    QString fetchUrl = getRequestUrlString();
    QNetworkProxy proxy = nc->getProxyByUrl(fetchUrl);
    networkManager->setProxy(proxy);
    ioLog.trace(fetchUrl);

    QUrl requestUrl(fetchUrl);
    downloadReply = networkManager->get(QNetworkRequest(requestUrl));
    connect(downloadReply, SIGNAL(error(QNetworkReply::NetworkError)),
        this, SLOT(sl_onError(QNetworkReply::NetworkError)));
    connect( downloadReply, SIGNAL(uploadProgress( qint64, qint64 )),
        this, SLOT(sl_uploadProgress(qint64,qint64)) );

    QTimer::singleShot(100, this, SLOT(sl_cancelCheck()));
    QTimer::singleShot(60000, this, SLOT(sl_timeout()));

    loop->exec();
    disconnect(0, 0, this, 0);
    if (isCanceled() || hasError()) {
        return;
    }

    ioLog.trace("Download finished.");

    QByteArray result = downloadReply->readAll();

    //parse output
    XmlPicrParser parser(XmlPicrParser::Accession);
    parser.parse(result);
    if (!parser.getError().isEmpty()) {
        setError(parser.getError());
    } else {
        accNumber = getUniprotBestGuessResult(parser.getResults());
    }
}

const QString &ConvertDasIdTask::getSourceAccessionNumber() const {
    return resourceId;
}

const QString &ConvertDasIdTask::getAccessionNumber() const {
    return accNumber;
}

void ConvertDasIdTask::sl_replyFinished(QNetworkReply* reply) {
    Q_UNUSED(reply);
    loop->exit();
}

void ConvertDasIdTask::sl_onError(QNetworkReply::NetworkError error) {
    QNetworkReply *netReply = qobject_cast<QNetworkReply *>(sender());
    QString errorText;
    if (Q_LIKELY(NULL != netReply)) {
        errorText = netReply->errorString();
    } else {
        errorText = tr("undefined error (code %1)").arg(error);
    }
    stateInfo.setError(QString("Network error: %1").arg(errorText));
    loop->exit();
}

void ConvertDasIdTask::sl_uploadProgress(qint64 bytesSent, qint64 bytesTotal) {
    stateInfo.progress = bytesSent/ bytesTotal * 100;
}

void ConvertDasIdTask::sl_cancelCheck() {
    if (isCanceled()) {
        if (loop->isRunning()) {
            loop->exit();
        }
    } else {
        QTimer::singleShot(100, this, SLOT(sl_cancelCheck()));
    }
}

void ConvertDasIdTask::sl_timeout() {
    if (!hasError() && !isCanceled()) {
        setError(tr("Remote server does not respond"));
    }
    if (loop->isRunning()) {
        loop->exit();
    }
}

void ConvertDasIdTask::onProxyAuthenticationRequired(const QNetworkProxy &proxy, QAuthenticator *auth){
    auth->setUser(proxy.user());
    auth->setPassword(proxy.password());
    disconnect(this, SLOT(onProxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)));
}

QString ConvertDasIdTask::getRequestUrlString() {
    QString res = "";
    if (resourceId.isEmpty()) {
        return res;
    }

    res = baseUrl + "?accession=" + resourceId + databasePart + parametersPart;

    return res;
}

QString ConvertDasIdTask::getUniprotBestGuessResult(const QList<PicrElement> &results) {
    // UniProt 'best guess': In order of precedence: Swiss-Prot, TrEMBL, Swiss-Prot varsplic, and TrEMBL varsplic
    CHECK(!results.isEmpty(), "");
    PicrElement bestGuess = results.first();
    foreach (const PicrElement &result, results) {
        if (result > bestGuess) {
            bestGuess = result;
        }
    }
    return bestGuess.accessionNumber;
}


//////////////////////////////////////////////////////////////////////////
//GetDasIdsBySequenceTask
GetDasIdsBySequenceTask::GetDasIdsBySequenceTask(const QByteArray& _sequence) :
    Task(tr("Get an ID for the sequence"), TaskFlags_FOSCOE | TaskFlag_MinimizeSubtaskErrorText),
    sequence(_sequence),
    loop(NULL),
    downloadReply(NULL),
    networkManager(NULL) {
}

GetDasIdsBySequenceTask::~GetDasIdsBySequenceTask() {
    delete loop;
    delete networkManager;
}

void GetDasIdsBySequenceTask::run() {
    if (stateInfo.isCanceled()) {
        return;
    }
    stateInfo.progress = 0;
    ioLog.trace("Sending request to PICR...");

    loop = new QEventLoop;

    networkManager = new QNetworkAccessManager();
    connect(networkManager, SIGNAL(finished(QNetworkReply*)), this, SLOT(sl_replyFinished(QNetworkReply*)));
    connect(networkManager, SIGNAL(proxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)), this, SLOT(onProxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)));
    NetworkConfiguration* nc = AppContext::getAppSettings()->getNetworkConfiguration();

    ioLog.trace("Downloading xml file...");

    QString fetchUrl = getRequestUrlString();
    SAFE_POINT(!fetchUrl.isEmpty(), QString("Invalid fetch URL: %1").arg(fetchUrl), );
    QNetworkProxy proxy = nc->getProxyByUrl(fetchUrl);
    networkManager->setProxy(proxy);
    ioLog.trace(fetchUrl);

    QUrl requestUrl(fetchUrl);
    downloadReply = networkManager->get(QNetworkRequest(requestUrl));
    connect(downloadReply, SIGNAL(error(QNetworkReply::NetworkError)),
        this, SLOT(sl_onError(QNetworkReply::NetworkError)));
    connect( downloadReply, SIGNAL(uploadProgress( qint64, qint64 )),
        this, SLOT(sl_uploadProgress(qint64,qint64)) );

    loop->exec();
    ioLog.trace("Download finished.");

    QByteArray result = downloadReply->readAll();

    //parse output
    XmlPicrParser parser = getParser();
    parser.parse(result);
    if (!parser.getError().isEmpty()) {
        setError(parser.getError());
    } else {
        results = parser.getResults();
    }
}

QList<PicrElement> GetDasIdsBySequenceTask::getResults() {
    return results;
}

void GetDasIdsBySequenceTask::sl_replyFinished(QNetworkReply* reply) {
    Q_UNUSED(reply);
    loop->exit();
}

void GetDasIdsBySequenceTask::sl_onError(QNetworkReply::NetworkError error) {
    QNetworkReply *netReply = qobject_cast<QNetworkReply *>(sender());
    QString errorText;
    if (Q_LIKELY(NULL != netReply)) {
        errorText = netReply->errorString();
    } else {
        errorText = tr("undefined error (code %1)").arg(error);
    }
    stateInfo.setError(QString("Network error: %1").arg(errorText));
    loop->exit();
}

void GetDasIdsBySequenceTask::sl_uploadProgress(qint64 bytesSent, qint64 bytesTotal) {
    stateInfo.progress = bytesSent/ bytesTotal * 100;
}

void GetDasIdsBySequenceTask::onProxyAuthenticationRequired(const QNetworkProxy &proxy, QAuthenticator *auth){
    auth->setUser(proxy.user());
    auth->setPassword(proxy.password());
    disconnect(this, SLOT(onProxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)));
}

//////////////////////////////////////////////////////////////////////////
//GetDasIdsByExactSequenceTask
const QString GetDasIdsByExactSequenceTask::baseUrl = QString("http://www.ebi.ac.uk/Tools/picr/rest/getUPIForSequence");
const QString GetDasIdsByExactSequenceTask::databasePart = CommonDasSettings::databaseStr;
const QString GetDasIdsByExactSequenceTask::parametersPart = QString("&includeattributes=true");

GetDasIdsByExactSequenceTask::GetDasIdsByExactSequenceTask(const QByteArray& _sequence) :
    GetDasIdsBySequenceTask(_sequence) {
}

GetDasIdsByExactSequenceTask::~GetDasIdsByExactSequenceTask() {
}

QString GetDasIdsByExactSequenceTask::getRequestUrlString() {
    QString res;
    if (sequence.isEmpty()) {
        return res;
    }

    res = baseUrl + "?sequence=" + sequence + databasePart + parametersPart;
    return res;
}

XmlPicrParser GetDasIdsByExactSequenceTask::getParser() {
    return XmlPicrParser(XmlPicrParser::ExactSequence);
}

//////////////////////////////////////////////////////////////////////////
//PicrBlastSettings
const QString PicrBlastSettings::PROGRAM = "program";
const QString PicrBlastSettings::IDENTITY = "identity";
const QString PicrBlastSettings::MATRIX = "matrix";
const QString PicrBlastSettings::FILTER = "filter";
const QString PicrBlastSettings::GAP_OPEN = "gap_open";
const QString PicrBlastSettings::GAP_EXT = "gap_ext";
const QString PicrBlastSettings::DROP_OFF = "drop_off";
const QString PicrBlastSettings::GAP_ALIGN = "gap_align";

void PicrBlastSettings::insert(const QString& key, const QVariant& value) {
    settings.insert(key, value);
}

bool PicrBlastSettings::isValid() {
    bool allEntries = !settings.value(PROGRAM, QString()).toString().isEmpty() &&
                      !settings.value(IDENTITY, QString()).toString().isEmpty() &&
                      !settings.value(MATRIX, QString()).toString().isEmpty() &&
                      !settings.value(FILTER, QString()).toString().isEmpty() &&
                      !settings.value(GAP_OPEN, QString()).toString().isEmpty() &&
                      !settings.value(GAP_EXT, QString()).toString().isEmpty() &&
                      !settings.value(DROP_OFF, QString()).toString().isEmpty() &&
                      !settings.value(GAP_ALIGN, QString()).toString().isEmpty();
    if (false == allEntries) {
        return false;
    }

    bool ok;
    settings.value(GAP_OPEN).toInt(&ok);
    if (false == ok) {
        return false;
    }
    settings.value(GAP_EXT).toInt(&ok);
    if (false == ok) {
        return false;
    }

    return true;
}

QString PicrBlastSettings::value(const QString& key) {
    return settings.value(key).toString();
}


//////////////////////////////////////////////////////////////////////////
//GetDasIdsByBlastTask
const QString GetDasIdsByBlastTask::baseUrl = QString("http://www.ebi.ac.uk/Tools/picr/rest/getUPIForBLAST");
const QString GetDasIdsByBlastTask::databasePart = CommonDasSettings::databaseStr;
const QString GetDasIdsByBlastTask::parametersPart = QString("&program=%1&filtertype=IDENTITY&identityvalue=%2&includeattributes=true&matrix=%3&filter=%4&gapopen=%5&gapext=%6&dropoff=%7&gapalign=%8");

GetDasIdsByBlastTask::GetDasIdsByBlastTask(const QByteArray& _sequence, const PicrBlastSettings& _settings) :
    GetDasIdsBySequenceTask(_sequence),
    settings(_settings) {
}

GetDasIdsByBlastTask::~GetDasIdsByBlastTask() {
}

QString GetDasIdsByBlastTask::getRequestUrlString() {
    QString res;
    if (sequence.isEmpty() || !settings.isValid()) {
        return res;
    }

    res = baseUrl + "?blastfrag=" + sequence + databasePart + parametersPart.arg(settings.value(PicrBlastSettings::PROGRAM)).
          arg(settings.value(PicrBlastSettings::IDENTITY)).
          arg(settings.value(PicrBlastSettings::MATRIX)).
          arg(settings.value(PicrBlastSettings::FILTER)).
          arg(settings.value(PicrBlastSettings::GAP_OPEN)).
          arg(settings.value(PicrBlastSettings::GAP_EXT)).
          arg(settings.value(PicrBlastSettings::DROP_OFF)).
          arg(settings.value(PicrBlastSettings::GAP_ALIGN));

    return res;
}

XmlPicrParser GetDasIdsByBlastTask::getParser() {
    return XmlPicrParser(XmlPicrParser::Blast);
}

}   // namespace
