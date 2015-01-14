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
#include <QtCore/QUrl>

#include <QtNetwork/QAuthenticator>
#include <QtNetwork/QNetworkReply>

#include <QtXml/QDomDocument>

#include <U2Core/AnnotationTableObject.h>
#include <U2Core/AppContext.h>
#include <U2Core/AppSettings.h>
#include <U2Core/DNAAlphabet.h>
#include <U2Core/DNASequence.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/DocumentModel.h>
#include <U2Core/GObjectRelationRoles.h>
#include <U2Core/LoadDocumentTask.h>
#include <U2Core/MultiTask.h>
#include <U2Core/PicrApiTask.h>
#include <U2Core/SaveDocumentTask.h>
#include <U2Core/U2AlphabetUtils.h>
#include <U2Core/U2DbiRegistry.h>
#include <U2Core/U2SequenceUtils.h>

#include "LoadDASDocumentTask.h"

namespace U2 {

LoadDasDocumentTask::LoadDasDocumentTask( const QString& accId, const QString& _fullPath, const DASSource& _referenceSource, const QList<DASSource>& _featureSources )
: BaseLoadRemoteDocumentTask(_fullPath, QVariantMap(), TaskFlags(TaskFlag_NoRun | TaskFlag_MinimizeSubtaskErrorText))
,accNumber(accId)
,featureSources(_featureSources)
,referenceSource(_referenceSource)
,loadSequenceTask(NULL)
,saveDocumentTask(NULL)
,seq(NULL)
{

}

void LoadDasDocumentTask::prepare(){
    BaseLoadRemoteDocumentTask::prepare();
    if (!isCached()){
        //load sequence
        loadSequenceTask = new LoadDasObjectTask(accNumber, referenceSource, DASSequence);
        addSubTask(loadSequenceTask);

        //load annotations
        foreach(const DASSource& s, featureSources){
            LoadDasObjectTask* featureTask = new LoadDasObjectTask(accNumber, s, DASFeatures);
            loadFeaturesTasks.append(featureTask);
            addSubTask(featureTask);
        }
    }
}

QString LoadDasDocumentTask::getFileFormat( const QString & dbid ){
    Q_UNUSED(dbid);
    return GENBANK_FORMAT;
}

GUrl LoadDasDocumentTask::getSourceUrl(){
    return GUrl();
}

QString LoadDasDocumentTask::getFileName(){
    format = getFileFormat("");
    accNumber.replace(";",",");
    QStringList accIds = accNumber.split(",");
    if (accIds.size() == 1 ) {
        return accNumber + "_das" +"." + format;
    } else if (accIds.size() > 1) {
        return accIds.first() + "_das_misc." + format;
    }

    return "";
}

QList<Task*> LoadDasDocumentTask::onSubTaskFinished( Task *subTask ){
    QList<Task *> subTasks;
    if ( isCanceled() || hasError() ) {
        return subTasks;
    }

    if ( subTask == loadDocumentTask ) {
        if ( subTask->hasError( ) ) {
            setError( tr( "Cannot load cached document: %1" ).arg( accNumber ) );
            return subTasks;
        }
        resultDocument = loadDocumentTask->takeDocument( );
    } else if ( subTask == saveDocumentTask ) {
        if ( saveDocumentTask->hasError( ) ) {
            setError( tr( "Cannot save document: %1" ).arg( accNumber ) );
            return subTasks;
        }
        if ( !subTask->isCanceled( ) ) {
            RecentlyDownloadedCache * cache = AppContext::getRecentlyDownloadedCache( );
            if ( NULL != cache ) {
                cache->append( fullPath );
            }
        }
    } else {
        if ( subTask == loadSequenceTask ) {
            if ( loadSequenceTask->hasError( ) ) {
                setError( tr( "Cannot find DAS reference sequence: %1" ).arg( accNumber ) );
                return subTasks;
            }
            if ( !isCanceled( ) ) {
                seq = loadSequenceTask->getSequence( );
            }

            loadSequenceTask = NULL;
        } else {
            LoadDasObjectTask *ftask = qobject_cast<LoadDasObjectTask *>( subTask );
            if ( NULL != ftask && loadFeaturesTasks.contains( ftask ) ) {
                const int idx = loadFeaturesTasks.indexOf( ftask );
                if ( idx == -1 ) {
                    return subTasks;
                }

                loadFeaturesTasks.removeAt( idx );

                if ( ftask->hasError( ) ) {
                    ioLog.info( tr( "Cannot find DAS features for '%1' on %2" ).arg( accNumber ).arg( ftask->getSource( ).getName( ) ) );
                } else {
                    //merge features
                    if ( !isCanceled( ) ) {
                        mergeFeatures( ftask->getAnnotationData( ) );
                    }
                }
            }
        }

        if ( isAllDataLoaded( ) ) {
            AnnotationTableObject *annotationTableObject = NULL;
            if ( !annotationData.isEmpty( ) ) {
                U2DbiRegistry *dbiReg = AppContext::getDbiRegistry();
                SAFE_POINT_EXT(NULL != dbiReg, setError("NULL DBI registry"), subTasks);

                U2DbiRef dbiRef = dbiReg->getSessionTmpDbiRef(stateInfo);
                CHECK_OP(stateInfo, subTasks);

                annotationTableObject = new AnnotationTableObject( "das_annotations", dbiRef );

                foreach ( const QString &grname, annotationData.keys( ) ) {
                    const QList<AnnotationData> sdata = annotationData[grname];
                    if ( !sdata.isEmpty( ) ) {
                        foreach ( AnnotationData d, sdata ) {
                            //setRegion
                            if ( NULL != seq ) {
                                if ( d.location->isSingleRegion( )
                                    && d.location->regions.first() == U2_REGION_MAX)
                                {
                                    U2Location newLoc = d.location;
                                    newLoc->regions.clear( );
                                    newLoc->regions.append( U2Region( 0, seq->length( ) ) );
                                    d.location = newLoc;
                                }
                            }
                            annotationTableObject->addAnnotation( d, grname );
                        }
                    }
                }
            }
            if ( NULL != seq ) {
                createLoadedDocument( );
                if ( NULL == resultDocument ) {
                    return subTasks;
                }

                U2EntityRef seqRef = U2SequenceUtils::import( resultDocument->getDbiRef( ), *seq, stateInfo );
                if ( stateInfo.isCoR( ) ) {
                    return subTasks;
                }
                U2SequenceObject *danseqob = new U2SequenceObject( seq->getName( ), seqRef );
                resultDocument->addObject( danseqob );

                if ( NULL != annotationTableObject ) {
                    annotationTableObject->addObjectRelation( GObjectRelation( danseqob,
                        ObjectRole_Sequence ) );
                    resultDocument->addObject( annotationTableObject );

                }

                saveDocumentTask = new SaveDocumentTask( resultDocument );
                subTasks.append( saveDocumentTask );
            }
         }
    }
    return subTasks;

}

bool LoadDasDocumentTask::isAllDataLoaded( ) {
    return ( NULL == loadSequenceTask && loadFeaturesTasks.isEmpty( ) );
}

void LoadDasDocumentTask::mergeFeatures( const QMap<QString, QList<AnnotationData> > &newAnnotations ) {
    const QStringList &keys =  newAnnotations.keys( );
    foreach ( const QString &key, keys ) {
        if ( annotationData.contains(key ) ) {
            const QList<AnnotationData> &curList = annotationData[key];
            const QList<AnnotationData> &tomergeList = newAnnotations[key];
            foreach ( const AnnotationData &d, tomergeList ) {
                if ( !curList.contains(d ) ) {
                    annotationData[key].append( d );
                }
            }
        } else {
            annotationData.insert( key, newAnnotations[key] );
        }
    }

}



//////////////////////////////////////////////////////////////////////////
//LoadDASObjectTask
LoadDasObjectTask::LoadDasObjectTask( const QString& accId, const DASSource& _source, DASObjectType objType)
    :Task(tr("Load DAS data for '%1' from %2").arg(accId).arg(_source.getName()), TaskFlags_FOSCOE | TaskFlag_MinimizeSubtaskErrorText)
,accNumber(accId)
,source(_source)
,objectType(objType)
,loop(NULL)
,downloadReply(NULL)
,networkManager(NULL)
,seq(NULL)
{

}

LoadDasObjectTask::~LoadDasObjectTask(){
    delete loop;
    delete networkManager;
}

void LoadDasObjectTask::run(){
    if (stateInfo.isCanceled()){
        return;
    }
    stateInfo.progress = 0;
    ioLog.trace("Start loading data from DAS...");

    loop = new QEventLoop;

    networkManager = new QNetworkAccessManager();
    connect(networkManager, SIGNAL(finished(QNetworkReply*)), this, SLOT(sl_replyFinished(QNetworkReply*)));
    connect(networkManager, SIGNAL(proxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)), this, SLOT(onProxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)));
    NetworkConfiguration* nc = AppContext::getAppSettings()->getNetworkConfiguration();

    ioLog.trace("Downloading xml file...");

    QString fetchUrl = DASSourceRegistry::getRequestURLString(source, accNumber, objectType);
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
    if ( ( result.size() < 100 ) && result.contains("Nothing has been found")) {
        setError(tr("Sequence with ID=%1 is not found.").arg(accNumber));
        return;
    }

    //parse output
    if (objectType == DASSequence){
        XMLDASSequenceParser parser;
        parser.parse(result);
        if (!parser.getError().isEmpty()){
            setError(parser.getError());
        }else{
            seq = parser.getSequence();
        }
    }else if(objectType == DASFeatures){
        XMLDASFeaturesParser parser;
        parser.parse(result);
        if (!parser.getError().isEmpty()){
            setError(parser.getError());
        }else{
            annotationData = parser.getAnnotationData();
        }
    }

}

DNASequence* LoadDasObjectTask::getSequence() {
    return seq;
}

const QString& LoadDasObjectTask::getAccession() const {
    return accNumber;
}

const DASSource& LoadDasObjectTask::getSource() const {
    return source;
}

const QMap<QString, QList<AnnotationData> >& LoadDasObjectTask::getAnnotationData( ) const {
    return annotationData;
}

void LoadDasObjectTask::sl_replyFinished( QNetworkReply* reply ) {
    Q_UNUSED(reply);
    loop->exit();
}

void LoadDasObjectTask::sl_onError( QNetworkReply::NetworkError error ){
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

void LoadDasObjectTask::sl_uploadProgress( qint64 bytesSent, qint64 bytesTotal ){
    stateInfo.progress = bytesSent/ bytesTotal * 100;
}

void LoadDasObjectTask::sl_cancelCheck() {
    if (isCanceled()) {
        if (loop->isRunning()) {
            loop->exit();
        }
    } else {
        QTimer::singleShot(100, this, SLOT(sl_cancelCheck()));
    }
}

void LoadDasObjectTask::sl_timeout() {
    if (!hasError() && !isCanceled()) {
        setError(tr("Remote server does not respond"));
    }
    if (loop->isRunning()) {
        loop->exit();
    }
}

void LoadDasObjectTask::onProxyAuthenticationRequired(const QNetworkProxy &proxy, QAuthenticator *auth){
    auth->setUser(proxy.user());
    auth->setPassword(proxy.password());
    disconnect(this, SLOT(onProxyAuthenticationRequired(const QNetworkProxy&, QAuthenticator*)));
}

//////////////////////////////////////////////////////////////////////////
//LoadDasFeaturesTask
LoadDasFeaturesTask::LoadDasFeaturesTask(const QStringList& accId, const QList<DASSource>& source)
: Task(tr("Load DAS annotations for current sequence"), TaskFlags(TaskFlag_NoRun) | TaskFlag_CancelOnSubtaskCancel | TaskFlag_ReportingIsSupported | TaskFlag_ReportingIsEnabled ), featureSources(source), accessionNumbers(accId) {
}

const QMap<QString, QList<AnnotationData> >& LoadDasFeaturesTask::getAnnotationData( ) const {
    return annotationData;
}

void LoadDasFeaturesTask::prepare() {
    foreach(const QString &accNumber, accessionNumbers){
        foreach (DASSource featureSource, featureSources) {
            LoadDasObjectTask * loadAnnotationsTask = new LoadDasObjectTask(accNumber, featureSource, DASFeatures);
            addSubTask(loadAnnotationsTask);
        }
    }
}

QList<Task*> LoadDasFeaturesTask::onSubTaskFinished(Task* subTask) {
    QList<Task*> res;
    LoadDasObjectTask* loadDasObjectTask = dynamic_cast<LoadDasObjectTask*>(subTask);
    SAFE_POINT(NULL != loadDasObjectTask, "Incorrect subtask in LoadDasObjectsTask", res);
    if(loadDasObjectTask->hasError()) {
        reports += "<font size=\"5\" color=\"orange\">";
        reports += tr("Can not receive response from the server \"") + loadDasObjectTask->getSource().getName() + "\"</font><br>";
    }
    else {
        QMap<QString, QList<AnnotationData> > data = loadDasObjectTask->getAnnotationData();
        int annotationsNumber = 0;
        foreach(const QString& key, data.keys()) {
            annotationsNumber += data[key].size();
        }
        reports += tr("<font size=\"5\" color=\"green\">Received %1 annotations from the server \"%2\"</font><br>").arg(annotationsNumber).arg(loadDasObjectTask->getSource().getName());
        mergeFeatures(loadDasObjectTask->getAnnotationData());
    }
    return res;
}

void LoadDasFeaturesTask::mergeFeatures(const QMap<QString, QList<AnnotationData> >& newAnnotations) {
    const QStringList& keys =  newAnnotations.keys();
    foreach (const QString& key, keys) {
        if (annotationData.contains(key)) {
            const QList<AnnotationData>& curList = annotationData[key];
            const QList<AnnotationData>& tomergeList = newAnnotations[key];
            foreach ( const AnnotationData &d, tomergeList ) {
                if (!curList.contains(d)) {
                    annotationData[key].append(d);
                }
            }
        } else {
            annotationData.insert(key, newAnnotations[key]);
        }
    }
}

QString LoadDasFeaturesTask::generateReport() const {
    return reports;
}

//////////////////////////////////////////////////////////////////////////
//ConvertAndLoadDASDocumentTask
ConvertIdAndLoadDasDocumentTask::ConvertIdAndLoadDasDocumentTask(const QString& accId,
                                                                 const QString& _fullPath,
                                                                 const DASSource& _referenceSource,
                                                                 const QList<DASSource>& _featureSources,
                                                                 bool _convertId) :
    Task(QString("Convert ID and load DAS document for: %1").arg(accId), TaskFlags(TaskFlag_NoRun | TaskFlag_CancelOnSubtaskCancel | TaskFlag_MinimizeSubtaskErrorText)),
    convertDasIdTask(NULL),
    loadDasDocumentTask(NULL),
    accessionNumber(accId),
    fullPath(_fullPath),
    referenceSource(_referenceSource),
    featureSources(_featureSources),
    convertId(_convertId)
{

}

QString ConvertIdAndLoadDasDocumentTask::getConvertedAccessionNumber() const {
    CHECK(NULL != convertDasIdTask && convertDasIdTask->isFinished(), "");
    CHECK(!convertDasIdTask->isCanceled() && !convertDasIdTask->hasError(), "");
    return convertDasIdTask->getAccessionNumber();
}

void ConvertIdAndLoadDasDocumentTask::prepare() {
    if (convertId) {
        convertDasIdTask = new ConvertDasIdTask(accessionNumber);
        addSubTask(convertDasIdTask);
    } else {
        loadDasDocumentTask = new LoadDasDocumentTask(accessionNumber, fullPath, referenceSource, featureSources);
        addSubTask(loadDasDocumentTask);
    }
}

QList<Task*> ConvertIdAndLoadDasDocumentTask::onSubTaskFinished(Task *subTask) {
    QList<Task*> subTasks;

    if (subTask->isCanceled()) {
        return subTasks;
    }

    if (subTask == convertDasIdTask) {
        if (!convertDasIdTask->getAccessionNumber().isEmpty() && !convertDasIdTask->hasError()) {
            ioLog.details(QString("\"%1\" was converted into \"%2\"").
                          arg(accessionNumber).
                          arg(convertDasIdTask->getAccessionNumber()));
            accessionNumber = convertDasIdTask->getAccessionNumber();
        }
        loadDasDocumentTask = new LoadDasDocumentTask(accessionNumber, fullPath, referenceSource, featureSources);
        subTasks << loadDasDocumentTask;
    }
    if (subTask == loadDasDocumentTask && loadDasDocumentTask->hasError()) {
        setError(loadDasDocumentTask->getError());
    }

    return subTasks;
}

//////////////////////////////////////////////////////////////////////////
//XMLDASSequenceParser
XMLDASSequenceParser::XMLDASSequenceParser(){
    seq = NULL;
}

#define DAS_SEQ_DASSEQUENCE "DASSEQUENCE"
#define DAS_SEQ_SEQUENCE "SEQUENCE"
#define DAS_SEQ_ID "id"
void XMLDASSequenceParser::parse( const QByteArray& data ){
    //http://www.biodas.org/documents/spec-1.6.html
    QDomDocument pDoc;
    pDoc.setContent(data);

    QDomElement dasSeq = pDoc.documentElement();
    if(dasSeq.tagName() != DAS_SEQ_DASSEQUENCE){
        setError(QString("No %1 tag").arg(DAS_SEQ_DASSEQUENCE));
        return;
    }
    //here may be multiple sequence tags, but we take only the first one
    QDomNode seqNode = dasSeq.firstChild();
    if (!seqNode.isNull()){
        QDomElement seqElement = seqNode.toElement();
        if(seqElement.tagName() != DAS_SEQ_SEQUENCE){
            setError(QString("No %1 tag").arg(DAS_SEQ_SEQUENCE));
            return;
        }
        QString sequenceId = seqElement.attribute(DAS_SEQ_ID).trimmed();
        QByteArray sequence = seqElement.text().toLatin1().trimmed();

        const DNAAlphabet* a = U2AlphabetUtils::findBestAlphabet(sequence.data(), sequence.size());

        seq = new DNASequence(sequenceId, sequence, a);

    }else{
        setError(QString("No %1 tag").arg(DAS_SEQ_SEQUENCE));
        return;
    }

}

//////////////////////////////////////////////////////////////////////////
//XMLDASFeaturesParser
XMLDASFeaturesParser::XMLDASFeaturesParser( ) {

}

#define DAS_FEATURE_DASGFF "DASGFF"
#define DAS_FEATURE_GFF "GFF"
#define DAS_FEATURE_SEGMENT "SEGMENT"
#define DAS_FEATURE_ID "id"
#define DAS_FEATURE_START "start"
#define DAS_FEATURE_STOP  "stop"
#define DAS_FEATURE_FEATURE "FEATURE"
#define DAS_FEATURE_LABEL "label"
#define DAS_FEATURE_HREF "href"

#define DAS_FEATURE_TYPE "TYPE"
#define DAS_FEATURE_METHOD "METHOD"
#define DAS_FEATURE_START_POS "START"
#define DAS_FEATURE_END_POS "END"
#define DAS_FEATURE_SCORE "SCORE"
#define DAS_FEATURE_ORIENTATION "ORIENTATION"
#define DAS_FEATURE_PHASE "PHASE"
#define DAS_FEATURE_NOTE "NOTE"
#define DAS_FEATURE_LINK "LINK"
#define DAS_FEATURE_TARGET "TARGET"
#define DAS_FEATURE_PARENT "PARENT"
#define DAS_FEATURE_PART "PART"

void XMLDASFeaturesParser::parse( const QByteArray& data ){
    //http://www.biodas.org/documents/spec-1.6.html
    QDomDocument pDoc;
    pDoc.setContent(data);

    //DASGFF
    QDomElement dasGff = pDoc.documentElement();
    if(dasGff.tagName() != DAS_FEATURE_DASGFF){
        setError(QString("No %1 tag").arg(DAS_FEATURE_DASGFF));
        return;
    }

    //GFF
    QDomNode gff = dasGff.firstChild();
    if (!gff.isNull()){
        QDomElement gffElement = gff.toElement();
        if(gffElement.tagName() != DAS_FEATURE_GFF){
            setError(QString("No %1 tag").arg(DAS_FEATURE_GFF));
            return;
        }

        //SEGMENT
        QDomNode gffSegment = gffElement.firstChild();
        while (!gffSegment.isNull()){
            QDomElement gffSegmentElement = gffSegment.toElement();
            if(gffSegmentElement.tagName() != DAS_FEATURE_SEGMENT){
                //no annotations
                return;
            }
            QString sequenceId = gffSegmentElement.attribute(DAS_FEATURE_ID);
            qint64 start = gffSegmentElement.attribute(DAS_FEATURE_START).toInt();
            qint64 stop = gffSegmentElement.attribute(DAS_FEATURE_STOP).toInt();

            //FEATURE
            QDomNode featureSegment = gffSegmentElement.firstChild();
            while (!featureSegment.isNull()){
                QDomElement featureSegmentElement = featureSegment.toElement();
                if(featureSegmentElement.tagName() != DAS_FEATURE_FEATURE){
                    setError(QString("No %1 tag").arg(DAS_FEATURE_FEATURE));
                    return;
                }
                //annotation data
                QString featureLabel = featureSegmentElement.attribute(DAS_FEATURE_LABEL);
                QString featureId = featureSegmentElement.attribute(DAS_FEATURE_ID);

                QString groupName = "";
                QString groupId = "";
                QString methodQual = "";
                qint64 startPos = -1;
                qint64 endPos = -1;
                float score = -1.0f;
                bool complemented = false;
                QString note = "";
                QString link = "";
                QString target = "";

                //FEATURE_ATTRS
                QDomNode featureAttrSegment = featureSegmentElement.firstChild();
                while (!featureAttrSegment.isNull()){
                    QDomElement featureAttrElement = featureAttrSegment.toElement();

                    QString tagName = featureAttrElement.tagName();
                    if (tagName == DAS_FEATURE_TYPE){
                        //group name
                        groupName = featureAttrElement.text();
                        groupId = featureAttrElement.attribute(DAS_FEATURE_ID);
                        if (groupName.isEmpty()){
                            groupName = groupId;
                        }
                    }else if (tagName == DAS_FEATURE_METHOD){
                        methodQual = featureAttrElement.text();
                        if (methodQual.isEmpty()){
                            methodQual = featureAttrElement.attribute(DAS_FEATURE_ID);
                        }

                    }else if (tagName == DAS_FEATURE_START_POS){
                        QString startText = featureAttrElement.text();
                        if (!startText.isEmpty()){
                            startPos = startText.toInt();
                        }

                    }else if (tagName == DAS_FEATURE_END_POS){
                        QString endText = featureAttrElement.text();
                        if (!endText.isEmpty()){
                            endPos = endText.toInt();
                        }
                    }else if (tagName == DAS_FEATURE_SCORE){
                        QString scoreText = featureAttrElement.text();
                        if (!scoreText.isEmpty() && scoreText != "-"){
                            score = scoreText.toFloat();
                        }

                    }else if (tagName == DAS_FEATURE_ORIENTATION){
                        QString complText = featureAttrElement.text();
                        if (complText == "-"){
                            complemented = true;
                        }

                    }else if (tagName == DAS_FEATURE_PHASE){
                        //skip

                    }else if (tagName == DAS_FEATURE_NOTE){
                        note = featureAttrElement.text();

                    }else if (tagName == DAS_FEATURE_LINK){
                        link = featureAttrElement.attribute(DAS_FEATURE_HREF);

                    }else if (tagName == DAS_FEATURE_TARGET){
                        target = featureAttrElement.text();
                        if (target.isEmpty()){
                            target = featureAttrElement.attribute(DAS_FEATURE_ID);
                        }

                    }else if (tagName == DAS_FEATURE_PARENT){
                        //skip

                    }else if (tagName == DAS_FEATURE_PART){
                        //skip

                    }


                    featureAttrSegment = featureAttrSegment.nextSibling();
                }

                featureSegment = featureSegment.nextSibling();

                AnnotationData data;
                data.name = featureId.simplified( );
                U2Region reg; //1-based start


                if (startPos == -1 || endPos == -1){//non-positional
                    data.qualifiers.append(U2Qualifier("non-positional", "yes"));
                    if ((start == 0 && stop == 0) || (start > stop)){
                        reg = U2_REGION_MAX;
                    }else{
                        reg = U2Region(start - 1, stop - start + 1);
                    }
                }else{
                    if ((startPos == 0 && endPos == 0) || (startPos > endPos)){
                        reg = U2_REGION_MAX;
                    }else{
                        reg = U2Region(startPos - 1, endPos - startPos + 1);
                    }
                }
                data.location->regions << reg;
                data.setStrand(complemented ? U2Strand::Complementary : U2Strand::Direct);

                if (!methodQual.isEmpty()){
                    data.qualifiers.append(U2Qualifier("method", methodQual.simplified()));
                }

                if (score != -1.0f){
                    data.qualifiers.append(U2Qualifier("score", QString::number(score)));
                }

                if (!note.isEmpty()){
                    data.qualifiers.append(U2Qualifier("note", note.simplified()));
                }

                if (!link.isEmpty()){
                    data.qualifiers.append(U2Qualifier("link", link.simplified()));
                }

                if (!target.isEmpty()){
                    data.qualifiers.append(U2Qualifier("target", target.simplified()));
                }

                if (!sequenceId.isEmpty()){
                    data.qualifiers.append(U2Qualifier("seq_id", sequenceId.simplified()));
                }

                if (groupName.isEmpty()){
                    groupName = "das_features";
                }

                data.qualifiers.append(U2Qualifier("type", groupName.simplified()));

                if (!groupId.isEmpty()){
                    data.qualifiers.append(U2Qualifier("type_id", groupId.simplified()));
                }

                data.qualifiers.append(U2Qualifier("feature_id", featureId.simplified()));
                if (!featureLabel.isEmpty()){
                    data.qualifiers.append(U2Qualifier("feature_label", featureLabel.simplified()));
                }

                annotationData[groupName].append(data);
            }


            gffSegment = gffSegment.nextSibling();
        }

    }else{
        setError(QString("No %1 tag").arg(DAS_FEATURE_GFF));
        return;
    }
}

ConvertIdAndLoadDasFeaturesTask::ConvertIdAndLoadDasFeaturesTask(const QStringList &accessionNumbers, const QList<DASSource> &featureSources, bool convertId) :
    Task(tr("Convert ID and load DAS features for: %1").arg(accessionNumbers.join(", ")), TaskFlags(TaskFlag_CancelOnSubtaskCancel | TaskFlag_NoRun | TaskFlag_MinimizeSubtaskErrorText)),
    convertDasIdTasks(NULL),
    loadDasFeaturesTask(NULL),
    accessionNumbers(accessionNumbers),
    featureSources(featureSources),
    convertId(convertId)
{
}

void ConvertIdAndLoadDasFeaturesTask::prepare() {
    if (convertId) {
        QList<Task *> convertTasks;
        foreach (const QString &accessionNumber, accessionNumbers) {
            convertTasks << new ConvertDasIdTask(accessionNumber);
        }

        convertDasIdTasks = new MultiTask(tr("Convert IDs task"), convertTasks, TaskFlags(TaskFlag_CancelOnSubtaskCancel));
        addSubTask(convertDasIdTasks);
    } else {
        loadDasFeaturesTask = new LoadDasFeaturesTask(accessionNumbers, featureSources);
        addSubTask(loadDasFeaturesTask);
    }
}

QList<Task *> ConvertIdAndLoadDasFeaturesTask::onSubTaskFinished(Task *subTask) {
    QList<Task *> subTasks;

    if (subTask->isCanceled()) {
        return subTasks;
    }

    if (subTask == convertDasIdTasks) {
        QStringList convertedAccessionNumbers;
        foreach (Task *convertTask, convertDasIdTasks->getSubtasks()) {
            ConvertDasIdTask *convertDasIdTask = qobject_cast<ConvertDasIdTask *>(convertTask);
            if (!convertDasIdTask->getAccessionNumber().isEmpty() && !convertDasIdTask->hasError()) {
                ioLog.details(tr("\"%1\" was converted into \"%2\"").
                              arg(convertDasIdTask->getSourceAccessionNumber()).
                              arg(convertDasIdTask->getAccessionNumber()));
                convertedAccessionNumbers << convertDasIdTask->getAccessionNumber();
            }
        }
        loadDasFeaturesTask = new LoadDasFeaturesTask(convertedAccessionNumbers, featureSources);
        subTasks << loadDasFeaturesTask;
    }

    if (subTask == loadDasFeaturesTask && loadDasFeaturesTask->hasError()) {
        setError(loadDasFeaturesTask->getError());
    }

    return subTasks;
}

const QMap<QString, QList<AnnotationData> > &ConvertIdAndLoadDasFeaturesTask::getAnnotationData() const {
    return loadDasFeaturesTask->getAnnotationData();
}

} //namespace
