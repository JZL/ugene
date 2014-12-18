/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2014 UniPro <ugene@unipro.ru>
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

#ifndef _U2_SNP_EFF_TASK_H_
#define _U2_SNP_EFF_TASK_H_

#include <U2Core/Task.h>
#include <U2Core/ExternalToolRunTask.h>

namespace U2 {


class SnpEffSetting{
public:
    SnpEffSetting(): inputUrl(""), outDir(""), inFormat(""), outFormat(""),genome(""), updownLength(""),homohetero(""), seqChange(""), filterOut(""), chrPos(""){}

    QString inputUrl;
    QString outDir;
    QString inFormat;
    QString outFormat;
    QString genome;
    QString updownLength;
    QString homohetero;
    QString seqChange;
    QString filterOut;
    QString chrPos;
};

class SnpEffTask : public Task {
    Q_OBJECT
public:
    SnpEffTask(const SnpEffSetting &settings);

    void prepare();
    void run();

    QString getResult(){return resultUrl;}
    QString getSummaryUrl();
    QString getResFileUrl();

protected:
    QString getDataPath();
    QStringList getParameters(U2OpStatus& os);

protected:
    SnpEffSetting settings;
    QString resultUrl;
};


class SnpEffParser : public ExternalToolLogParser {
public:
    SnpEffParser();

    void parseOutput(const QString& partOfLog);
    void parseErrOutput(const QString& partOfLog);

private:
    QString lastErrLine;
};

}//namespace

#endif // _U2_SNP_EFF_H_

