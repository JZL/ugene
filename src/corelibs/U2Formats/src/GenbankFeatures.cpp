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

#include "GenbankFeatures.h"

#include <QtCore/QHash>
#include <QtCore/QMutex>
#include <QtCore/QMutexLocker>

#include <U2Core/FeatureColors.h>

namespace U2 {

QMutex GBFeatureUtils::allKeys_mutex;
QMutex GBFeatureUtils::getKeyGroups_mutex;
QMutex GBFeatureUtils::getKey_mutex;

const QByteArray GBFeatureUtils::QUALIFIER_AMINO_STRAND("ugene_amino_strand");
const QByteArray GBFeatureUtils::QUALIFIER_AMINO_STRAND_YES("yes");
const QByteArray GBFeatureUtils::QUALIFIER_AMINO_STRAND_NO("no");

const QByteArray GBFeatureUtils::QUALIFIER_NAME("ugene_name");
const QByteArray GBFeatureUtils::QUALIFIER_GROUP("ugene_group");

const QString GBFeatureUtils::DEFAULT_KEY = GBFeatureUtils::getKeyInfo(GBFeatureKey_misc_feature).text;

const QString GBFeatureUtils::QUALIFIER_CUT = "cut";

static QColor cl(const QString& txt) {
    QColor res;
    if (txt != "000000") {
        res.setNamedColor("#"+txt);
    }
    return res;
}


#define FKE(key, text, color, amino, desc, quals) \
    features[key] = GBFeatureKeyInfo(key, text, color.isValid() ? color : FeatureColors::genLightColor(text), amino, desc); \
    if (strlen(quals) > 0) { \
        features[key].namingQuals = QString(quals).split(",", QString::SkipEmptyParts);\
    }

#define FK(key, text, color, amino, desc) \
    FKE(key, text, color, amino, desc, "label")


const QVector<GBFeatureKeyInfo>& GBFeatureUtils::allKeys() {
    QMutexLocker locker(&allKeys_mutex);
    static QVector<GBFeatureKeyInfo> features(GBFeatureKey_NUM_KEYS);
    static bool inited = false;
    if (inited) {
        return features;
    }
    inited = true;
    FK(GBFeatureKey_assembly_gap,     "assemlby_gap",     cl("000000"), false, QObject::tr("Gap between two components of a genome or transcriptome assembly"));
    FK(GBFeatureKey_attenuator,       "attenuator",       cl("000000"), false, QObject::tr("Sequence related to transcription termination"));
    FK(GBFeatureKey_bond,             "Bond",             cl("000000"), false, QObject::tr("Describes disulfide bonds (for protein files)"));
    FK(GBFeatureKey_C_region,         "C_region",         cl("000000"), false, QObject::tr("Span of the C immunological feature"));
    FK(GBFeatureKey_CAAT_signal,      "CAAT_signal",      cl("000000"), false, QObject::tr("`CAAT box' in eukaryotic promoters"));
    FKE(GBFeatureKey_CDS,             "CDS",             cl("9bffff"), true,  QObject::tr("Sequence coding for amino acids in protein (includes stop codon)"),
        "label,protein_id,locus_tag,gene,function,product");
    FK(GBFeatureKey_conflict,         "conflict",         cl("000000"), false, QObject::tr("Independent sequence determinations differ"));
    FK(GBFeatureKey_centromere,       "centromere",       cl("000000"), false, QObject::tr("Region of biological interest identified as a centromere and which has been experimentally characterized"));
    FK(GBFeatureKey_D_loop,           "D-loop",           cl("000000"), false, QObject::tr("Displacement loop"));
    FK(GBFeatureKey_D_segment,        "D_segment",        cl("000000"), false, QObject::tr("Span of the D immunological feature"));
    FK(GBFeatureKey_enhancer,         "enhancer",         cl("000000"), false, QObject::tr("Cis-acting enhancer of promoter function"));
    FK(GBFeatureKey_exon,             "exon",             cl("000000"), false, QObject::tr("Region that codes for part of spliced mRNA"));
    FK(GBFeatureKey_gap,              "gap",              cl("000000"), false, QObject::tr("Gap in the sequence"));
    FKE(GBFeatureKey_gene,            "gene",             cl("00ffc8"), false, QObject::tr("Region that defines a functional gene, possibly including upstream (promotor, enhancer, etc) and downstream control elements, and for which a name has been assigned."),
        "label,gene,locus_tag,product,function");
    FK(GBFeatureKey_GC_signal,        "GC_signal",        cl("000000"), false, QObject::tr("`GC box' in eukaryotic promoters"));
    FK(GBFeatureKey_iDNA,             "iDNA",             cl("000000"), false, QObject::tr("Intervening DNA eliminated by recombination"));
    FK(GBFeatureKey_intron,           "intron",           cl("000000"), false, QObject::tr("Transcribed region excised by mRNA splicing"));
    FK(GBFeatureKey_J_region,         "J_region",         cl("000000"), false, QObject::tr("Span of the J immunological feature"));
    FK(GBFeatureKey_LTR,              "LTR",              cl("000000"), false, QObject::tr("Long terminal repeat"));
    FK(GBFeatureKey_mat_peptide,      "mat_peptide",      cl("000000"), true,  QObject::tr("Mature peptide coding region (does not include stop codon)"));
    FK(GBFeatureKey_misc_binding,     "misc_binding",     cl("000000"), false, QObject::tr("Miscellaneous binding site"));
    FK(GBFeatureKey_misc_difference,  "misc_difference",  cl("000000"), false, QObject::tr("Miscellaneous difference feature"));
    FKE(GBFeatureKey_misc_feature,    "misc_feature",     cl("000000"), false, QObject::tr("Region of biological significance that cannot be described by any other feature")
        , "label,note");
    FK(GBFeatureKey_misc_recomb,      "misc_recomb",      cl("000000"), false, QObject::tr("Miscellaneous, recombination feature"));
    FK(GBFeatureKey_misc_RNA,         "misc_RNA",         cl("000000"), false, QObject::tr("Miscellaneous transcript feature not defined by other RNA keys"));
    FK(GBFeatureKey_misc_signal,      "misc_signal",      cl("000000"), false, QObject::tr("Miscellaneous signal"));
    FK(GBFeatureKey_misc_structure,   "misc_structure",   cl("000000"), false, QObject::tr("Miscellaneous DNA or RNA structure"));
    FK(GBFeatureKey_mobile_element,   "mobile_element",   cl("000000"), false, QObject::tr("Region of genome containing mobile elements"));
    FK(GBFeatureKey_modified_base,    "modified_base",    cl("000000"), false, QObject::tr("The indicated base is a modified nucleotide"));
    FK(GBFeatureKey_mRNA,             "mRNA",             cl("000000"), false, QObject::tr("Messenger RNA"));
    FK(GBFeatureKey_ncRNA,            "ncRNA",            cl("000000"), false, QObject::tr("A non-protein-coding gene, other than ribosomal RNA and transfer RNA, the functional molecule of which is the RNA transcrip"))
    FK(GBFeatureKey_N_region,         "N_region",         cl("000000"), false, QObject::tr("Span of the N immunological feature"));
    FK(GBFeatureKey_old_sequence,     "old_sequence",     cl("000000"), false, QObject::tr("Presented sequence revises a previous version"));
    FK(GBFeatureKey_operon,           "operon",           cl("000000"), false, QObject::tr("Region containing polycistronic transcript including a cluster of genes that are under the control of the same regulatory sequences/promotor and in the same biological pathway"));
    FK(GBFeatureKey_oriT,             "oriT",             cl("000000"), false, QObject::tr("Origin of transfer; region of a DNA molecule where transfer is initiated during the process of conjugation or mobilization"));
    FK(GBFeatureKey_polyA_signal,     "polyA_signal",     cl("000000"), false, QObject::tr("Signal for cleavage & polyadenylation"));
    FK(GBFeatureKey_polyA_site,       "polyA_site",       cl("000000"), false, QObject::tr("Site at which polyadenine is added to mRNA"));
    FK(GBFeatureKey_precursor_RNA,    "precursor_RNA",    cl("000000"), false, QObject::tr("Any RNA species that is not yet the mature RNA product"));
    FK(GBFeatureKey_prim_transcript,  "prim_transcript",  cl("000000"), false, QObject::tr("Primary (unprocessed) transcript"));
    FK(GBFeatureKey_primer,           "primer",           cl("000000"), false, QObject::tr("Primer binding region used with PCR"));
    FK(GBFeatureKey_primer_bind,      "primer_bind",      cl("000000"), false, QObject::tr("Non-covalent primer binding site"));
    FK(GBFeatureKey_promoter,         "promoter",         cl("000000"), false, QObject::tr("A region involved in transcription initiation"));
    FK(GBFeatureKey_protein_bind,     "protein_bind",     cl("000000"), false, QObject::tr("Non-covalent protein binding site on DNA or RNA"));
    FK(GBFeatureKey_RBS,              "RBS",              cl("000000"), false, QObject::tr("Ribosome binding site"));
    FK(GBFeatureKey_rep_origin,       "rep_origin",       cl("000000"), false, QObject::tr("Replication origin for duplex DNA"));
    FK(GBFeatureKey_repeat_region,    "repeat_region",    cl("ccccff"), false, QObject::tr("Sequence containing repeated subsequences"));
    FK(GBFeatureKey_repeat_unit,      "repeat_unit",      cl("ccccff"), false, QObject::tr("One repeated unit of a repeat_region"));
    FK(GBFeatureKey_rRNA,             "rRNA",             cl("000000"), false, QObject::tr("Ribosomal RNA"));
    FK(GBFeatureKey_S_region,         "S_region",         cl("000000"), false, QObject::tr("Span of the S immunological feature"));
    FK(GBFeatureKey_satellite,        "satellite",        cl("000000"), false, QObject::tr("Satellite repeated sequence"));
    FK(GBFeatureKey_scRNA,            "scRNA",            cl("000000"), false, QObject::tr("Small cytoplasmic RNA"));
    FK(GBFeatureKey_sig_peptide,      "sig_peptide",      cl("000000"), false, QObject::tr("Signal peptide coding region"));
    FK(GBFeatureKey_snRNA,            "snRNA",            cl("000000"), false, QObject::tr("Small nuclear RNA"));
    FK(GBFeatureKey_source,           "source",           cl("cccccc"), false, QObject::tr("Identifies the biological source of the specified span of the sequence"));
    FK(GBFeatureKey_stem_loop,        "stem_loop",        cl("000000"), false, QObject::tr("Hair-pin loop structure in DNA or RNA"));
    FK(GBFeatureKey_STS,              "STS",              cl("00dcdc"), false, QObject::tr("Sequence Tagged Site; operationally unique sequence that identifies the combination of primer spans used in a PCR assay"));
    FK(GBFeatureKey_TATA_signal,      "TATA_signal",      cl("000000"), false, QObject::tr("`TATA box' in eukaryotic promoters"));
    FK(GBFeatureKey_telomere,         "telomere",         cl("000000"), false, QObject::tr("Region of biological interest identified as a telomere and which has been experimentally characterized"));
    FK(GBFeatureKey_terminator,       "terminator",       cl("000000"), false, QObject::tr("Sequence causing transcription termination"));
    FK(GBFeatureKey_tmRNA,            "tmRNA",            cl("000000"), false, QObject::tr("Transfer messenger RNA; tmRNA acts as a tRNA first, and then as an mRNA that encodes a peptide tag; the ribosome translates this mRNA region of tmRNA and attaches the encoded peptide tag to the C-terminus of the unfinished protein; this attached tag targets the protein for destruction or proteolysis"));
    FK(GBFeatureKey_transit_peptide,  "transit_peptide",  cl("000000"), false, QObject::tr("Transit peptide coding region"));
    FK(GBFeatureKey_transposon,       "transposon",       cl("000000"), false, QObject::tr("Transposable element (TN)"));
    FK(GBFeatureKey_tRNA,             "tRNA",             cl("c8fac8"), false, QObject::tr("Transfer RNA"));
    FK(GBFeatureKey_unsure,           "unsure",           cl("000000"), false, QObject::tr("Authors are unsure about the sequence in this region"));
    FK(GBFeatureKey_V_region,         "V_region",         cl("000000"), false, QObject::tr("Span of the V immunological feature"));
    FK(GBFeatureKey_V_segment,        "V_segment",        cl("000000"), false, QObject::tr("Variable segment of immunoglobulin light and heavy chains, and T-cell receptor alpha, beta, and gamma chains; codes for most of the variable region (V_region) and the last few amino acids of the leader peptide"));
    FK(GBFeatureKey_variation,        "variation",        cl("ffff9b"), false, QObject::tr("A related population contains stable mutation"));
    FK(GBFeatureKey__10_signal,       "-10_signal",       cl("000000"), false, QObject::tr("`Pribnow box' in prokaryotic promoters"));
    FK(GBFeatureKey__35_signal,       "-35_signal",       cl("000000"), false, QObject::tr("`-35 box' in prokaryotic promoters"));
    FK(GBFeatureKey_3_clip,           "3'clip",           cl("000000"), false, QObject::tr("3'-most region of a precursor transcript removed in processing"));
    FK(GBFeatureKey_3_UTR,            "3'UTR",            cl("ffcde6"), false, QObject::tr("3' untranslated region (trailer)"));
    FK(GBFeatureKey_5_clip,           "5'clip",           cl("000000"), false, QObject::tr("5'-most region of a precursor transcript removed in processing"));
    FK(GBFeatureKey_5_UTR,            "5'UTR",            cl("ffc8c8"), false, QObject::tr("5' untranslated region (leader)"));
    FK(GBFeatureKey_Protein,          "Protein",          cl("000000"), false, QObject::tr("'Protein' feature key"));
    FK(GBFeatureKey_Region,           "Region",           cl("000000"), false, QObject::tr("'Region' feature key"));
    FK(GBFeatureKey_Site,             "Site",             cl("000000"), false, QObject::tr("'Site' feature key"));

#ifdef _DEBUG
    for (int i=0; i<features.size(); i++) {
        assert(features[i].id != GBFeatureKey_UNKNOWN && features[i].id == i);
        assert(!features[i].text.isEmpty());
        assert(!features[i].desc.isEmpty());
    }
#endif

    return features;
}


const QMultiMap<QString, GBFeatureKey>& GBFeatureUtils::getKeyGroups() {
    QMutexLocker locker(&getKeyGroups_mutex);
    static QMultiMap<QString, GBFeatureKey> groups;

    if (groups.isEmpty()) {
        QString genes = QObject::tr("Genes");
        groups.insert(genes, GBFeatureKey_CDS);
        groups.insert(genes, GBFeatureKey_exon);
        groups.insert(genes, GBFeatureKey_gene);
        groups.insert(genes, GBFeatureKey_intron);
        groups.insert(genes, GBFeatureKey_mRNA);
        groups.insert(genes, GBFeatureKey_polyA_site);
        groups.insert(genes, GBFeatureKey_precursor_RNA);
        groups.insert(genes, GBFeatureKey_prim_transcript);
        groups.insert(genes, GBFeatureKey_promoter);
        groups.insert(genes, GBFeatureKey_3_clip);
        groups.insert(genes, GBFeatureKey_3_UTR);
        groups.insert(genes, GBFeatureKey_5_clip);
        groups.insert(genes, GBFeatureKey_5_UTR);

        QString signls = QObject::tr("Signals");
        groups.insert(signls, GBFeatureKey_attenuator);
        groups.insert(signls, GBFeatureKey_CAAT_signal);
        groups.insert(signls, GBFeatureKey_GC_signal);
        groups.insert(signls, GBFeatureKey_enhancer);
        groups.insert(signls, GBFeatureKey_mat_peptide);
        groups.insert(signls, GBFeatureKey_misc_signal);
        groups.insert(signls, GBFeatureKey_oriT);
        groups.insert(signls, GBFeatureKey_rep_origin);
        groups.insert(signls, GBFeatureKey_polyA_site);
        groups.insert(signls, GBFeatureKey_polyA_signal);
        groups.insert(signls, GBFeatureKey_promoter);
        groups.insert(signls, GBFeatureKey_RBS);
        groups.insert(signls, GBFeatureKey_sig_peptide);
        groups.insert(signls, GBFeatureKey_terminator);
        groups.insert(signls, GBFeatureKey_transit_peptide);
        groups.insert(signls, GBFeatureKey_TATA_signal);
        groups.insert(signls, GBFeatureKey__35_signal);
        groups.insert(signls, GBFeatureKey__10_signal);

        QString binding = QObject::tr("Binding");
        groups.insert(binding, GBFeatureKey_misc_binding);
        groups.insert(binding, GBFeatureKey_primer_bind);
        groups.insert(binding, GBFeatureKey_protein_bind);
        groups.insert(binding, GBFeatureKey_bond);

        QString variation = QObject::tr("Variation");
        groups.insert(variation, GBFeatureKey_conflict);
        groups.insert(variation, GBFeatureKey_misc_difference);
        groups.insert(variation, GBFeatureKey_modified_base);
        groups.insert(variation, GBFeatureKey_old_sequence);
        groups.insert(variation, GBFeatureKey_STS);
        groups.insert(variation, GBFeatureKey_unsure);
        groups.insert(variation, GBFeatureKey_variation);

        QString repeats = QObject::tr("Repeats");
        groups.insert(repeats, GBFeatureKey_LTR);
        groups.insert(repeats, GBFeatureKey_repeat_region);
        groups.insert(repeats, GBFeatureKey_repeat_unit);
        groups.insert(repeats, GBFeatureKey_satellite);
        groups.insert(repeats, GBFeatureKey_transposon); //TODO: recheck grouping

        QString rna = QObject::tr("RNA");
        groups.insert(rna, GBFeatureKey_ncRNA);
        groups.insert(rna, GBFeatureKey_misc_RNA);
        groups.insert(rna, GBFeatureKey_mRNA);
        groups.insert(rna, GBFeatureKey_rRNA);
        groups.insert(rna, GBFeatureKey_scRNA);
        groups.insert(rna, GBFeatureKey_snRNA);
        groups.insert(rna, GBFeatureKey_tRNA);
        groups.insert(rna, GBFeatureKey_tmRNA);

        QString misc = QObject::tr("Misc");
        groups.insert(misc, GBFeatureKey_assembly_gap);
        groups.insert(misc, GBFeatureKey_centromere);
        groups.insert(misc, GBFeatureKey_D_loop);
        groups.insert(misc, GBFeatureKey_gap);
        groups.insert(misc, GBFeatureKey_iDNA);
        groups.insert(misc, GBFeatureKey_misc_binding);
        groups.insert(misc, GBFeatureKey_misc_feature);
        groups.insert(misc, GBFeatureKey_misc_recomb);
        groups.insert(misc, GBFeatureKey_misc_structure);
        groups.insert(misc, GBFeatureKey_mobile_element);
        groups.insert(misc, GBFeatureKey_operon);
        groups.insert(misc, GBFeatureKey_primer); //TODO: recheck grouping
        groups.insert(misc, GBFeatureKey_source);
        groups.insert(misc, GBFeatureKey_stem_loop);
        groups.insert(misc, GBFeatureKey_telomere);
        groups.insert(misc, GBFeatureKey_Protein);
        groups.insert(misc, GBFeatureKey_Region);
        groups.insert(misc, GBFeatureKey_Site);

        QString spans = QObject::tr("Spans");
        groups.insert(spans, GBFeatureKey_C_region);
        groups.insert(spans, GBFeatureKey_D_segment);
        groups.insert(spans, GBFeatureKey_J_region);
        groups.insert(spans, GBFeatureKey_N_region);
        groups.insert(spans, GBFeatureKey_S_region);
        groups.insert(spans, GBFeatureKey_V_region);
        groups.insert(spans, GBFeatureKey_V_segment);

#ifdef _DEBUG
        //check that no feature lost
        QVector<bool> featureInGroup(GBFeatureKey_NUM_KEYS, false);
        foreach(const QString& groupName, groups.keys()) {
            QList<GBFeatureKey> values = groups.values(groupName);
            foreach (GBFeatureKey k, values) {
                featureInGroup[k] = true;
            }
        }
        for (int i=0; i<GBFeatureKey_NUM_KEYS; i++) {
            GBFeatureKey fk = GBFeatureKey(i);
            assert(featureInGroup[fk]);
        }
#endif
    }
    return groups;
}

bool GBFeatureUtils::isFeatureHasNoValue(const QString& featureName)
{
    if (featureName == "pseudo")
    {
        return true;
    }
    return false;
}

GBFeatureKey GBFeatureUtils::getKey(const QString& text) {
    QMutexLocker locker(&getKey_mutex);
    static QHash<QString, GBFeatureKey> keysByText;
    if (keysByText.isEmpty()) {
        foreach(const GBFeatureKeyInfo& ki, allKeys()) {
            keysByText[ki.text] = ki.id;
        }
    }
    return keysByText.value(text, GBFeatureKey_UNKNOWN);
}

}//namespace
