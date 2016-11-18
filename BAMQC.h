//Author: Sebastian Roskosch <Sebastian.Roskosch[at]bihealth.de>
#include <iostream>
#include <seqan/bam_io.h>
#include <parse.h>
#include <seqan/find.h>

using namespace seqan;

////////////////////General functions//////////////////////////
//checks the flags and mapping quality of a record and returns false if any undesired flag is found, true otherwise.
inline bool checkRecord (const BamAlignmentRecord & record, const ProgramOptions & options)
{
    if(hasFlagDuplicate(record) ||
        hasFlagQCNoPass(record) ||
        hasFlagSecondary(record)||
        hasFlagSupplementary(record) ||
        hasFlagUnmapped(record) ||
        record.mapQ < options.minMapQ)
        return false;
    else return true;
}
////////////////////Insert-size distribution functions//////////////////////////
//Get Distance from leftmost base of first mate to rightmost base of right mate 
//(template length) from one record and add it to counter in TInsertDistr
inline int countInsertSize(TInsertDistr & counts, const BamAlignmentRecord & record, const ProgramOptions & options)
{
    int32_t insertSize = record.tLen;
    if (insertSize >= 0 && insertSize <= options.maxInsert)
    {
        ++counts[insertSize];                     //Increase counter
        return 0;
    }
    else return 1;          //return 1 if record was right mate or longer than maxInsert (not counted)
}
//wrapper for counting the insert sizes of the whole Bam-file
inline bool wrapCountInsertSize(TInsertDistr & counts, BamFileIn & bamFile, ProgramOptions & options)
{
    resize(counts, options.maxInsert + 1, 0);
    BamAlignmentRecord record;
    try
    {
        while (!atEnd(bamFile))
        {
            readRecord(record, bamFile);
            if (checkRecord(record, options))
                countInsertSize(counts, record, options);
        }
    }
    catch (Exception const & e)
    {
        std::cerr << "Error: "  << e.what() << std::endl;
        return false;
    }
    return true;
}
//////////////////Artifact-check functions/////////////////////
//scans a String for occurences of "needle" using shift-or algorithm, saves occurences in occ return true if found
inline bool findTriplet(String<unsigned> & occ, Dna5String haystack, Dna5String & needle)
{
    Finder<Dna5String> finder(haystack);
    Pattern<Dna5String, ShiftOr> pattern(needle);
    while (find(finder, pattern))
        appendValue(occ, beginPosition(finder));
    return (length(occ) > 0);
}
//Get forward and reverse needle.
//return 1 if everything is ok, and 0 if the flags are not properly set
inline int getNeedles(Dna5String & needle, Dna5String & revNeedle, BamAlignmentRecord & record)
{
    if (hasFlagFirst(record))                           //Select proper pattern for first/last mate of read pair
    {
        if (!hasFlagRC(record))                         //Also consider if read is reverse complement
        {
            needle = (Dna5String)"CTG";
            revNeedle = (Dna5String)"CAG";
        }
        else
        {
            needle = (Dna5String)"CAG";
            revNeedle = (Dna5String)"CTG";
        }
        return 1;
    }
    else if (hasFlagLast(record))
    {
        if (!hasFlagRC(record))
        {
            needle = (Dna5String)"CAG";
            revNeedle = (Dna5String)"CTG";
        }
        else
        {
            needle = (Dna5String)"CTG";
            revNeedle = (Dna5String)"CAG";
        }
        return 1;
    }
    else
        return 0;
}
//Return true, if the record contains a relevant triplet, false otherwise
inline bool findNextTriplet(String<unsigned> & occ,
                           String<unsigned> & nocc,
                           BamAlignmentRecord & record)
{
    unsigned c = 0;
    unsigned nc = 0;
    Dna5String needle;
    Dna5String revNeedle;
    if (!getNeedles(needle, revNeedle, record))         //Select proper pattern for first/last mate of read pair
        return false;                                       //Skip record without proper first/second mate flag.
    clear(occ);
    clear(nocc);
    c = findTriplet(occ, (Dna5String)record.seq, needle);
    for (unsigned i = 0; i < length(occ); ++i)          //Correct pos in occ for starting position of alignment
        occ[i] += record.beginPos;
    nc = findTriplet(nocc, (Dna5String)record.seq, revNeedle);
    for (unsigned i = 0; i < length(nocc); ++i)             //Correct pos in occ for starting position of alignment
        nocc[i] += record.beginPos;
    return (c || nc);                                           //If there were occurences of the pattern...
}
//Check each sequence using findTriplet() with the repective pattern until the pattern occurs.
//Additionally performs insert-size counting
//Return true if the current record has a relevant triplet, false otherwise
inline bool findNextTripletAndCountInsert(String<unsigned> & occ,
                           String<unsigned> & nocc,
                           TInsertDistr & InsertCounts,
                           BamAlignmentRecord & record,
                           const ProgramOptions & options)
{
    Dna5String needle;
    Dna5String revNeedle;
    if (!getNeedles(needle, revNeedle, record))         //Select proper pattern for first/last mate of read pair
        return false;                                   //Skip record without proper first/second mate flag.
    unsigned c = 0;
    unsigned nc = 0;
    clear(occ);
    clear(nocc);
    c = findTriplet(occ, (Dna5String)record.seq, needle);
    for (unsigned i = 0; i < length(occ); ++i)          //Correct pos in occ for starting position of alignment
        occ[i] += record.beginPos;
    nc = findTriplet(nocc, (Dna5String)record.seq, revNeedle);
    for (unsigned i = 0; i < length(nocc); ++i)             //Correct pos in occ for starting position of alignment
        nocc[i] += record.beginPos;
    countInsertSize(InsertCounts, record, options);
    return (c || nc);                                            //If there were occurences of the pattern...
}
//Takes a sequence id (chr) and a position and returns the triplet of the reference genome at that position +- 1
inline int getRefAt (Dna5String & ref, FaiIndex & faiIndex, CharString id, unsigned pos)
{
    unsigned idx = 0;                                           //Holder for position of id in index
    ref = "";
    getIdByName(idx, faiIndex, id);                              //Get position of sequence by its ID and save it in idx
    if (pos + 1 > sequenceLength(faiIndex, idx))       //Make sure the position lies within the boundaries of the index
        pos = sequenceLength(faiIndex, idx) - 2;
    readRegion(ref, faiIndex, idx, pos, pos + 3);  //Get infix
    return 0;
}
//Takes the infix from the reference and and compares it to the artifact context.
//Return true if it is a CCG > CAG (on first mate) or CGG > CTG (second mate) change, false otherwise
inline bool checkContext(const Dna5String & ref, bool firstMate, bool rc)
{
    if (firstMate)
    {
        if (!rc)
            return (ref == "CGG");
        else
            return (ref == "CCG");
    }
    else
    {
        if (!rc)
            return (ref == "CCG");
        else
            return (ref == "CGG");
    } 
}
//Vice-versa than check Context.
inline bool checkNAContext(const Dna5String & ref, bool firstMate, bool rc)
{
    if (firstMate)
    {
        if (!rc)
            return (ref == "CCG");
        else
            return (ref == "CGG");
    }
    else
    {
        if (!rc)
            return (ref == "CGG");
        else
            return (ref == "CCG");
    } 
}
//Check the record and return false if it should be skippe, true otherwise.
inline bool checkAndSkip(CharString & contig,
                        CharString & previousContig,
                        BamAlignmentRecord & record,
                        BamFileIn & bamFile,
                        FaiIndex & faiIndex,
                        const ProgramOptions & options)
{
    if (!checkRecord(record, options))
        return false;
    unsigned idx = 0;
    contig = getContigName(record, bamFile);
    if(!getIdByName(idx, faiIndex, contig))
    {
        if(contig != previousContig)
        {
            std::cout << "WARNING: Cannot find contig " << contig << " in index. Skipping..." << std::endl;
            previousContig = contig;
        }
        return false;
    }
    return true;
}
//Wrapper for calling all functions necessary for finding all CCG > CAG or CGG > CTG occurences and non-artifacts.
//Returns the number of occurences as i1 and number of non-artifactual conversions as i2
inline bool getArtifactCount(unsigned (& artifactConv) [2][2],
                                                 unsigned (& normalConv) [2][2], 
                                                 BamFileIn & bamFile,
                                                 FaiIndex & faiIndex,
                                                 const ProgramOptions & options)
{
    unsigned hits = 0;                  //counter for verified hits.
    unsigned nonHits = 0;               //counter for non artifactual conversions
    BamAlignmentRecord record;
    String<unsigned> occ = "";
    reserve(occ, 5);
    String<unsigned> nocc = "";
    reserve(occ, 5);
    Dna5String ref = "";               //Will hold triplet of reference after call of getRefAt
    CharString previousContig = "";
    CharString contig = "";
    try
    {
        while (!atEnd(bamFile))
        {
            readRecord(record, bamFile);
            if(!checkAndSkip(contig, previousContig, record, bamFile, faiIndex, options))
                continue;
            if(!findNextTriplet(occ, nocc, record))
                continue;
            CharString id = getContigName(record, bamFile);
            bool isFirst = hasFlagFirst(record);
            bool isRC = hasFlagRC(record);
            for (unsigned i = 0; i < length(occ); ++i)
            {
                getRefAt(ref, faiIndex,  id, occ[i]);
                if (checkContext(ref, isFirst, isRC))
                {
                    ++hits;
                    ++artifactConv[isFirst][isRC];
                }
            }
            for (unsigned j = 0; j < length(nocc); ++j)
            {
                getRefAt(ref, faiIndex,  id, nocc[j]);
                if (checkNAContext(ref, isFirst, isRC))
                {
                    ++nonHits;
                    ++normalConv[isFirst][isRC];
                }
            }
        }
        return true;
    }
    catch (Exception const & e)
    {
        std::cerr << "Error: "  << e.what() << std::endl;
        return false;
    }
}
///////////////Wrapper function for calling both checks in one run/////////////
//Return 1 on error, 0 otherwise
inline bool wrapDoAll(unsigned (& artifactConv) [2][2],
                     unsigned (& normalConv) [2][2],
                     TInsertDistr & InsertCounts, 
                     BamFileIn & bamFile,
                     ProgramOptions & options)
{
    FaiIndex faiIndex;
    if(!loadRefIdx(faiIndex, toCString(options.refPath)))
        return 0;
    BamAlignmentRecord record;
    unsigned hits = 0;                  //counter for verified hits.
    unsigned nonHits = 0;               //counter for non artifactual conversions
    String<unsigned> occ = "";          //positions of artifacual triplets
    reserve(occ, 5);
    String<unsigned> nocc = "";         //positions of non-artifactual triplets
    reserve(occ, 5);
    Dna5String ref = "";               //Will hold triplet of reference after call of getRefAt
    CharString previousContig = "";
    CharString contig = "";
    try
    {
        while (!atEnd(bamFile))
        {
            readRecord(record, bamFile);
            if(!checkAndSkip(contig, previousContig, record, bamFile, faiIndex, options))
                continue;
            if(!findNextTripletAndCountInsert(occ, nocc, InsertCounts, record, options))
                continue;
            contig = getContigName(record, bamFile);
            bool isFirst = hasFlagFirst(record);
            bool isRC = hasFlagRC(record);
            for (unsigned i = 0; i < length(occ); ++i)
            {
                getRefAt(ref, faiIndex, contig, occ[i]);
                if (checkContext(ref, isFirst, isRC))
                {
                    ++hits;
                    ++artifactConv[isFirst][isRC];
                }
            }
            for (unsigned j = 0; j < length(nocc); ++j)
            {
                getRefAt(ref, faiIndex, contig, nocc[j]);
                if (checkNAContext(ref, isFirst, isRC))
                {
                    ++nonHits;
                    ++normalConv[isFirst][isRC];
                }
            }
        }
        return true;
    }
    catch (Exception const & e)
    {
        std::cerr << "Error: "  << e.what() << std::endl;
        return false;
    }
}