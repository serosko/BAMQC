#include<BAMQC.h>
#include <seqan/find.h>

using namespace seqan;

//Maybe filter for primary alignments and quality, duplicates etc.
//Check if read contains an CAG (or CTG respectively if second mate)
//Check CIGAR, if middle pos is a (mis)match (maybe drop this step)
//Compare, if Ref has an C (or G repsectively if second mate) at that position.

//checks the flags and mapping quality of a record and returns false if any undesired flag is found, true otherwise.
inline bool checkRecordCAGT (const BamAlignmentRecord & record, const ProgramOptions & options)
{
    if(hasFlagDuplicate(record) ||
        hasFlagQCNoPass(record) ||
        hasFlagSecondary(record)||
        hasFlagSupplementary(record) ||
        hasFlagUnmapped(record) ||
        record.mapQ < options.minMapQCAGT)
        return false;
    else return true;
}

//scans a String for occurences of "needle" using shift-or algorithm, saves occurences in occ and return their ammount
inline int findTriplet(String<unsigned> & occ, Dna5String haystack, Dna5String & needle)
{
    Finder<Dna5String> finder(haystack);
    Pattern<Dna5String, ShiftOr> pattern(needle);
    unsigned c = 0;
    while (find(finder, pattern))
    {
            appendValue(occ, beginPosition(finder));
            ++c;
    }
    return c;
}

//Check each sequence using findTriplet() with the repective pattern until the pattern occurs.
//Call checkRecordCAGT before looking for triplet.
inline int findNextTriplet(String<unsigned> & occ,
                           String<unsigned> & nocc,
                           BamFileIn & bamFile,
                           BamAlignmentRecord & record,
                           const ProgramOptions & options)
{
    unsigned c = 0;
    unsigned nc = 0;
    while (!atEnd(bamFile))
    {
        readRecord(record, bamFile);
        if (checkRecordCAGT(record, options))
        {
            Dna5String needle;
            Dna5String revNeedle;
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
            }
            else
            {
                continue;                                       //Skip record without proper first/second mate flag.
            }
            c = findTriplet(occ, (Dna5String)record.seq, needle);
            for (unsigned i = 0; i < length(occ); ++i)          //Correct pos in occ for starting position of alignment
                occ[i] += record.beginPos;
            nc = findTriplet(nocc, (Dna5String)record.seq, revNeedle);
            for (unsigned i = 0; i < length(nocc); ++i)             //Correct pos in occ for starting position of alignment
                nocc[i] += record.beginPos;
            if (c || nc)                                            //If there were occurences of the pattern...
                return c;
        }
    }
    return c;
}

//Takes a sequence id (chr) and a position and returns the triplet of the reference genome at that position +- 1
inline int getRefAt (Dna5String & ref, FaiIndex & faiIndex, CharString id, unsigned pos)
{
    unsigned idx = 0;                                           //Holder for position of id in index
    ref = "";
    if (!getIdByName(idx, faiIndex, id))                        //Get position of sequence by its ID and save it in idx
    {
        std::cout << "WARNING: Cannot find ID " << id << " in index. Skipping...\n";
        return 1;
    }
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

//Wrapper for calling all functions necessary for finding all CCG > CAG or CGG > CTG occurences and non-artifacts.
//Returns the number of occurences as i1 and number of non-artifactual conversions as i2
inline Pair<unsigned, unsigned> getArtifactCount(BamFileIn & bamFile, FaiIndex & faiIndex, const ProgramOptions & options)
{
    unsigned hits = 0;                  //counter for verified hits.
    unsigned nonHits = 0;               //counter for non artifactual conversions
    BamAlignmentRecord record;
    String<unsigned> occ = "";
    reserve(occ, 5);
    String<unsigned> nocc = "";
    reserve(occ, 5);
    Dna5String ref = "";                //Will hold triplet of reference after call of getRefAt
    while (!atEnd(bamFile))
    {
        clear(occ);
        clear(nocc);
        findNextTriplet(occ, nocc, bamFile, record, options);
        for (unsigned i = 0; i < length(occ); ++i)
        {
            getRefAt(ref, faiIndex,  getContigName(record, bamFile), occ[i]);
            if (checkContext(ref, hasFlagFirst(record), hasFlagRC(record)))
                ++hits;
        }
        for (unsigned j = 0; j < length(nocc); ++j)
        {
            getRefAt(ref, faiIndex,  getContigName(record, bamFile), nocc[j]);
            if (checkNAContext(ref, hasFlagFirst(record), hasFlagRC(record)))
                ++nonHits;
        }
    }
    return Pair<unsigned, unsigned> (hits, nonHits);
}