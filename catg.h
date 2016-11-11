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
//Return number of hits in next sequence with hits
inline int findNextTriplet(String<unsigned> & occ,
                    BamFileIn & bamFile,
                    BamAlignmentRecord & record,
                    const ProgramOptions & options)
{
    while (!atEnd(bamFile))
    {
        readRecord(record, bamFile);
        if (checkRecordCAGT(record, options))
        {
            Dna5String needle;
            if (hasFlagFirst(record))                           //Select proper pattern for first/last mate
                needle = (Dna5String)"CTG";
            else if (hasFlagLast(record))
                needle = (Dna5String)"CAG";
            else continue;                                      //Skip record without proper first/second mate flag.
            unsigned c = findTriplet(occ, (Dna5String)record.seq, needle);
            if (c != 0)                                         //If there were occurences of the pattern...
                return c;
        }
    }
    return 0;
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
inline bool checkContext(const Dna5String & ref, bool firstMate)
{
    if (firstMate)
        return (ref == "CGG");
    else
        return (ref == "CCG");
}

//Wrapper for calling all functions necessary for finding all CCG > CAG or CGG > CTG occurences.
//Returns the number of occurences.
inline unsigned getArtifactCount(BamFileIn & bamFile, FaiIndex & faiIndex, const ProgramOptions & options)
{
    unsigned hits = 0;                  //counter for verified hits.
    BamAlignmentRecord record;
    String<unsigned> occ = "";
    reserve(occ, 5);
    Dna5String ref = "";                //Will hold triplet of reference after call of getRefAt
    while (!atEnd(bamFile))
    {
        unsigned c = 0;                 //counter for number of occurences in next record having any occurences
        c = findNextTriplet(occ, bamFile, record, options);
        if (c > 0)
            for (unsigned i = 0; i < c; ++i)
            {
                getRefAt(ref, faiIndex,  getContigName(record, bamFile), occ[i]);
                std::cout << hasFlagFirst(record) << ":";               //Todo Look why only CGG is reported as ref, even in reverse strand
                std::cout << ref << std::endl;
                if (checkContext(ref, hasFlagFirst(record)))
                    ++hits;
            }
    }
    return hits;
}
