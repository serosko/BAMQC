#include <iostream>
#include <seqan/bam_io.h>
#include <parse.h>

using namespace seqan;

//Get Distance from leftmost base of first mate to rightmost base of right mate 
//(template length) from one record and add it to counter in TInsertDistr
inline int countInsertSize(TInsertDistr & counts, const BamAlignmentRecord & record, ProgramOptions & options)
{
    int32_t insertSize = record.tLen;
    if (insertSize >= 0 && insertSize <= options.maxInsert)
    {
        ++counts[insertSize];                     //Increase counter
        return 0;
    }
    else return 1;          //return 1 if record was right mate or longer than maxInsert (not counted)
}

//checks the flags and mapping quality of a record and returns false if any undesired flag is found, true otherwise.
inline bool checkRecord (const BamAlignmentRecord & record, const ProgramOptions & options)
{
    if(hasFlagDuplicate(record) ||
        hasFlagQCNoPass(record) ||
        hasFlagSecondary(record)||
        hasFlagSupplementary(record) ||
        record.mapQ < options.minMapQ)
        return false;
    else return true;
}

//wrapper for counting the insert sizes of the whole Bam-file
inline int wrapCountInsertSize(TInsertDistr & counts, BamFileIn & bamFile, ProgramOptions & options)
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
        std::cout << "Error: "  << e.what() << std::endl;
        return 1;
    }
    return 0;
}