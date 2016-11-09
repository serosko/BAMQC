#include <iostream>
#include <seqan/bam_io.h>

//TODO: Do this for coordinate sorted BAM files.

using namespace seqan;

//String holding the number of inserts of each length. Index 1 holds the number
//of inserts with length 1, index 2 holds the number of inserts with length 2...
//index 0 holds the number of segments without a mapped partner or the i
//information is not available
typedef String<unsigned> TInsertDistr;

//Get Distance from leftmost base of first mate to rightmost base of right mate 
//(template length) from one record and add it to counter in TInsertDistr
int countInsertSize(TInsertDistr & counts, const BamAlignmentRecord & record)
{
    int32_t insertSize = record.tLen;
    if (insertSize >= 0)
    {
        if ((unsigned)insertSize > length(counts)) //Check for sufficient space
            resize(counts, insertSize * 2, 0);
        ++counts[insertSize];                     //Increase counter
        return 0;
    }
    else return 1;          //return 1 if record was right mate (not counted)
}
//wrapper for counting the insert sizes of the whole Bam-file
int wrapCountInsertSize(TInsertDistr & counts, BamFileIn & bamFile)
{
    BamAlignmentRecord record;
    try
    {
        while (!atEnd(bamFile))
        {
            readRecord(record, bamFile);
            countInsertSize(counts, record);
        }
    }
    catch (Exception const & e)
    {
        std::cout << "Error: "  << e.what() << std::endl;
        return 1;
    }
    return 0;
}
