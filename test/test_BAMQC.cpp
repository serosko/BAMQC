#undef SEQAN_ENABLE_TESTING
#define SEQAN_ENABLE_TESTING 1

#include "../BAMQC.h"
#include <seqan/basic.h>

using namespace seqan;

// -----------------------------------------------------------------------------
// Unit Tests for BAMQC
// -----------------------------------------------------------------------------

SEQAN_DEFINE_TEST(test_checkRecord)
{
    ProgramOptions options;
    options.minMapQ = 30;
    BamAlignmentRecord unmapped;
    BamAlignmentRecord firstMate;
    BamAlignmentRecord lowQual;
    unmapped.flag = 4;              //read unmapped
    unmapped.mapQ = 30;
    firstMate.flag = 64;            //first in pair
    firstMate.mapQ = 30;
    lowQual.flag = 64;              //first in pair
    lowQual.mapQ = 29;              //lowQuality
    SEQAN_ASSERT_EQ(checkRecord(unmapped, options), false);
    SEQAN_ASSERT_EQ(checkRecord(firstMate, options), true);
    SEQAN_ASSERT_EQ(checkRecord(lowQual, options), false);
}
SEQAN_DEFINE_TEST(test_countInsertSize)
{
    TInsertDistr counts;
    resize(counts, 100, 0);
    BamAlignmentRecord record;
    ProgramOptions options;
    options.maxInsert = 50;
    record.tLen = 51;
    SEQAN_ASSERT_EQ(countInsertSize(counts, record, options), 1);
    record.tLen = 0;
    SEQAN_ASSERT_EQ(countInsertSize(counts, record, options), 0);
    record.tLen = 50;
    SEQAN_ASSERT_EQ(countInsertSize(counts, record, options), 0);
    for (unsigned i = 0 ; i < length(counts); ++i)
    {
        if(i!=50 && i != 0)
            SEQAN_ASSERT_EQ(counts[i], 0u);
        else
            SEQAN_ASSERT_EQ(counts[i], 1u);
    }
}
SEQAN_DEFINE_TEST(test_findTriplet)
{
    Dna5String haystack = "CNNNNCTCGTAAACTGAAAATCGTCGAAAACT"; //AAAA at 16 and 26
    Dna5String needle1 = "CCCC";            //No hits
    Dna5String needle2 = "AAAA";            //Hits at 16 and 26
    String<unsigned>  occ = "";
    SEQAN_ASSERT_EQ(findTriplet(occ, haystack, needle1), false);
    SEQAN_ASSERT_EQ(length(occ), 0u);
    SEQAN_ASSERT_EQ(findTriplet(occ, haystack, needle2), true);
    SEQAN_ASSERT_EQ(length(occ), 2u);
    SEQAN_ASSERT_EQ(occ[0], 16u);
    SEQAN_ASSERT_EQ(occ[1], 26u);
}
SEQAN_DEFINE_TEST(test_getNeedles)
{
    Dna5String dummy1 = "";
    Dna5String dummy2 = "";
    Dna5String needleFirst = "";
    Dna5String revNeedleFirst = "";
    Dna5String needleSecond = "";
    Dna5String revNeedleSecond = "";
    Dna5String needleFirstRC = "";
    Dna5String revNeedleFirstRC = "";
    Dna5String needleSecondRC = "";
    Dna5String revNeedleSecondRC = "";
    BamAlignmentRecord  nothing;
    BamAlignmentRecord  first;
    BamAlignmentRecord  second;
    BamAlignmentRecord  firstRC;
    BamAlignmentRecord  secondRC;
    nothing.flag = 0;
    first.flag = 65;
    second.flag = 129;
    firstRC.flag = 81;
    secondRC.flag = 145;
    SEQAN_ASSERT_EQ(getNeedles(dummy1, dummy2, nothing), 0);
    SEQAN_ASSERT_EQ(getNeedles(needleFirst, revNeedleFirst, first), 1);
    SEQAN_ASSERT_EQ(getNeedles(needleSecond, revNeedleSecond, second), 1);
    SEQAN_ASSERT_EQ(getNeedles(needleFirstRC, revNeedleFirstRC, firstRC), 1);
    SEQAN_ASSERT_EQ(getNeedles(needleSecondRC, revNeedleSecondRC, secondRC), 1);
    SEQAN_ASSERT_EQ(needleFirst[1], (Dna5)'T');
    SEQAN_ASSERT_EQ(revNeedleFirst[1], (Dna5)'A');
    SEQAN_ASSERT_EQ(needleSecond[1], (Dna5)'A');
    SEQAN_ASSERT_EQ(revNeedleSecond[1], (Dna5)'T');
    SEQAN_ASSERT_EQ(needleFirstRC[1], (Dna5)'A');
    SEQAN_ASSERT_EQ(revNeedleFirstRC[1], (Dna5)'T');
    SEQAN_ASSERT_EQ(needleSecondRC[1], (Dna5)'T');
    SEQAN_ASSERT_EQ(revNeedleSecondRC[1], (Dna5)'A');
}
SEQAN_DEFINE_TEST(test_checkContext)
{
    Dna5String cgg = "CGG";
    Dna5String ccg = "CCG";
    Dna5String nnn = "NNN";
    SEQAN_ASSERT_EQ(checkContext(cgg, true, false), true);
    SEQAN_ASSERT_EQ(checkContext(cgg, false, false), false);
    SEQAN_ASSERT_EQ(checkContext(cgg, true, true), false);
    SEQAN_ASSERT_EQ(checkContext(cgg, false, true), true);
    SEQAN_ASSERT_EQ(checkContext(ccg, true, false), false);
    SEQAN_ASSERT_EQ(checkContext(ccg, false, false), true);
    SEQAN_ASSERT_EQ(checkContext(ccg, true, true), true);
    SEQAN_ASSERT_EQ(checkContext(ccg, false, true), false);
    SEQAN_ASSERT_EQ(checkContext(nnn, true, false), false);
    SEQAN_ASSERT_EQ(checkContext(nnn, false, false), false);
    SEQAN_ASSERT_EQ(checkContext(nnn, true, true), false);
    SEQAN_ASSERT_EQ(checkContext(nnn, false, true), false);
}

SEQAN_BEGIN_TESTSUITE(test_BAMQC)
{
    SEQAN_CALL_TEST(test_checkRecord);
    SEQAN_CALL_TEST(test_countInsertSize);
    SEQAN_CALL_TEST(test_findTriplet);
    SEQAN_CALL_TEST(test_getNeedles);
    SEQAN_CALL_TEST(test_checkContext);
}
SEQAN_END_TESTSUITE