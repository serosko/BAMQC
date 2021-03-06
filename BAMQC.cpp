//Author: Sebastian Roskosch <Sebastian.Roskosch[at]bihealth.de>
#include "BAMQC.h"

int main(int argc, char const ** argv)
{
    ProgramOptions options;
    //Initialize parser
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options,
                                                              argc,
                                                              argv);
    if (res != seqan::ArgumentParser::PARSE_OK) 
        return res == seqan::ArgumentParser::PARSE_ERROR;               //Terminate on parsing errors
    if (inputCheck(options))                                            //Terminate if check of parameters fails.
        return 1;
    feedBack(options);
    BamFileIn bamFile;                                                  //Prepare and load BAM-file
    if (!loadBAM(bamFile, options.inPath))
        return 1;
    BamHeader header;                                                   //Read header to get to right position in file
    readHeader(header, bamFile);
    if (options.conv && options.insDist)                                //Perform both checks in one run  to reduce I/O
    {
        TInsertDistr InsertCounts = "";
        resize(InsertCounts, options.maxInsert + 1, 0);
        unsigned artifactConv [2][2] = {0};                              //table for all artifacual conversions
        unsigned normalConv [2][2] = {0};                                //table for all non-artifactual conversions
        if (!wrapDoAll(artifactConv, normalConv, InsertCounts, bamFile, options))
            return 1;
        if (!wrapOutputInserts(InsertCounts, options))
            return 1;
        if (!wrapOutputArtifacts(artifactConv, normalConv, options))
            return 1;
    }
    else
    {
        if (options.insDist)                                            //Only InsertSize Distribution
        {
        TInsertDistr InsertCounts = "";
        if (!wrapCountInsertSize(InsertCounts, bamFile, options))
            return 1;
        if (!wrapOutputInserts(InsertCounts, options))
            return 1;
        }
        else if (options.conv)                                          //Only CCG > CAG or CGG > CTG Artifact check
        {
            FaiIndex faiIndex;
            if (!loadRefIdx(faiIndex, toCString(options.refPath)))
                return 1;
            unsigned artifactConv [2][2] = {0};
            unsigned normalConv [2][2] = {0};
            if (!getArtifactCount(artifactConv, normalConv, bamFile, faiIndex, options))
                return 1;
            if (!wrapOutputArtifacts(artifactConv, normalConv, options))
                return 1;
        }
    }
    return 0;
}