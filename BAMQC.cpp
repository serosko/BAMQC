#include <catg.h>

int main(int argc, char const ** argv)
{
    ProgramOptions options;
    //Initialize parser
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options,
                                                              argc,
                                                              argv);
    if (res != seqan::ArgumentParser::PARSE_OK) 
        return res == seqan::ArgumentParser::PARSE_ERROR;               //Terminate on parsing errors
    if(inputCheck(options))                                             //Terminate if check of parameters fails.
        return 1;
    BamFileIn bamFile;                                                  //Prepare and load BAM-file
    loadBAM(bamFile, options.inPath);
    BamHeader header;                                                   //Read header to get to right position in file
    readHeader(header, bamFile);
    if (options.insDist)                                                //InsertSize Distribution
    {
        TInsertDistr InsertCounts = "";
        wrapCountInsertSize(InsertCounts, bamFile, options);
        wrapOutput(InsertCounts, options);
    }
    if (options.catg)                                                   // CCG > CAG or CGG > CTG Artifact check
    {
        FaiIndex faiIndex;
        loadRefIdx(faiIndex, options.refPath);
        unsigned oxoCounts = getArtifactCount(bamFile, faiIndex, options);
        std::cout << oxoCounts << std::endl;
    }
    return 0;
}