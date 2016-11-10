#include <BAMQC.h>

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
    if (options.insDist)
    {
        TInsertDistr counts = "";
        wrapCountInsertSize(counts, bamFile, options);
        wrapOutput(counts, options);
    }
    if (options.catg)
    {
        FaiIndex faiIndex;
        loadRefIdx(faiIndex, options.refPath);

    }
    return 0;
}