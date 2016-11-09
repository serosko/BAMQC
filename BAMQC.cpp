#include <BAMQC.h>

int main(int argc, char const ** argv)
{
    ProgramOptions options;
    //Initialize parser
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options,
                                                              argc,
                                                              argv);
    //Check if input was correctly parsed
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    BamFileIn bamFile;                              //Prepare and load BAM-file
    loadBAM(bamFile, options.inPath);

    BamHeader header;
    readHeader(header, bamFile);
//     typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
//     TBamContext const & bamContext = context(bamFile);
//     for (unsigned i = 0; i < length(contigNames(bamContext)); ++i)
//         std::cout << contigNames(bamContext)[i] << '\t'
//                   << contigLengths(bamContext)[i] << '\n';
    TInsertDistr counts = "";
    resize(counts, 500, 0);
    wrapCountInsertSize(counts, bamFile, options);
    wrapOutput(counts, options);
    return 0;
}