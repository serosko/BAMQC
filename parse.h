#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>

using namespace seqan;

struct ProgramOptions //Struct holding all program options.
{
    CharString inPath;
};
//Parse the command line and/or display help message.
ArgumentParser::ParseResult parseCommandLine(ProgramOptions & options,
                                                    int argc,
                                                    char const ** argv)
{
    ArgumentParser parser("BAMQC");
    setShortDescription(parser, "Simple quality-control for BAM-files.");
    addUsageLine(parser, "BAM_FILE1 [BAM_FILE2] [OPTIONS] [-o Output_Prefix]");
    setDate(parser, __DATE__);
    setVersion(parser, "0.0.1");
    ArgParseArgument fileArg(ArgParseArgument::INPUT_FILE, "FILES", true);
    setValidValues(fileArg, "bam");
    addArgument(parser, fileArg);
    //Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    //Terminate on error
    if (res != ArgumentParser::PARSE_OK)
        return res;
    //Get path of input BAM.
    getArgumentValue(options.inPath, parser, 0);
    return ArgumentParser::PARSE_OK;
}
//Load BAM-file.
int loadBAM(BamFileIn& bamFile, const CharString bamFileName)
{
    if (!open(bamFile, toCString(bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
        return 1;
    }
    return 0;
}