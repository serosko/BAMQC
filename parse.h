#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>


using namespace seqan;

//String holding the number of inserts of each length. Index 1 holds the number
//of inserts with length 1, index 2 holds the number of inserts with length 2...
//index 0 holds the number of segments without a mapped partner or the i
//information is not available
typedef String<unsigned> TInsertDistr;

struct ProgramOptions //Struct holding all program options.
{
    CharString inPath;
    CharString outPath;
    int maxInsert;
    unsigned minMapQ;
};
//Parse the command line and/or display help message.
ArgumentParser::ParseResult parseCommandLine(ProgramOptions & options,
                                                    int argc,
                                                    char const ** argv)
{
    ArgumentParser parser("BAMQC");
    setShortDescription(parser, "Simple quality-control for BAM-files.");
    addUsageLine(parser, "BAM_FILE1 [BAM_FILE2 ...] [OPTIONS] [-O Output_File]");
    setDate(parser, __DATE__);
    setVersion(parser, "0.0.1");
    
    ArgParseArgument fileArg(ArgParseArgument::INPUT_FILE, "FILES", true);
    setValidValues(fileArg, "bam");
    addArgument(parser, fileArg);
    
    addOption(parser, seqan::ArgParseOption(
    "O", "output-file", "Path to the output file.",
    seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));
    
    addOption(parser, seqan::ArgParseOption(
    "m", "max-insert", "Maximum insert size.",
    seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "max-insert", "10000");
    setMinValue(parser, "max-insert", "100");
    
        addOption(parser, seqan::ArgParseOption(
    "mmq", "min-mapq", "Minimum mapping quality.",
    seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "min-mapq", "25");
    setMinValue(parser, "min-mapq", "0");
    setMaxValue(parser, "min-mapq", "244");
    
    //Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    //Terminate on error
    if (res != ArgumentParser::PARSE_OK)
        return res;
    //Get path of input BAM.
    getArgumentValue(options.inPath, parser, 0);
    getOptionValue(options.outPath, parser, "output-file");
    getOptionValue(options.maxInsert, parser, "max-insert");
    getOptionValue(options.minMapQ, parser, "min-mapq");
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
/////////////////////Output functions////////////////////////
//Dertermine first an last non-zero insert size for cleaner output
Pair<unsigned, unsigned> getFirstLast (const TInsertDistr & counts)
{
    Pair<unsigned, unsigned> firstLast = Pair<unsigned, unsigned> (0,0);
    for(unsigned i = 1; i < length(counts); ++i)
    {
        if (counts[i] != 0)
        {
            firstLast.i1 = i;
            break;
        }
    }
    for(unsigned j = length(counts) - 1; j >= firstLast.i1; --j)
    {
        if (counts[j] != 0)
        {
            firstLast.i2 = j;
            break;
        }
    }
    return firstLast;
}
//Format output stats
void formatStats(std::stringstream & out, const TInsertDistr & counts, const Pair<unsigned, unsigned> & firstLast)
{
    unsigned outputLength = firstLast.i2 - firstLast.i1 + 1;
    reserve(out, outputLength * 10);
    for (unsigned i = firstLast.i1; i <= firstLast.i2; ++i)
    {
        out << i << '\t' << counts[i] << '\n';    //Todo fix bug
    }
}
//Write stats to file or std::out
void writeStats(const std::stringstream & out, const ProgramOptions & options)
{
    if (empty(options.outPath))
    {
        std::cout << out.str();
    }
    else
    {
        std::ofstream outStream;
        outStream.open(toCString(options.outPath));
        if (outStream)
        {
            outStream << out.str();
            outStream.close();
            std::cout << "Output written to " << options.outPath << std::endl;
        }
        else
        {
            std::cout << "Error while writing output-file.\n";
        }
    }
}
//Wrapper for calling getFirstLast, formatStats and writeStats
void wrapOutput (const TInsertDistr & counts, const ProgramOptions & options)
{
    Pair<unsigned, unsigned> firstLast = getFirstLast(counts); //get borders of distribution for clean output
    std::stringstream out;
    formatStats(out, counts, firstLast);
    writeStats(out, options);
}