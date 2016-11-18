//Author: Sebastian Roskosch <Sebastian.Roskosch[at]bihealth.de>
#include <seqan/arg_parse.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

using namespace seqan;

/////////////////////Typedefs////////////////////////
//String holding the number of inserts of each length. Index 1 holds the number
//of inserts with length 1, index 2 holds the number of inserts with length 2...
//index 0 holds the number of segments without a mapped partner or the
//information is not available
typedef String<unsigned> TInsertDistr;

struct ProgramOptions //Struct holding all program options.
{
    CharString inPath;
    CharString refPath;
    CharString outPathInserts;
    CharString outPathArtifacts;
    bool insDist = false;
    int maxInsert;
    unsigned minMapQ;
    bool conv = false;
    unsigned verbosity = 1;
};
/////////////////////Parsing functions////////////////////////
//Parse the command line and/or display help message.
ArgumentParser::ParseResult parseCommandLine(ProgramOptions & options,
                                                    int argc,
                                                    char const ** argv)
{
    ArgumentParser parser("BAMQC");
    setShortDescription(parser, "Simple quality-control for (single-sample) BAM-files.");
    addUsageLine(parser, "BAM_FILE [OPTIONS]");
    setDate(parser, __DATE__);
    setVersion(parser, "1.0.0");
    
    addSection(parser, "I/O Options");
    ArgParseArgument fileArg(ArgParseArgument::INPUT_FILE, "FILE", false);
    setValidValues(fileArg, "bam sam");
    addArgument(parser, fileArg);

    addOption(parser, seqan::ArgParseOption(
    "r", "reference", "Path to reference genome. Required for C>A/G>T-Artifact-check.",
    seqan::ArgParseArgument::INPUT_FILE, "IN"));
    setValidValues(parser, "reference", "fasta fa fastq fq fasta.gz fa.gz fastq.gz fq.gz fasta.bz2 fa.bz2 fastq.bz2 fq.bz2");

    addOption(parser, seqan::ArgParseOption(
    "oi", "output-file-inserts", "Path to output file for the insert-size distribution.",
    seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));

    addOption(parser, seqan::ArgParseOption(
    "oc", "output-file-conversions", "Path to output file for the C>A/G>T-Artifact-check.",
    seqan::ArgParseArgument::OUTPUT_FILE, "OUT"));

    addSection(parser, "General Options");
    addOption(parser, seqan::ArgParseOption(
    "mmq", "min-mapq", "Minimum mapping quality.",
    seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "min-mapq", "25");
    setMinValue(parser, "min-mapq", "0");
    setMaxValue(parser, "min-mapq", "244");

    addOption(parser, seqan::ArgParseOption("v0", "no-verbosity", "Disable parameter feedback."));

    addSection(parser, "Insert-size-distribution Options");
    addOption(parser, seqan::ArgParseOption(
              "i", "insert-size-distribution",
              "Counts the insert size of each valid read-pair. Output to standard output if -oi with path is not specified."));

    addOption(parser, seqan::ArgParseOption(
    "m", "max-insert", "Maximum insert size. Sizes above will be ignored.",
    seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(parser, "max-insert", "1000");
    setMinValue(parser, "max-insert", "100");

    addSection(parser, "C>A/G>T-Artifact Options");
    addOption(parser, seqan::ArgParseOption(
              "c", "conversion-artifact",
              "Perform check for C>A/G>T artifacts induced during sample preparation (Costello et al. (2013))."
              "Requires reference genome. Output to standard output if -oc with path is not specified."));

    addTextSection(parser, "Examples");
    addListItem(parser,
            "\\fBBAMQC\\fP \\fBfile.bam\\fP \\fB-i\\fP",
            "Calculate insert-size distribution of mapped reads in \"file.bam\" and print output to standard output.");
    addListItem(parser,
            "\\fBBAMQC\\fP \\fBfile.bam\\fP \\fB-r reference.fa\\fP \\fB-i\\fP \\fB-c\\fP",
            "Check if reads in \"file.bam\" show an increased ammount of strand specific C>A/G>T conversions. and "
            "print output to standard output.");
     addListItem(parser,
        "\\fBBAMQC\\fP \\fBfile.bam\\fP \\fB-r reference.fa\\fP \\fB-oi inserts.txt\\fP \\fB-oc conversions.txt\\fP "
        "\\fB-mmq 30\\fP \\fB-m 500\\fP",
        "Perform both of the above checks in one run, but only consider alignments with mapping-quality of at least phred 30 (both checks) "
        "and a maximum considered insert-size of 500 (insert-size calculation). Output is saved in \"inserts.txt\" and \"conversions.txt.\"");

    //Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    //Terminate on error
    if (res != ArgumentParser::PARSE_OK)
        return res;
    //Get path of input BAM.
    getArgumentValue(options.inPath, parser, 0);
    getOptionValue(options.refPath, parser, "reference");
    getOptionValue(options.outPathInserts, parser, "output-file-inserts");
    getOptionValue(options.outPathArtifacts, parser, "output-file-conversions");
    options.insDist = isSet(parser, "insert-size-distribution");
    getOptionValue(options.maxInsert, parser, "max-insert");
    getOptionValue(options.minMapQ, parser, "min-mapq");
    options.conv = isSet(parser, "conversion-artifact");
    options.verbosity = !isSet(parser, "no-verbosity");
    return ArgumentParser::PARSE_OK;
}
//Check parameters for consistency. Return 1 on inconsistencies and 0 on pass.
inline int inputCheck(ProgramOptions & options)
{
    if (!empty(options.outPathInserts))
        options.insDist = true;
    if (!empty(options.outPathArtifacts))
        options.conv = true;
    if (!(options.insDist || options.conv))
    {
        std::cerr << "Error: No checks selected. Nothing to be done. Terminating.\n";
        return 1;
    }
    //check if both or none of conversion-flags and reference genome are given.
    if (options.conv && empty(options.refPath))
    {
        std::cerr << "Error: Missing reference genome for C>A/G>T artifact-check. Terminating.\n";
        return 1;
    }
    else if (!options.conv && !empty(options.refPath))
    {
        std::cerr << "Error: Reference genome given, but no required (consider setting the -i flag or giving a path "
        "for the output using -oc option). Terminating.\n";
        return 1;
    }
    return 0; //all go
}
//Return how the parameters are interpreted
inline void feedBack(const ProgramOptions & options)
{
    if (options.verbosity == 0) return;
    std::cout << "Parameters as interpreted:" << std::endl
              << "BAM-File: " << options.inPath << std::endl;
    if (!empty(options.refPath))
        std::cout << "Reference Genome: " << options.refPath << std::endl;
    std::cout << "Determine Insert-Size Distribution: ";
    if (options.insDist)
    {
        std::cout << "Yes" << std::endl
                  << "Output for Insert-Size Distribution: ";
        if (empty(options.outPathInserts))
            std::cout << "Standard Output" <<std::endl;
        else
            std::cout << options.outPathInserts << std::endl;
    }
    else
    {
        std::cout << "No" << std::endl;
    }
    std::cout << "Perform C>A/G>T-Conversion Check: ";
    if (options.conv)
    {
        std::cout << "Yes" << std::endl
                  << "Output for Conversion Check: ";
        if (empty(options.outPathArtifacts))
            std::cout << "Standard Output" <<std::endl;
        else
            std::cout << options.outPathArtifacts << std::endl;
    }
    else
    {
        std::cout << "No" << std::endl;
    }
    if (options.insDist)
        std::cout << "Maximum Considered Insert-Size: " << options.maxInsert << std::endl;
    std::cout << "Minimum Mapping-Quality: " << options.minMapQ << std::endl
              << "Verbosity: " << options.verbosity << std::endl
              << std::endl
              << "Perfoming selected tasks..." << std::endl
              << std::endl;
}
/////////////////////Input-file functions////////////////////////
//Load BAM-file.
inline bool loadBAM(BamFileIn & bamFile, const CharString & bamFileName)
{
    if (!open(bamFile, toCString(bamFileName)))
    {
        std::cerr << "ERROR: Could not open " << bamFileName << std::endl;
        return false;
    }
    return true;
}
//Load index of reference genome. If not available, build it on the fly.
inline bool loadRefIdx(FaiIndex & faiIndex, const CharString & refFileName)
{
    if (!open(faiIndex, toCString(refFileName)))
    {
        if (!build(faiIndex, toCString(refFileName)))
        {
            std::cerr << "ERROR: Index could not be loaded or built.\n";
            return false;
        }
        if (!save(faiIndex))
        {
            std::cerr << "ERROR: Index could not be written do disk.\n";
            return false;
        }
    }
    return true;
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
inline void formatStats(std::stringstream & out, const TInsertDistr & counts, const Pair<unsigned, unsigned> & firstLast)
{
    if (firstLast.i1 == 0 && firstLast.i2 == 0)
    {
        out << "No valid inserts detected.";
        return;
    }
    unsigned outputLength = firstLast.i2 - firstLast.i1 + 1;
    reserve(out, outputLength * 10);
    for (unsigned i = firstLast.i1; i <= firstLast.i2; ++i)
    {
        out << i << '\t' << counts[i] << '\n';
    }
}
inline void formatArtifacts(std::stringstream & out, unsigned (& artifactConv) [2][2], unsigned (& normalConv) [2][2])
{
    unsigned hits = artifactConv[0][0] + artifactConv[0][1] + artifactConv[1][0] + artifactConv[1][1];
    unsigned nonHits = normalConv[0][0] + normalConv[0][1] + normalConv[1][0] + normalConv[1][1];
    out << "Artifact-like Conversions: " << hits << " total" << std::endl
        << "\tForward\tReverse" << std::endl
        << "1st\t" << artifactConv[1][0] << "\t" << artifactConv[1][1] << std::endl
        << "2nd\t" << artifactConv[0][0] << "\t" << artifactConv[0][1] << std::endl << std::endl
        << "Other Conversions: " << nonHits << " total" << std::endl
        << "\tForward\tReverse" << std::endl
        << "1st\t" << normalConv[1][0] << "\t" << normalConv[1][1] << std::endl
        << "2nd\t" << normalConv[0][0] << "\t" << normalConv[0][1] << std::endl << std::endl
        << "Fraction of artifact-like conversions: " << (double)hits / double(hits + nonHits) << std::endl;
}
//Write stats to file or std::out
inline bool writeStats(const std::stringstream & out, const CharString & outPath)
{
    if (empty(outPath))
    {
        std::cout << "\n" << out.str();
        return true;
    }
    else
    {
        std::ofstream outStream;
        outStream.open(toCString(outPath));
        if (outStream)
        {
            outStream << out.str();
            outStream.close();
            std::cout << "Output written to " << outPath << std::endl;
            return true;
        }
        else
        {
            std::cerr << "Error while writing output-file.\n";
        }
    }
    return false;
}
//Wrapper for calling getFirstLast, formatStats and writeStats (=Wrtingin insert distribution to file)
inline bool wrapOutputInserts (const TInsertDistr & counts, const ProgramOptions & options)
{
    Pair<unsigned, unsigned> firstLast = getFirstLast(counts); //get borders of distribution for clean output
    std::stringstream out;
    formatStats(out, counts, firstLast);
    if (!writeStats(out, options.outPathInserts))
        return false;
    else return true;
}
//Wrapper for writing the conversions to file
inline bool wrapOutputArtifacts (unsigned (& artifactConv) [2][2],
                                 unsigned (& normalConv) [2][2],
                                 const ProgramOptions & options)
{
    std::stringstream out;
    formatArtifacts(out, artifactConv, normalConv);
    if (!writeStats(out, options.outPathArtifacts))
        return false;
    else return true;
}