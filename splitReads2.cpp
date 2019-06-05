#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include<iostream>
#include<fstream>
#include <seqan/find.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>

using namespace seqan;
using namespace std;





int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");

//     addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
//     setRequired(parser, "bam");

    addOption(parser, ArgParseOption("f", "foundBarcodes", "Path to the flexbar result", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "foundBarcodes");
//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output files prefix", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

//     addOption(parser, ArgParseOption("bL", "barcodeL", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));
//     setRequired(parser, "barcodeL");
    addOption(parser, ArgParseOption("bS", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString dictPath, flexPath, bamPath, outputPath;
    int batchSize1 = 100000, barcodeLength;

//     getOptionValue(bamPath, parser, "bam");
    getOptionValue(flexPath, parser, "foundBarcodes");
    getOptionValue(outputPath, parser, "output");
//     getOptionValue(barcodeLength, parser, "barcodeL");
    getOptionValue(batchSize1, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");



    //prepare flex result fasta
    StringSet<CharString> idsReads;
    SeqFileIn seqFileInFlex(toCString(flexPath));


    CharString outputRevcomp = outputPath;
    outputPath += "_left_tail_trimmed.fasta";
    outputRevcomp += "_right_tail_trimmed.fasta";

    SeqFileOut seqFileOut(toCString(outputPath));
    SeqFileOut seqFileOutRevcomp(toCString(outputRevcomp));

    while (!atEnd(seqFileInFlex))
    {
        Dna5String read;
        CharString id;
        try
        {

            readRecord(id, read, seqFileInFlex);

            Finder<CharString> finder(id);
            Pattern<CharString, Horspool> pattern("_Flexbar_removal_");
            find(finder, pattern);
            int test = beginPosition(finder);
            if(test > 0){
                Finder<CharString> finder(id);
                Pattern<CharString, Horspool> pattern(" revcomp");
                find(finder, pattern);
                int end = beginPosition(finder);
    //             std::cout << "end: " << end << "\n";
                bool doRC = end > 0;
                if (doRC)
                    writeRecord(seqFileOutRevcomp, id, read);
                else
                    writeRecord(seqFileOut, id, read);
            }
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    close(seqFileInFlex);


    std::cout << "Finished!\n";

    return 0;
}



