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

    addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "bam");

    addOption(parser, ArgParseOption("f", "foundBarcodes", "Path to the flexbar result", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "foundBarcodes");
//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output files prefix", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

    addOption(parser, ArgParseOption("v", "verbose", "Check for if Sequence is located on the correct side according to the orientation"));

//     addOption(parser, ArgParseOption("bL", "barcodeL", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));
//     setRequired(parser, "barcodeL");
    addOption(parser, ArgParseOption("bS", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString dictPath, flexPath, bamPath, outputPath;
    int batchSize1 = 100000, barcodeLength;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(flexPath, parser, "foundBarcodes");
    getOptionValue(outputPath, parser, "output");
//     getOptionValue(barcodeLength, parser, "barcodeL");
    getOptionValue(batchSize1, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");



    //prepare flex result fasta
    StringSet<CharString> idsReads;
    SeqFileIn seqFileInFlex(toCString(flexPath));
    {
        StringSet<Dna5String> reads;
        readRecords(idsReads, reads, seqFileInFlex);
    }

    int f = 0;
    std::map<CharString, bool> umiBarMap;
    std::cout << "N idsReads " << length(idsReads) << "\n";
    while(f < length(idsReads)){
        CharString & id = idsReads[f];

        if(verbose){
            std::cout << toCString(id) << "\n";
            std::cout << length(id) << "\n";
        }
        Finder<CharString> finder(id);
        Pattern<CharString, Horspool> pattern("_Flexbar_removal_");
        find(finder, pattern);
        int begin = endPosition(finder);
        if(verbose)
            std::cout << "begin: " << begin << "\n";

        if(begin > 0){
            int endName = beginPosition(finder);
            CharString readName = prefix(id, endName);
            if(verbose)
                std::cout << toCString(readName) << "\n";
            Finder<CharString> finder2(id);
            Pattern<CharString, Horspool> pattern2(" revcomp");
            find(finder2, pattern2);
            int end = beginPosition(finder2);
            if(verbose)
                std::cout << "end: " << end << "\n";
            bool doRC = end > 0;

            umiBarMap[readName] = doRC;
        }
        else
        {
            umiBarMap[id] = false;
        }
        ++f;



    }

    //prepare Bam
    BamFileIn bamFileIn;
    BamHeader header;
    if (!open(bamFileIn, toCString(bamPath)))
    {
        std::cerr << "ERROR: could not open input file " << bamPath << ".\n";
        return 1;
    }

    //write Header
    CharString outputPathReverse = outputPath;
    outputPath += "_left_tail_trimmed.bam";
    outputPathReverse += "_right_tail_trimmed.bam";
    BamFileOut bamFileOut(context(bamFileIn), toCString(outputPath));
    BamFileOut bamFileOutRevcomp(context(bamFileIn), toCString(outputPathReverse));
    try
    {
        readHeader(header, bamFileIn);
        writeHeader(bamFileOut, header);
        writeHeader(bamFileOutRevcomp, header);
    }
    catch (IOError const & e)
    {
        std::cerr << "ERROR: could not copy header. " << e.what() << "\n";
    }

    //modify bam (Add UMI and Barcode (determined from flexbar) to each Alignment record)
    BamAlignmentRecord record;
    CharString id;
    Dna5String seq;
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
            if(verbose){
                std::cout << "Accession in Bam " << toCString(record.qName) << "\n";
            }

            auto search = umiBarMap.find(record.qName);
            if (search != umiBarMap.end()) {
                if(verbose)
                    std::cout << "Found " << search->first << " " << search->second << '\n';
            } else {
                std::cout << "Not found\n";
                std::cout << "Accession from bam file not found in Flexbar result file\n";
//                 writeRecord(bamFileOut, record);
//                 continue;
                exit(0);
            }
                //TODO compare accessions
            //umi and barcode are found
            if(search->second){
                if(verbose)
                    std::cout << "RevComp\n";
                writeRecord(bamFileOutRevcomp, record);
            }
            else
            {
                writeRecord(bamFileOut, record);
            }

        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
    }

    close(seqFileInFlex);
    close(bamFileOut);


    std::cout << "Finished!\n";

    return 0;
}



