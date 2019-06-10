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

    addOption(parser, ArgParseOption("a", "dict", "Path to the dictionary", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "dict");

    addOption(parser, ArgParseOption("f", "foundBarcodes", "Path to the flexbar result", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "foundBarcodes");
//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

    addOption(parser, ArgParseOption("bL", "barcodeL", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));
    setRequired(parser, "barcodeL");
    addOption(parser, ArgParseOption("bS", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString dictPath, flexPath, bamPath, outputPath;
    int batchSize1 = 100000, barcodeLength;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(dictPath, parser, "dict");
    getOptionValue(flexPath, parser, "foundBarcodes");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(barcodeLength, parser, "barcodeL");
    getOptionValue(batchSize1, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");


    //Load umi + barcorde dictionary
    StringSet<CharString> idsDict;
    StringSet<Dna5String> umi;
    SeqFileIn seqFileInReads(toCString(dictPath));
    readRecords(idsDict, umi, seqFileInReads);


    //prepare dictionary
    std::map<CharString,Dna5String> myMap;
    int m = 0;
    while(m < length(umi)){
        myMap[idsDict[m]] = umi[m];
        ++m;
    }

    //prepare flex result fasta
    StringSet<CharString> idsReads;
    SeqFileIn seqFileInFlex(toCString(flexPath));
    {
        StringSet<Dna5String> reads;
        readRecords(idsReads, reads, seqFileInFlex);
    }

    int f = 0;
    std::map<CharString,std::pair<Dna5String, Dna5String> > umiBarMap;
    std::cout << "N idsReads " << length(idsReads) << "\n";
    while(f < length(idsReads)){
        CharString & id = idsReads[f];

        if(verbose){
//             std::cout << "BamName " << toCString(record.qName) << "\n";
            std::cout << toCString(id) << "\n";
            std::cout << length(id) << "\n";
        }
        Finder<CharString> finder(id);
        Pattern<CharString, Horspool> pattern("_Flexbar_removal_");
        find(finder, pattern);
        int begin = endPosition(finder);
        std::cout << "begin: " << begin << "\n";




        if(begin > 0){
            int endName = beginPosition(finder);
            CharString readName = prefix(id, endName);
            std::cout << toCString(readName) << "\n";
            Finder<CharString> finder2(id);
            Pattern<CharString, Horspool> pattern2(" revcomp");
            find(finder2, pattern2);
            int end = beginPosition(finder2);
            std::cout << "end: " << end << "\n";
            bool doRC = end > 0;
            CharString acc;
            // if a reverse complement was found
            if(doRC)
                acc = infix(id, begin, end);
            else
                acc = suffix(id, begin);

            if(verbose)
                std::cout << "Search acc: " << toCString(acc) << "\n";

            // look up dictionary to find corresponding umi and barcode
            auto search = myMap.find(acc);
            if (search != myMap.end()){
                if(verbose)
                    std::cout << "Found " << search->first << " " << search->second << '\n';
            } else {
                std::cout << "Not found\n";
                std::cout << "Accession from bam file not found in Dictionary\n";
//                     writeRecord(bamFileOut, record);
//                     continue;
                exit(0);
            }
            Dna5String mySeq = search->second;
            Dna5String mybar = prefix(search->second, barcodeLength);
            Dna5String myumi = suffix(search->second, barcodeLength);
//             myMap[idsDict[m]] = umi[m];
            umiBarMap[readName] = std::make_pair(mybar, myumi);
        }
        else
        {
            umiBarMap[id] = std::make_pair("", "");
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
    BamFileOut bamFileOut(context(bamFileIn), toCString(outputPath));
    try
    {
        readHeader(header, bamFileIn);
        writeHeader(bamFileOut, header);
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
                    std::cout << "Found " << search->first << " " << toCString(search->second.first) << '\n';
            } else {
                std::cout << "Not found\n";
                std::cout << "Accession from bam file not found in Dictionary\n";
                writeRecord(bamFileOut, record);
                continue;
            }
                //TODO compare accessions
            //umi and barcode are found
            if(length(search->second.first) > 0){
                Dna5String mybar = search->second.first;
                Dna5String myumi = search->second.second;
                //add umi and barcode as tags to bam
                BamTagsDict tags(record.tags);
                setTagValue(tags, "CR", mybar);
                setTagValue(tags, "UR", myumi);
                record.tags = host(tags);
                writeRecord(bamFileOut, record);
            }
            else
            {
                if(verbose)
                    std::cout << "Skip\n";
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



