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



void writeRead(SeqFileOut & seqFileOut, string & last_id, Dna5String & end1, Dna5String & end2, bool rc, bool verbose = false)
{
    if(rc){
        Dna5StringReverseComplement rcend1(end1);
        Dna5StringReverseComplement rcend2(end2);

        if(verbose){
            std::cout << "Use reverse Complement of reads to restore original orientation.\n";
            std::cout << "before: " << end1 << "\n";
            std::cout << "before2: " << end2 << "\n";
        }
        Dna5String tmpend1 = rcend1;
        end1 = rcend2;
        end2 = tmpend1;
        if(verbose){
            std::cout << "after:  " << end1 << "\n";
            std::cout << "after2:  " << end2 << "\n";
        }
    }
    else
    {
        if(verbose)
            std::cout << "Is in original orientation\n";
    }

    CharString idend1 = last_id;
    CharString idend2 = last_id;
    idend1 += "_end1";
    idend2 += "_end2";

    writeRecord(seqFileOut, idend1, end1);
    writeRecord(seqFileOut, idend2, end2);
}


bool checkIfSameOrientation(Dna5String & originalRead, Dna5String & bamRead){

    bool reverseComplement;
    Dna5StringReverseComplement rcBamRead(bamRead);
//     appendValue(rcReads, myModifier);
    for(int i = 0; i < length(bamRead); ++i)
    {
//         std::cout << i << ": " << originalRead[i] << "\t" << bamRead[i] << "\t" << rcBamRead[i] << "\n";

        if(originalRead[i] != bamRead[i] && originalRead[i] != rcBamRead[i]){
            std::cout << "Bam and fasta reads are not the same.\n";
            exit(0);
        }
        if(originalRead[i] != bamRead[i]){
            reverseComplement = true;
            break;
        }
        if(originalRead[i] != rcBamRead[i]){
            reverseComplement = false;
            break;
        }
    }
/*
    if(reverseComplement)
        std::cout << "reverseComplement\n";
    else
        std::cout << "not\n";*/

    return reverseComplement;
}

Dna5String & getOriginalRead(auto & readMap, string & readID)
{
    CharString creadID = readID;
    auto search = readMap.find(creadID);
    if (search == readMap.end())
    {
        std::cout << "Did not find " << readID << " in fasta file." << "\n";
        exit(0);
    }
    Dna5String & originalRead = search->second;

    return originalRead;
}

int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");

    addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "bam");
/*
    addOption(parser, ArgParseOption("p", "first", "First p reads", ArgParseArgument::INTEGER, "INT"));
    hideOption(getOption(parser, "first"));

    addOption(parser, ArgParseOption("mt", "threshold", "Number of matches required to define start or end of the cDNA part", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));*/

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString bamPath;
    int batchSize1 = 100000;

    int barcode_umi_length = 30;
    int threshold = 5;
    int first = 9999999;

    getOptionValue(bamPath, parser, "bam");
//     getOptionValue(barcodeLength, parser, "barcodeL");
//     getOptionValue(first, parser, "first");
//     getOptionValue(threshold, parser, "threshold");
//     bool verbose = isSet(parser, "verbose");
    
    
        //Read first time
    std::map<CharString, CharString> readMap;
    {
        StringSet<CharString> readIDs;
        StringSet<CharString> reads;
        BamAlignmentRecord record;

        BamFileIn bamFileIn;
        if (!open(bamFileIn, toCString(bamPath)))
        {
            std::cerr << "ERROR: could not open input file " << bamPath << ".\n";
            return 1;
        }
        BamHeader header;
        readHeader(header, bamFileIn);

        try
        {
            while (!atEnd(bamFileIn))
            {
                readRecord(record, bamFileIn);
                appendValue(readIDs, record.qName);
                appendValue(reads, record.seq);
            }
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
            exit(0);
        }
        //prepare dictionary
        int m = 0;
        while(m < length(readIDs)){
            if(length(reads[m]) > 0)
                readMap[readIDs[m]] = reads[m];
            ++m;
        }
        close(bamFileIn);
    }

    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamPath)))
    {
        std::cerr << "ERROR: Could not open " << bamPath << std::endl;
        return 1;
    }
    // Open output file, BamFileOut accepts also an ostream and a format tag.
    BamFileOut bamFileOut(context(bamFileIn), std::cout, Sam());

    try
    {
        // Copy header.
        BamHeader header;
        readHeader(header, bamFileIn);
        writeHeader(bamFileOut, header);

        // Copy records.
        BamAlignmentRecord record;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            CharString id = record.qName;
            if(length(record.seq) == 0){
                auto search = readMap.find(id);
                if (search == readMap.end())
                {
                    std::cerr << "Did not find " << id << " in bam file. There is no sequence for this readid denoted in the bam file" << "\n";
                    exit(0);
                }
                CharString seq = search->second;
                record.seq = seq;
            }
            writeRecord(bamFileOut, record);
        }
    }
    catch (Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}



