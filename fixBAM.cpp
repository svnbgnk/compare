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
    
    CigarElement<> hard = seqan::CigarElement<>('H', 1);
    CigarElement<> soft = seqan::CigarElement<>('S', 1);
    auto op_hard = hard.operation;
    auto op_soft = soft.operation;
    
        //Read first time
    std::map<CharString, std::tuple<CharString, CharString, uint16_t, uint16_t> > readMap;
    {
        StringSet<CharString> readIDs;
        StringSet<CharString> reads;
        StringSet<CharString> quals;
        StringSet<String<CigarElement<> >> cigars;
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
                appendValue(quals, record.qual);
                appendValue(cigars, record.cigar);
            }
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
            exit(0);
        }
        //prepare dictionary
        int m = 0;
        //TODO RC and middle part is clipped
        while(m < length(readIDs)){
            
            if(length(reads[m]) > 0){
                int16_t len = length(cigars[m]) - 1;
                auto search = readMap.find(readIDs[m]);
                if (search != readMap.end())
                {
//                     std::cerr << "Found id again\n";
                    std::tuple<CharString, CharString, uint16_t, uint16_t> res = search->second;
                    CharString & seq = std::get<0>(res);
                    CharString & qual = std::get<1>(res);
                    uint16_t ofrontH = std::get<2>(res);
                    uint16_t oendH = std::get<3>(res);
                    uint16_t fchardClips = (op_hard == cigars[m][0].operation || op_soft == cigars[m][0].operation) ? cigars[m][0].count : 0;
                    uint16_t bchardClips = (op_hard == cigars[m][len].operation || op_soft == cigars[m][len].operation) ? cigars[m][len].count : 0;
                    auto newseq = seq;
                    auto newqual = qual;
                    if(fchardClips < ofrontH){
                        newseq = infix(newseq, 0, ofrontH - fchardClips);
                        newqual = infix(newqual, 0, ofrontH - fchardClips);
                        append(newseq, reads[m]);
                        append(newqual, qual[m]);
                        std::cerr << "m: " << m << "\t" << readIDs[m] << "\tbetter read found updating";
                        ofrontH = fchardClips;
                    }
                    if(bchardClips < oendH){
                        seq = infix(seq, length(seq) + oendH - bchardClips - 1, length(seq) - 1);
                        qual = infix(qual, length(seq) + oendH - bchardClips - 1, length(qual) - 1);
                        append(newseq, seq);
                        append(newqual, qual);
                        std::cerr << "m: " << m << "\t" << readIDs[m] << "\tbetter read found updating";
                        oendH = bchardClips;
                    }
                    search->second = std::make_tuple(newseq, quals[m], ofrontH, oendH);
                }else{
//                     std::cerr << "Found new id\n";
                    uint16_t fchardClips = (op_hard == cigars[m][0].operation || op_soft == cigars[m][0].operation) ? cigars[m][0].count : 0;
                    uint16_t bchardClips = (op_hard == cigars[m][len].operation || op_soft == cigars[m][len].operation) ? cigars[m][len].count : 0;
//                     std::cerr << "Clips\n";
                    readMap[readIDs[m]] = std::make_tuple(reads[m], quals[m], fchardClips, bchardClips);
//                     std::cerr << "fin\n";
                }
            }
            ++m;
        }
        close(bamFileIn);
    }

    std::cerr << "Finished Map\n";
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

                std::tuple<CharString, CharString, uint16_t, uint16_t> res = search->second;
                CharString & seq = std::get<0>(res);
                CharString & qual = std::get<1>(res);
                int ofrontH = std::get<2>(res);
                int oendH = std::get<3>(res);
                
                auto cigar = record.cigar;
                int len = length(cigar) - 1;
                int fchardClips = (op_hard == cigar[0].operation || op_soft == cigar[0].operation) ? cigar[0].count : 0;
                int bchardClips = (op_hard == cigar[len].operation || op_soft == cigar[len].operation) ? cigar[len].count : 0;
                
                record.seq = infix(seq, fchardClips - ofrontH, length(seq) - bchardClips + oendH - 1);
                record.qual = infix(qual, fchardClips - ofrontH, length(seq) - bchardClips + oendH - 1);
                std::cerr << "id: " << id << "\t" << length(seq) << ":" << length(record.seq) << "\t" << ofrontH << ":" << oendH << "\t"  << fchardClips << ":" << bchardClips << "\n" << toCString(seq) << "\n\n";
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



