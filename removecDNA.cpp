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

    addOption(parser, ArgParseOption("o", "output", "Path to fasta output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("mt", "threshold", "Number of matches required to start or end cDNA part", ArgParseArgument::INTEGER, "INT"));
    
    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString bamPath, outputPath;
    int batchSize1 = 100000;
    
    int barcode_umi_length = 30;
    int threshold = 5;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(outputPath, parser, "output");
//     getOptionValue(barcodeLength, parser, "barcodeL");
    getOptionValue(threshold, parser, "threshold");
    bool verbose = isSet(parser, "verbose");

    //prepare Bam
    BamFileIn bamFileIn;
    BamHeader header;
    if (!open(bamFileIn, toCString(bamPath)))
    {
        std::cerr << "ERROR: could not open input file " << bamPath << ".\n";
        return 1;
    }
    readHeader(header, bamFileIn);

    SeqFileOut seqFileOut(toCString(outputPath));

    //modify bam (Add UMI and Barcode (determined from flexbar) to each Alignment record)
    CigarElement<> soft = seqan::CigarElement<>('S', 1);
    CigarElement<> match = seqan::CigarElement<>('M', 1);
    CigarElement<> del = seqan::CigarElement<>('D', 1);
    auto op_soft = soft.operation;
    auto op_match = match.operation;
    auto op_del = del.operation;
    BamAlignmentRecord record;
    string last_id;
    int best_score = 0;
    Dna5String best_end1;
    Dna5String best_end2;
    bool start = true;
    
    while (!atEnd(bamFileIn))
    {
        try
        {
            readRecord(record, bamFileIn);
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
            continue;
        }
            
            Dna5String end1 = record.seq;
            Dna5String end2 = record.seq;
            CharString id = record.qName;
            
            if(start){
                last_id = toCString(id);
                start = false;
            }
            
            if(verbose){
                std::cout << "Accession in Bam " << toCString(id) << "\n";
                std::cout << "Read length: " << length(end1) /*<< ": " << record.seq*/ << "\n";
            }
            
            String<CigarElement<> > cigar = record.cigar;
            //take front
            int pos = 0;
            int match_count = 0;
            bool found = false;
            for(int i = 0; i < length(cigar); ++i){
//                  std::cout << "C: " << cigar[i].count << cigar[i].operation << "\n";
//                 std::cout << match_count << "\n";
//                 std::cout << pos << "\n";
//                 std::cout << (cigar[i].operation == match.operation) << "\n";
                if (cigar[i].operation == op_match)
                    match_count += cigar[i].count;
                
                if (match_count > threshold){
//                     std::cout << "found pos: " << pos << "\n";
                    end1 = prefix(end1, pos + 1);
                    
                    found = true;
                    break;
                }
                if(cigar[i].operation == op_match || cigar[i].operation == op_soft || cigar[i].operation == op_del)
                    pos += cigar[i].count;
            }
            if(!found){
                std::cout << "Take whole read since it was not aligned\n";
                writeRecord(seqFileOut, id, end2);
                continue;
            }
                
//             std::cout << "final Match count: " << match_count << " found: " << found << "\n";
//             std::cout << "Cut end" << "\n";
            pos = 0;
            int match_count2 = 0;
            bool found2 = false;
            for(int i = length(cigar) - 1; i >= 0 ; --i){
//                 std::cout << "C: " << cigar[i].count << cigar[i].operation << "\n";
//                 std::cout << pos << "\n";
                if (cigar[i].operation == op_match)
                    match_count2 += cigar[i].count;
                
                if (match_count2 > threshold){
//                     std::cout << "found pos: " << pos << "\n";
                    end2 = suffix(end2, length(end2) - pos - 1);
                    found2 = true;
                    break;
                }
                
                if(cigar[i].operation == op_match || cigar[i].operation == op_soft || cigar[i].operation == op_del)
                    pos += cigar[i].count;
            }
            
            if(verbose){
                std::cout << "lenght end1: " << length(end1) << "\tScore: " << match_count << "\n";
                
                std::cout << "lenght end2: " << length(end2) << "\tScore: " << match_count2 << "\n";
                
                std::cout << "Score: " << (match_count + match_count2) << "\n";
            }
            
            
        int score = match_count + match_count2;
        
        if(last_id.compare(toCString(id)) == 0){
            if(score > best_score){
                best_score = score;
                best_end1 = end1;
                best_end2 = end2;
            }
        }
        else
        {
            CharString idend1 = last_id;
            CharString idend2 = last_id;
            idend1 += "_end1";
            idend2 += "_end2";
            
            if(verbose){
                std::cout << "Select this one: " << last_id << "\n";
                std::cout << "Use alignment with Score: " << best_score << "\n";
            }
            
            writeRecord(seqFileOut, best_end1, end1);
            writeRecord(seqFileOut, best_end2, end2);
            
            best_score = score;
            best_end1 = end1;
            best_end2 = end2;
            last_id = toCString(id);
        }
        
        if(atEnd(bamFileIn)){
            CharString idend1 = last_id;
            CharString idend2 = last_id;
            idend1 += "_end1";
            idend2 += "_end2";
            writeRecord(seqFileOut, best_end1, end1);
            writeRecord(seqFileOut, best_end2, end2);
        }

    }


    std::cout << "Finished!\n";

    return 0;
}



