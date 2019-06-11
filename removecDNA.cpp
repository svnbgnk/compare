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

    addOption(parser, ArgParseOption("f", "fasta", "Path to the fasta with original reads to recover orientation", ArgParseArgument::INPUT_FILE, "IN"));

    addOption(parser, ArgParseOption("o", "output", "Path to fasta output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
    setRequired(parser, "output");

    addOption(parser, ArgParseOption("p", "first", "First p reads", ArgParseArgument::INTEGER, "INT"));
    hideOption(getOption(parser, "first"));

    addOption(parser, ArgParseOption("mt", "threshold", "Number of matches required to define start or end of the cDNA part", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString bamPath, outputPath;
    CharString fastaPath = "";
    int batchSize1 = 100000;

    int barcode_umi_length = 30;
    int threshold = 5;
    int first = 9999999;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(fastaPath, parser, "fasta");
    getOptionValue(outputPath, parser, "output");
//     getOptionValue(barcodeLength, parser, "barcodeL");
    getOptionValue(first, parser, "first");
    getOptionValue(threshold, parser, "threshold");
    bool verbose = isSet(parser, "verbose");

    bool checkOrientation = false;

    //Load original reads

    std::map<CharString,Dna5String> readMap;

    if(length(fastaPath) > 0){
        StringSet<CharString> readIDs;
        StringSet<Dna5String> reads;
        checkOrientation = true;
        SeqFileIn seqFileInReads(toCString(fastaPath));
        readRecords(readIDs, reads, seqFileInReads);
        //prepare dictionary
        int m = 0;
        while(m < length(readIDs)){
            readMap[readIDs[m]] = reads[m];
            ++m;
        }
    }

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
    Dna5String last_read;
    int best_score = 0;
    Dna5String best_end1;
    Dna5String best_end2;

    bool start = true;


    int y = 0;
    while (!atEnd(bamFileIn))
    {
        ++y;

        try
        {
            readRecord(record, bamFileIn);
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
            continue;
        }

        if (y > first)
            continue;

            Dna5String end1 = record.seq;
            Dna5String end2 = record.seq;
            CharString id = record.qName;

            if(start){
                last_id = toCString(id);
                last_read = record.seq;
                start = false;
            }


            if(length(record.seq) == 0){
                std::cout << toCString(id) << "\n";
                std::cout << "Read could not be loaded\n";
                continue;
            }

            if(verbose){
                std::cout << "Accession in Bam " << toCString(id) << "\n";
                std::cout << "Read length: " << length(record.seq) << ": " << record.seq << "\n";
            }



            String<CigarElement<> > cigar = record.cigar;
            //take front
            int pos = 0;
            //count amount of overall matches in cigar string
            int match_count = 0;
            bool found = false;
            for(int i = 0; i < length(cigar); ++i){

                if (cigar[i].operation == op_match)
                    match_count += cigar[i].count;

                //taake prefix of the read as soon as more than threshold many matches are found
                if (match_count > threshold && !found){
//                     std::cout << "found pos: " << pos << "\n";
                    end1 = prefix(end1, pos + 1);
                    found = true;
                }
                //calculate end position of prefix according to cigar string
                if(!found && (cigar[i].operation == op_match || cigar[i].operation == op_soft || cigar[i].operation == op_del))
                    pos += cigar[i].count;
            }

            int score = 0;


            if(!found){
                std::cout << "Take whole read since it was not aligned\n";
                writeRecord(seqFileOut, id, end2);
                continue;
            }

                pos = 0;
                int match_count2 = 0;
                bool found2 = false;
                for(int i = length(cigar) - 1; i >= 0 ; --i){

                    if (cigar[i].operation == op_match)
                        match_count2 += cigar[i].count;

                    //since the overall match count is know we can stop as soon as we determined the start position of the prefix
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
                    std::cout << "lenght end1: " << length(end1) << "\n";

                    std::cout << "lenght end2: " << length(end2) << "\n";

                    std::cout << "Score: " << (match_count) << "\n";
                }


                score = match_count;

        if(last_id.compare(toCString(id)) == 0){
            if(score > best_score){
                best_score = score;
                best_end1 = end1;
                best_end2 = end2;
            }
        }
        else
        {

            if(verbose){

                std::cout << "New Read Id, therefore write alignment with best score of: " << best_score << "\n";
                std::cout << "from read id: " << last_id << "\n";
            }

            bool rc = false;
            if(checkOrientation){

                Dna5String & originalRead = getOriginalRead(readMap, last_id);
//                 Dna5String bamRead = record.seq;
//                 std::cout << "ori: " << originalRead << "\n";

                rc = !checkIfSameOrientation(originalRead, last_read);
            }

            writeRead(seqFileOut, last_id, best_end1, best_end2, rc, verbose);

            best_score = score;
            best_end1 = end1;
            best_end2 = end2;
            last_id = toCString(id);
            last_read = record.seq;
        }

        if(atEnd(bamFileIn)){/*
            bool rc = false;
            if(checkOrientation){
                Dna5String & originalRead = getOriginalRead(readMap, last_id);
                std::cout << "ori: " << originalRead << "\n";
                rc = !checkIfSameOrientation(originalRead, record.seq);
            }*/
            writeRead(seqFileOut, last_id, best_end1, best_end2, false);
        }



    }
    close(seqFileOut);


    std::cout << "Finished!\n";

    return 0;
}



