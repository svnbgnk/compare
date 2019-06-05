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
#include <string>

using namespace seqan;
using namespace std;

//TODO star did not write down unmapped reads to to rerun script

struct Stats
{
    uint32_t reads_count = 0;
    
    uint32_t ambiguous_count = 0;
    uint32_t ambiguous1_count = 0;
    uint32_t ambiguous2_count = 0;
    
    uint32_t unmapped_count = 0;
    uint32_t unmapped1_count = 0;
    uint32_t unmapped2_count = 0;
    uint32_t similar_count = 0;
    std::vector<uint32_t> differences; 
    
    void print()
    {
        std::cout << "Overall reads: " << 2*reads_count << "\n";
        
        std::cout << "Global ambiguous: " << ambiguous_count << "\n";
        std::cout << "Global unmapped: " << unmapped_count << "\n";
        
        std::cout << "1 ambiguous: " << ambiguous1_count << "\n";
        std::cout << "1 unmapped: " << unmapped1_count << "\n";
        
        std::cout << "2 ambiguous: " << ambiguous2_count << "\n";
        std::cout << "2 unmapped: " << unmapped2_count << "\n";
        
        std::cout << "Similar: " << similar_count << "\n";
    }
    
    void sim_diff()
    {
        std::cout << "Cigar of similar alignments differ by: \n";
        for(int i = 0; i < differences.size(); ++i){
           std::cout << differences[i] << ", ";
           if(i % 10 == 0 && i != 0)
               std::cout << "\n";
        }
    }
};


void load_alignments_of_read(BamFileIn & bamFileIn1,
                             BamFileIn & bamFileIn2,
                             auto & records1,
                             auto & records2,
                             BamAlignmentRecord & record1,
                             BamAlignmentRecord & record2,
                             bool start_sam,
                             bool verbose = false)
{
    bool empty_last_name = false;
    if(start_sam)
        empty_last_name = true;
    
    string last_name{};
    if(!start_sam){
        last_name = toCString(record1.qName);
        if(hasFlagFirst(record1)){
            records1.first.push_back(record1);
        }else{
            records1.second.push_back(record1);
        }
        if(hasFlagFirst(record2)){
            records2.first.push_back(record2);
        }else{
            records2.second.push_back(record2);
        }
    }else{
        start_sam = false;
    }
    
    while (!atEnd(bamFileIn1))
    {
        try
        {
            clear(record1);
            readRecord(record1, bamFileIn1);
            if(verbose)
                std::cout << toCString(record1.qName) << "\n";
            if (!empty_last_name){
                if (last_name.compare(toCString(record1.qName)) != 0)
                    break;
            }
            else
            {
                last_name = toCString(record1.qName);
                empty_last_name = false;
            }
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
        if(hasFlagFirst(record1)){
            records1.first.push_back(record1);
        }
        else
        {
            records1.second.push_back(record1);
        }
    }
    
    while (!atEnd(bamFileIn2))
    {
        try
        {
            clear(record2);
            readRecord(record2, bamFileIn2); 
            if (last_name.compare(toCString(record2.qName)) != 0)
                    break;
        }
        catch (IOError const & e)
        {
            std::cerr << "ERROR: could not copy record. " << e.what() << "\n";
        }
        if(hasFlagFirst(record2))
        {
            records2.first.push_back(record2);
        }
        else
        {
            records2.second.push_back(record2);
        }
    }
}

int main(int argc, char const * argv[])
{
    ArgumentParser parser("Search");

    addOption(parser, ArgParseOption("s1", "sam1", "Path to the sam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "sam1");

    addOption(parser, ArgParseOption("s2", "sam2", "Path to the sam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "sam2");
    
//     setRequired(parser, "reads1");

    addOption(parser, ArgParseOption("o", "output", "Path to output file", ArgParseArgument::OUTPUT_FILE, "OUT"));
//     setRequired(parser, "output");

    addOption(parser, ArgParseOption("bS", "batchSize", "Number of reads with are loaded at the same Time. Default = 100000", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString samPath1, samPath2, outputPath;
    int batchSize = 100000, barcodeLength;

    getOptionValue(samPath1, parser, "sam1");
    getOptionValue(samPath2, parser, "sam2");
    getOptionValue(outputPath, parser, "output");
    getOptionValue(batchSize, parser, "batchSize");
    bool verbose = isSet(parser, "verbose");

    //prepare loading sam
    BamFileIn bamFileIn1;
    BamHeader header1;
    if (!open(bamFileIn1, toCString(samPath1)))
    {
        std::cerr << "ERROR: could not open input file " << samPath1 << ".\n";
        return 1;
    }
    readHeader(header1, bamFileIn1);
//     std::cout << "check2\n";

    //prepare loading second sam
    BamFileIn bamFileIn2;
    BamHeader header2;
    if (!open(bamFileIn2, toCString(samPath2)))
    {
        std::cerr << "ERROR: could not open input file " << samPath2 << ".\n";
        return 1;
    }
    readHeader(header2, bamFileIn2);


    Stats stats;
    
    
    BamAlignmentRecord lrecord1;
    BamAlignmentRecord lrecord2;


    
    
    
    
/*
    //TODO if records > 2 than sort after  cordinate and seperate mates
    for(int i = 0; i < records1.first.size(); ++i){
        String<CigarElement<> > mycigar = records1.first[i].cigar;
        std::cout << toCString(records1.first[i].qName) << "\n";
        if (length(mycigar) == 0)
            std::cout << "No alignment found\n";
//         std::cout << "Length of Cigar: " << length(records1[i].cigar) << "\n";
        for(int j = 0; j < length(mycigar); ++j)
            std::cout << mycigar[j].count << mycigar[j].operation << " p";
         std::cout << "\n";
    }*/

    bool start_sam = true; 
    for(int r = 0; r < 16; ++r)
    {
        ++stats.reads_count;
    pair<std::vector<BamAlignmentRecord>, std::vector<BamAlignmentRecord> > records1;
    pair<std::vector<BamAlignmentRecord>, std::vector<BamAlignmentRecord> > records2;
    
    
    load_alignments_of_read(bamFileIn1, bamFileIn2, records1, records2, lrecord1, lrecord2, start_sam);
    if(start_sam)
        start_sam = false;
    
    std::cout << "Read: " << r << "\n";
    std::cout << "1: " << toCString(records1.first[0].qName) << "\n2: " <<  toCString(records2.first[0].qName) << "\n";
    std::cout << "1m: " << toCString(records1.second[0].qName) << "\n2m: " <<  toCString(records2.second[0].qName) << "\n";
    
//     std::cout << toCString(records1.first[0].qName) << "\n";
    std::cout << "Positions: " << records1.first[0].beginPos << "\t" << records2.first[0].beginPos << "\n";
    std::cout << "Positions: " << records1.second[0].beginPos << "\t" << records2.second[0].beginPos << "\n";

//     std::cout << "Cigar length: " << length(records1.first[0].cigar) << "\t" << length(records2.first[0].cigar) << "\n";
//     std::cout << "Cigar length: " << length(records1.second[0].cigar) << "\t" << length(records2.second[0].cigar) << "\n";
    
    /*
    string tmp = toCString(records1.first[0].qName);
    if(tmp.compare(toCString(records1.second[0].qName)) != 0){
        std::cout << "readnames differ\n";
        exit(0);
    }*/
    
    for(int i = 0; i < 2; ++i)
    {
        
        BamAlignmentRecord record1;
        BamAlignmentRecord record2;
        
        bool similar = false;
        bool unmapped = false, unmapped1 = false, unmapped2 = false;
        bool ambiguously = false, ambiguously1 = false, ambiguously2 = false;
        if(i == 0){
            record1 = records1.first[0];
            record2 = records2.first[0];
            
            ambiguously1 = records1.first.size() > 1;//hasFlagMultiple(record1); //flag is also set if alignments is not saved in file?
            ambiguously2 = records2.first.size() > 1;//hasFlagMultiple(record2);
        }
        else
        {
            record1 = records1.second[0];
            record2 = records2.second[0];
//             std::cout << toCString(record1.qName) << "\t" << hasFlagUnmapped(record1); << "--adfads\n";
            
            ambiguously1 = records1.second.size() > 1;//hasFlagMultiple(record1); //flag is also set if alignments is not saved in file?
            ambiguously2 = records2.second.size() > 1;//hasFlagMultiple(record2);
        }
        

        
        unmapped1 = length(record1.cigar) == 0;  //hasFlagUnmapped(record1);
        unmapped2 = length(record2.cigar) == 0;  //hasFlagUnmapped(record2);
        
        stats.ambiguous1_count += ambiguously1;
        stats.ambiguous2_count += ambiguously2;
        if(ambiguously1 && ambiguously2){
            ++stats.ambiguous_count;
        }
        
        if(unmapped1 && unmapped2)
            ++stats.unmapped_count;
        stats.unmapped1_count += unmapped1;
        stats.unmapped2_count += unmapped2;
        
        
        //check more than POS
        std::cout << "Status " << unmapped1 << unmapped2 << ambiguously1 << ambiguously2 << "\n";
        
        if(!unmapped1 && !unmapped2 && !ambiguously1 && !ambiguously2){
            std::cout << "Check Similar\n";
            std::cout << "1: " << record1.beginPos << "\t2:" << record2.beginPos << "\n";
            similar = (record1.beginPos < (record2.beginPos + 10) && record1.beginPos + 10 > (record2.beginPos));
        }
        
        stats.similar_count += similar;
        
        
        //compareCigarStrings only with uniquly mapping reads
        if(similar)
        {
            std::cout << "Check Cigar\n";
            String<CigarElement<> > cigar1 = record1.cigar;
            String<CigarElement<> > cigar2 = record2.cigar;
            int16_t pos1 = 0, pos2 = 0, opos1 = 0, opos2 = 0;
            
            uint32_t dif = 0;
            while(true)
            {
                if(pos1 >= length(cigar1) && pos2 >= length(cigar2))
                {
                    break;
                }else if(pos1 >= length(cigar1) || pos2 >= length(cigar2))
                {
                    ++dif;
                }
                else if (cigar1[pos1].operation != cigar2[pos2].operation)
                {
                    ++dif;
                }
                //TODO jump muliple position in a single operation
                ++opos1;
                ++opos2;
                if(opos1 >= cigar1[pos1].count){
                    ++pos1;
                    opos1 = 0;
                }
                
                if(opos2 >= cigar1[pos2].count){
                    ++pos2;
                    opos2 = 0;
                }
            }
            std::cout << "Diff: " << dif << "\n";
            stats.differences.push_back(dif);
        }
    }
//         std::cout << "loaded1: " << (records1.first.size() + records1.second.size()) << " alignments\n";
//         std::cout << "loaded2: " << (records2.first.size() + records2.second.size()) << " alignments\n";
    }
    std::cout << "\n";
    stats.print();
    stats.sim_diff();
    
    std::cout << "Finished!\n";

    return 0;
}



