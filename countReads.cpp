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
#include <iostream>
#include <fstream>

using namespace seqan;
using namespace std;



int main(int argc, char const * argv[])
{
    ArgumentParser parser("Count Reads");

    addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "bam");
    
//     addOption(parser, ArgParseOption("g", "gtf", "Path to the gtf file", ArgParseArgument::INPUT_FILE, "IN"));
//     setRequired(parser, "gtf");
    
    
//     addOption(parser, ArgParseOption("su", "suffix", "Suffix concatenated to ouput files (beside .bam)",
//                                      ArgParseOption::STRING));
    
//     addOption(parser, ArgParseOption("o", "output", "Path to the output prefix", ArgParseArgument::INPUT_FILE, "IN"));
    
//     addOption(parser, ArgParseOption("s", "step", "Number of exons in single bam", ArgParseArgument::INTEGER, "INT"));
    
//     addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
    
/*
    addOption(parser, ArgParseOption("p", "first", "First p reads", ArgParseArgument::INTEGER, "INT"));
    hideOption(getOption(parser, "first"));

    addOption(parser, ArgParseOption("mt", "threshold", "Number of matches required to define start or end of the cDNA part", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));*/

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString bamPath, gtfPath, outputPathPrefix = "reads", suffix = "";


    getOptionValue(bamPath, parser, "bam");


    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamPath)))
    {
        std::cerr << "ERROR: Could not open " << bamPath << std::endl;
        return 1;
    }
    
    uint64_t reads = 0;
    uint64_t unmapped = 0;
    uint64_t multi = 0;

    BamHeader header;
    std::vector<std::vector<BamAlignmentRecord > > recordtable;
    try
    {
        // Copy header.
        readHeader(header, bamFileIn);
        // Copy records.
        BamAlignmentRecord record;
//         recordtable.resize(table.size());
        

    
        string lastrecordName = "";
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            if(lastrecordName.compare(toCString(record.qName)) == 0){
                ++multi;
            }else{
                lastrecordName = toCString(record.qName);
            }
            ++reads;
            //use map to jump to correct chromosom //use start pos and length
            if(hasFlagUnmapped(record))
                ++unmapped;
        }
    }
    catch (Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
    
    std::cout << "Reads: " << reads << "\n";
    std::cout << "Unmmapped: " << unmapped << "\n";
    std::cout << "Multi: " << multi << "\n";

    return 0;
}




/*
    
    std::vector<ofstream> filestreams(table.size());
    
    std::vector<BamFileOut> bamfileouts; 
    
    string prefix = "read";
    for(int i = 0; i < table.size(); ++i){
//         ofstream mybamstream;
//         filestreams.push_back(mybamstream);
        string filename = "read.bam";//prefix + to_string(i) + ".bam";
//         mybamstream.open(filename);
        filestreams[i].open(filename);
//         BamFileOut bamfileout{context(bamFileIn), filestreams[i], Bam()};
        bamfileouts.push_back(BamFileOut{context(bamFileIn), filestreams[i], Bam()});
    }*/
