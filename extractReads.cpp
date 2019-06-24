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
    ArgumentParser parser("Extract Reads");

    addOption(parser, ArgParseOption("b", "bam", "Path to the bam file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "bam");
    
    addOption(parser, ArgParseOption("g", "gtf", "Path to the gtf file", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "gtf");
    
    addOption(parser, ArgParseOption("o", "output", "Path to the output prefix", ArgParseArgument::INPUT_FILE, "IN"));
    setRequired(parser, "gtf");
    
    addOption(parser, ArgParseOption("t", "threads", "Number of threads", ArgParseArgument::INTEGER, "INT"));
/*
    addOption(parser, ArgParseOption("p", "first", "First p reads", ArgParseArgument::INTEGER, "INT"));
    hideOption(getOption(parser, "first"));

    addOption(parser, ArgParseOption("mt", "threshold", "Number of matches required to define start or end of the cDNA part", ArgParseArgument::INTEGER, "INT"));

    addOption(parser, ArgParseOption("v", "verbose", ""));*/

    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    CharString bamPath, gtfPath, outputPathPrefix = "reads";
    int batchSize1 = 100000;

    int barcode_umi_length = 30;
    int threshold = 5;
    int first = 9999999;
    int threads = 1;

    getOptionValue(bamPath, parser, "bam");
    getOptionValue(gtfPath, parser, "gtf");
    getOptionValue(outputPathPrefix, parser, "output");
    getOptionValue(threads, parser, "threads");
//     getOptionValue(barcodeLength, parser, "barcodeL");
//     getOptionValue(first, parser, "first");
//     getOptionValue(threshold, parser, "threshold");
//     bool verbose = isSet(parser, "verbose");


    //Read gtf file
    ifstream inputFile(toCString(gtfPath));
    string line;
    vector<std::tuple<CharString, uint32_t, uint32_t> > table;
    if(inputFile.is_open()){
//         bool head = true;
        while (getline(inputFile, line, '\n'))
        {
            if(line.length() < 10)
                continue;
//             if(head){
//                 head = false;
//                 continue;
//             }
            std::tuple<CharString, uint32_t, uint32_t> row;
            std::istringstream sline(line);
            string element;
            std::vector<string> tmprow;
            int k = 0;
            while (getline(sline, element, '\t')){
//                 std::cout << "k: " << k << "\t" << element << "\n";
                if(k == 0 || k == 3 || k == 4)
                    tmprow.push_back(element);
                ++k;
            }
            CharString refid = tmprow[0];
            std::cout << tmprow[0] << "\t" << tmprow[1] << "\t" << tmprow[2] << "\n";
            int s = std::stoi(tmprow[1]);
            int e = std::stoi(tmprow[2]);
            table.push_back(std::make_tuple(refid, s, e));
        }
        std::cout << "Finished reading\n";/*
        for (unsigned i = 0; i < table.size(); i++){
                std::cout << std::get<0>(table[i]) << "\t" << std::get<1>(table[i]) << "\t" << std::get<2>(table[i]) << "\n";
        }
        std::cout << "\n";*/
        inputFile.close();
    }
    else
    {
        std::cout << "Could not open file: " << gtfPath << "\n";
    }

    // Open input file, BamFileIn can read SAM and BAM files.
    BamFileIn bamFileIn;
    if (!open(bamFileIn, toCString(bamPath)))
    {
        std::cerr << "ERROR: Could not open " << bamPath << std::endl;
        return 1;
    }

    BamHeader header;
    std::vector<std::vector<BamAlignmentRecord > > recordtable;
    try
    {
        // Copy header.
        readHeader(header, bamFileIn);
        // Copy records.
        BamAlignmentRecord record;
        recordtable.resize(table.size());
        
        string lastContig;
        int st = 0;
        int end = table.size() - 1;
        bool newChrom = true;
        bool stNew = true;
        bool endNew = false;
        while (!atEnd(bamFileIn))
        {
            readRecord(record, bamFileIn);
            //loop
            //use map to jump to correct chromosom //use start pos and length
            string recordContig = toCString(getContigName(record, bamFileIn));
            uint32_t recordBegin = record.beginPos;
            uint32_t recordEnd = recordBegin + length(record.seq);
            
            //determine search range
            if(lastContig.compare(recordContig) != 0){
                stNew = true;
                for(int i = 0; i < table.size(); ++i){
                    string rowContig = toCString(std::get<0>(table[i]));
                    std::cout << "Check rowContig: " << rowContig << "\t" << recordContig << "\n";
                    if(recordContig.compare(rowContig) == 0 && stNew){
                        stNew = false;
                        st = i;
                    }
                    if(!stNew && recordContig.compare(rowContig) != 0)
                    {
                        end = i;
                        break;
                    }
                    // in case last element matches
                    if(recordContig.compare(rowContig) == 0)
                        end = table.size();
                }
                
                lastContig = recordContig;
                std::cout << "Start: " << st << "\tEnd: " << end << "\n";
            }
            
            #pragma omp parallel for num_threads(threads) schedule(static)
            for(int i = st; i < end; ++i){
                //check row
                string rowContig = toCString(std::get<0>(table[i]));
                uint32_t rowBegin = std::get<1>(table[i]);
                uint32_t rowEnd = std::get<2>(table[i]);
                //std::cout << "recordBegin: " << recordBegin << "\t" << recordEnd << "\trow: " << rowBegin << "\t" << rowEnd << "\n";
                if((recordBegin > rowBegin && recordBegin < rowEnd) || (recordEnd > rowBegin && recordEnd < rowEnd))
                {
                    recordtable[i].push_back(record);
                }
            }
            
//             writeRecord(bamFileOut, record);
        }
        
        string prefix = toCString(outputPathPrefix);
        #pragma omp parallel for schedule(dynamic) num_threads(threads)
        for(int b = 0; b < recordtable.size(); ++b)
        {
            ofstream mybamstream;
            string bamName = prefix + to_string(b) + ".bam";
            mybamstream.open(bamName);
            // Open output file, BamFileOut accepts also an ostream and a format tag.
            BamFileOut bamFileOut(context(bamFileIn), mybamstream, Bam());
            writeHeader(bamFileOut, header);
            for(int i = 0; i < recordtable[b].size(); ++i){
                writeRecord(bamFileOut, recordtable[b][i]);
            }
            close(bamFileOut);
        }
    }
    
    
    catch (Exception const & e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }

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