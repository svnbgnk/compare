#include <seqan/sequence.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>

#include <iostream>
#include <seqan/align.h>
#include <seqan/bam_io.h>
#include <iostream>
#include <fstream>
#include <seqan/find.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <map>
#include <stdlib.h> 
#include <stdio.h>


using namespace seqan;
using namespace std;





int main(int argc, char const * argv[])
{
    string filepath{"short_example.txt"};
    ifstream inputFile(filepath);
    string line;
    vector<vector<string> > table;
    vector<double> probs;
    if(inputFile.is_open()){
        bool head = true;
        while (getline(inputFile, line, '\n'))
        {
            if(head){
                head = false;
                continue;
            }
                
            vector<string> row;
            std::istringstream sline(line);
            string element;
            int k = 0;
            while (getline(sline, element, '\t')){
                if(k == 0 || k == 1)
                    row.push_back(element);
                if(k == 1){
                    double tmp = atof(element.c_str());
                    probs.push_back(tmp);
                }
                ++k;
            }
            table.push_back(row);
        }
        
        std::cout << "Finished reading\n";
        for (unsigned i = 0; i < table.size(); i++){
            for(int j = 0; j < table[i].size(); ++j){
                std::cout << table[i][j] << "\t";
            }
            std::cout << "\n";
        }
        
        std::cout << "Probs\n";
        for(int i = 0; i < probs.size(); ++i){
            std::cout << probs[i] << "\t";
        }
        std::cout << "\n";
        
        inputFile.close();
    }
    else
    {
        std::cout << "Could not open file: " << filepath << "\n";
    }

    std::cout << "Finished!\n";

    return 0;
}



