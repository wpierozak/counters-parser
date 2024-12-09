#include<fstream>
#include<iostream>
#include<cstring>
#include<list>
#include<algorithm>
#include<cmath>
#include"variadic_table/include/VariadicTable.h"

constexpr std::string_view RatesLine = "RATES,";
constexpr std::string_view EmptyRatesLine = "RATES,-";

std::string TcmCounters[] = {"TRG 5", "TRG 4", "TRG 2", "TRG 1", "TGR 3", "BCK 0", "BCK 1",
"BCK 2", "BCK 3", "BCK 4", "BCK 5", "BCK 6", "BCK 7", "BCK 8", "BCK 9"};


template<size_t RatesNumber>
struct Row{
    double rates[RatesNumber];
};

template<size_t RatesNumber, size_t BufferSize>
std::list<Row<RatesNumber>> rateParser(const char* filename){
    std::list<Row<RatesNumber>> results;
    std::ifstream file(filename);
    if(file.is_open()==false){
        std::cerr<< "Failed to open " << filename << std::endl;
        abort();
    }
    char buffer[BufferSize];
    file.rdbuf()->pubsetbuf(buffer, BufferSize);
    std::string line; 
    line.reserve(2048);
    size_t lineN = 0;
    while(std::getline(file, line)){
        
        if(line.length() < EmptyRatesLine.length()){
            continue;
        }
        //std::cout << "Line: " << lineN++ << ": " << line << std::endl;
        if( (std::strncmp(line.c_str(), RatesLine.data(), RatesLine.length()) == 0) &&
            (std::strncmp(line.c_str(), EmptyRatesLine.data(), EmptyRatesLine.length()) != 0)){
            
            Row<RatesNumber> newRow;
            size_t pos = -1;
            size_t cnt = 0;
            while((pos = line.find(',', pos+1)) != std::string::npos){
                size_t end = line.find(',', pos+1);
                newRow.rates[cnt++] = std::stod(line.substr(pos+1, end));
            }
            if(cnt < RatesNumber){
                std::cerr << "Missing rates: received" << cnt << ", expected " << RatesNumber << std::endl;
            }
            results.emplace_back(std::move(newRow));
        }
    }
    return results;
}

template<size_t RatesNumber>
void parseRates(const char* fileName){
    std::list<Row<RatesNumber>> rates = rateParser<RatesNumber,(1u<<12)>(fileName);
    std::cout << "Parsed file" << std::endl;
    double mean[RatesNumber] = {0};
    double stddev[RatesNumber] = {0};
    double coeffVar[RatesNumber] = {0};
    double moreThanTwoStddev[RatesNumber]={0};
    double moreThanThreeStddev[RatesNumber]={0};

    std::for_each(rates.begin(), rates.end(), 
                [&mean](const Row<RatesNumber>& row){
                        for(int i= 0; i<RatesNumber;i++){
                            mean[i]+=row.rates[i];
                        }
                    }
                );
    for(int i = 0; i <RatesNumber; i++){
        mean[i] /= static_cast<double>(rates.size());
    }
    std::for_each(rates.begin(), rates.end(),
                [&stddev,&mean](const Row<RatesNumber>&row){
                        for(int i = 0; i <RatesNumber; i++){
                            stddev[i] += (row.rates[i]-mean[i])*(row.rates[i]-mean[i]);
                        }
                    }
                );
    
    for(int i = 0; i <RatesNumber; i++){
        stddev[i] =  sqrt(stddev[i]/static_cast<double>(rates.size()));
        coeffVar[i] = stddev[i]/mean[i];
    }

    std::for_each(rates.begin(), rates.end(),
                [mean,&stddev,&moreThanThreeStddev,&moreThanTwoStddev](const Row<RatesNumber>&row){
                        for(int i = 0; i <RatesNumber; i++){
                            if(abs(row.rates[i] - mean[i]) > 2*stddev[i]){
                                moreThanTwoStddev[i]++;
                            }
                            if(abs(row.rates[i] - mean[i]) > 3*stddev[i]){
                                moreThanThreeStddev[i]++;
                            }
                        }
                    }
                );

    VariadicTable<std::string,size_t,double,double,double,size_t,size_t> vt({"Rate", "Samples", "Mean", "Stddev", "Coefficient of variation", ">2*stddev", ">3*stddev"},15);
    
    for(int i = 0; i<RatesNumber; i++){
        if(RatesNumber == 15){
            vt.addRow(TcmCounters[i],rates.size(), mean[i],stddev[i],coeffVar[i],moreThanTwoStddev[i],moreThanThreeStddev[i]);
        }
        else{

            std::string rateName = "CH"+std::to_string(i/2+1);
            rateName += (i%2 == 0) ? "_CFD" : "_TRG";
            vt.addRow(rateName,rates.size(), mean[i],stddev[i],coeffVar[i],moreThanTwoStddev[i],moreThanThreeStddev[i]);
        }
    }
    vt.print(std::cout);
}

int main(int argc, const char** argv)
{
    if(argc != 3){
        abort();
    }
    if(strcmp(argv[1],"TCM")==0){
        parseRates<15>(argv[2]);
    }
    if(strcmp(argv[1],"PM")==0){
        parseRates<24>(argv[2]);
    }
    return 0;
}