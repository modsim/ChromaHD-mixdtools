#include "mixd.hpp"

void printUsage(char *binaryName)
{
    std::cout << "MIXD tool to extract a timestep from data.all." << std::endl;
    std::cout << "Might be better to just use `split -n <ts>/<nts> data.all > data.extracted`" << std::endl;
    std::cout << "Might be better to just use scripts/mslice" << std::endl;
    std::cout << "Usage: " << binaryName << " -f <data.all> -m <minf> -n <ndf> -t <timestep> -o <data.in>" << std::endl;
    std::cout << "<timestep> count starts with 0 as data.all stores initial state." << std::endl;
}

//Copied from: https://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
class InputParser{
    public:
        InputParser (int &argc, char **argv){
            for (int i=1; i < argc; ++i)
                this->tokens.push_back(std::string(argv[i]));
        }
        /// @author iain
        const std::string& getCmdOption(const std::string &option) const{
            std::vector<std::string>::const_iterator itr;
            itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
            if (itr != this->tokens.end() && ++itr != this->tokens.end()){
                return *itr;
            }
            static const std::string empty_string("");
            return empty_string;
        }
        /// @author iain
        bool cmdOptionExists(const std::string &option) const{
            return std::find(this->tokens.begin(), this->tokens.end(), option)
                   != this->tokens.end();
        }
    private:
        std::vector <std::string> tokens;
};


int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;

    /* if(argc != 3) */
    /* { */
    /*     printUsage(argv[0]); */
    /*     return 1; */
    /* } */

    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help")) {
        printUsage(argv[0]);
        return 1;
    }

    const std::string &data_filename          = input.getCmdOption("-f");
    const std::string &minf_filename          = input.getCmdOption("-m");
    const std::string &out_filename           = input.getCmdOption("-o");
    const std::string &ndf_str                = input.getCmdOption("-n");
    const std::string &ts_str                 = input.getCmdOption("-t");

    int ts = std::stoi(ts_str);
    int ndf = std::stoi(ndf_str);

    try
    {

        long ne, nn;
        readminf(minf_filename, &nn, &ne);

        MixdFile<double> data(data_filename, nn, ndf);
        data.read(0, false, (ts)*nn);

        data.fname(out_filename);
        data.write();

    }
    catch(mixd::MixdException e)
    {
        std::cout << e.msg() << std::endl;
    }

    return 0;
}
