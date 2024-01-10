#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <map>
#include <thread>
#include <boost/icl/interval_set.hpp>
#include <boost/icl/interval.hpp>
#include <omp.h>
#include <vector>
#include <unordered_set>

struct Tss
{
    std::string chromosome;
    std::string start;
    std::string end;
    std::string id;
    std::string strand;
    int counts = 1;
    double transcriptLengths = 0.0;
};

std::unordered_set<std::string> mappingReads(std::string const bedFilePath, std::unordered_map<std::string, std::unordered_map<std::string, Tss>> &chromosomeMap);
void analyzeChromosome(const std::unordered_map<std::string, Tss> &chromosomeData, const int tssRegion, std::vector<Tss> &outputVector, const int minTssCounts);
void analyzeStrandVectors(const std::vector<Tss> &tssStruct, const int tssRegion, std::vector<Tss> &results, const int minTssCounts);

// sort the vectors by counts, chromosome, and start position
bool tssComparator(const Tss &a, const Tss &b)
{
    // Sort by counts in descending order
    return a.counts > b.counts;
}

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " <bed file> <TSS region> <minimum counts per TSS> <yes/no for bed output for conversion to bedgraph" << std::endl;
        return 1;
    }

    std::string bedFilePath, bedOutput;
    try
    {
        bedFilePath = argv[1];
        bedOutput = argv[4];
    }
    catch (const std::invalid_argument &e)
    {
        std::cerr << "Invalid argument. Please provide a valid path to the bed file." << std::endl;
        return 1;
    }

    if (bedOutput != "yes" && bedOutput != "no")
    {
        std::cerr << "Invalid argument. Please provide yes or no for bed output." << std::endl;
        return 1;
    }

    int tssRegion, minTssCounts;
    try
    {
        tssRegion = std::stoi(argv[2]);
        minTssCounts = std::stoi(argv[3]);
    }
    catch (const std::invalid_argument &e)
    {
        std::cerr << "Invalid argument. Please provide valid integers for tssRegion and minTssCounts." << std::endl;
        return 1;
    }

    // Create a map  <chromosome,<key,value>>
    std::unordered_map<std::string, std::unordered_map<std::string, Tss>> chromosomeMap;
    std::unordered_set<std::string> uniqueChromosomes = mappingReads(bedFilePath, chromosomeMap);
    std::vector<std::string> chromosomeVector(uniqueChromosomes.begin(), uniqueChromosomes.end());

    // Set the number of threads to the number of available cores
    omp_set_num_threads(std::thread::hardware_concurrency());
    std::vector<std::vector<Tss>> allOutputs(uniqueChromosomes.size());

    // Parallelize the processing of each chromosome
#pragma omp parallel for
    for (size_t i = 0; i < chromosomeVector.size(); ++i)
    {
        const std::string &chromosomeKey = chromosomeVector[i];
        if (chromosomeMap.find(chromosomeKey) != chromosomeMap.end())
        {
            analyzeChromosome(chromosomeMap[chromosomeKey], tssRegion, allOutputs[i], minTssCounts);
        }
    }

    // Merge all vectors into one
    std::vector<Tss> finalOutputVector;
    for (const auto &vector : allOutputs)
    {
        finalOutputVector.insert(finalOutputVector.end(), vector.begin(), vector.end());
    }

    // Write finalOutput to a BED file
    std::ofstream outFile("tss-finder-output.bed");
    if (!outFile.is_open())
    {
        std::cerr << "Failed to open the output file for writing.\n";
        return 1;
    }
    outFile << "chromosome\tstart\tend\tcounts\tavg. transcript lengths\tstrand\n";
    for (const auto &tss : finalOutputVector)
    {
        int start = std::stoi(tss.start);
        outFile << tss.chromosome << "\t"
                << tss.start << "\t"
                << start + 1 << "\t"
                << tss.counts << "\t"
                << tss.transcriptLengths << "\t"
                << tss.strand << "\n";
    }
    // write to bed file for conversion to bedgraph
    if (bedOutput == "yes")
    {
        std::ofstream bedFile("tss-finder-bed-for-bg.bed");
        if (!bedFile.is_open())
        {
            std::cerr << "Failed to open the output file for writing.\n";
            return 1;
        }
        for (const auto &tss : finalOutputVector)
        {
            for (size_t i = 0; i < tss.counts; ++i)
            {
                int start = std::stoi(tss.start);
                bedFile << tss.chromosome << "\t"
                        << tss.start << "\t"
                        << start + 1 << "\t"
                        << tss.id << "\t"
                        << "."
                        << "\t"
                        << tss.strand << "\n";
            }
        }
        bedFile.close();
    }
    outFile.close();

    return 0;
}

std::unordered_set<std::string> mappingReads(std::string const bedFilePath, std::unordered_map<std::string, std::unordered_map<std::string, Tss>> &chromosomeMap)
{
    // set with chromosomes
    std::unordered_set<std::string> setChromosomes;
    // open bed file
    std::ifstream bed_file(bedFilePath);
    if (!bed_file.is_open())
    {
        std::cerr << "Failed to open the bed file: " << bedFilePath << std::endl;
        std::exit(1);
    }

    // loop through the bed file
    std::string chromosome, id, start, end, quality, strand;
    while (bed_file >> chromosome >> start >> end >> id >> quality >> strand)
    {
        // add the chromosome to the set
        setChromosomes.insert(chromosome);
        // <chromosome,<key,value>>
        std::string startSite;
        if (strand == "+")
        {
            startSite = start;
        }
        else
        {
            startSite = end;
        }

        std::string key = chromosome + startSite + strand;
        Tss tss;
        int startRead = std::stoi(start);
        int endRead = std::stoi(end);
        double length = (endRead - startRead) + 1.0;

        // if chromsome is not in the map, add it and add the tss_map to it
        if (chromosomeMap.find(chromosome) == chromosomeMap.end())
        {
            std::unordered_map<std::string, Tss> tss_map;
            chromosomeMap[chromosome] = tss_map;
            chromosomeMap[chromosome][key] = tss;
            chromosomeMap[chromosome][key].chromosome = chromosome;
            chromosomeMap[chromosome][key].start = startSite;
            chromosomeMap[chromosome][key].end = startSite;
            chromosomeMap[chromosome][key].id = id;
            chromosomeMap[chromosome][key].strand = strand;
            chromosomeMap[chromosome][key].transcriptLengths = length;
        }
        else if (chromosomeMap[chromosome].find(key) == chromosomeMap[chromosome].end())
        {
            chromosomeMap[chromosome][key] = tss;
            chromosomeMap[chromosome][key].chromosome = chromosome;
            chromosomeMap[chromosome][key].start = startSite;
            chromosomeMap[chromosome][key].end = startSite;
            chromosomeMap[chromosome][key].id = id;
            chromosomeMap[chromosome][key].strand = strand;
            chromosomeMap[chromosome][key].transcriptLengths = length;
        }
        else
        {
            chromosomeMap[chromosome][key].counts++;
            chromosomeMap[chromosome][key].transcriptLengths += length;
        }
    }
    // Calculate the average transcript length
    for (auto &chromosome : chromosomeMap)
    {
        for (auto &tss : chromosome.second)
        {
            tss.second.transcriptLengths = tss.second.transcriptLengths / tss.second.counts;
        }
    }
    // close bed file
    bed_file.close();
    return setChromosomes;
}

void analyzeChromosome(const std::unordered_map<std::string, Tss> &chromosomeData, const int tssRegion, std::vector<Tss> &outputVector, const int minTssCounts)
{
    // size of the map
    int size = chromosomeData.size();

    // create two maps, one for positive strand and one for negative strand
    std::vector<Tss> vectorPlusStrand;
    std::vector<Tss> vectorMinusStrand;

    // loop through the chromosome data
    for (const auto &entry : chromosomeData)
    {
        const std::string &key = entry.first;
        const Tss &tss = entry.second;

        // if the strand is positive, add it to the positive strand vector
        if (tss.strand == "+")
        {
            vectorPlusStrand.push_back(tss);
        }
        // if the strand is negative, add it to the negative strand vector
        else
        {
            vectorMinusStrand.push_back(tss);
        }
    }

    std::sort(vectorPlusStrand.begin(), vectorPlusStrand.end(), tssComparator);
    std::sort(vectorMinusStrand.begin(), vectorMinusStrand.end(), tssComparator);

    // Parallel processing of positive and negative strand vectors
    std::vector<Tss> outputPositiveStrand, outputNegativeStrand;
#pragma omp parallel sections
    {
#pragma omp section
        {
            analyzeStrandVectors(vectorPlusStrand, tssRegion, outputPositiveStrand, minTssCounts);
        }
#pragma omp section
        {
            analyzeStrandVectors(vectorMinusStrand, tssRegion, outputNegativeStrand, minTssCounts);
        }
    }

    // merge the two vectors
    outputVector.insert(outputVector.end(), outputPositiveStrand.begin(), outputPositiveStrand.end());
    outputVector.insert(outputVector.end(), outputNegativeStrand.begin(), outputNegativeStrand.end());
}

void analyzeStrandVectors(const std::vector<Tss> &tssStruct, const int tssRegion, std::vector<Tss> &results, const int minTssCounts)
{
    // create a map of intervals
    typedef boost::icl::discrete_interval<int> GenomicInterval;
    boost::icl::interval_set<int> tree;

    // Loop through the vector
    for (const auto &tss : tssStruct)
    {
        int start = std::stoi(tss.start);
        // Create an interval to check
        GenomicInterval interval = GenomicInterval::closed(start, start);

        // Check if the interval is not in the tree
        if ((!boost::icl::intersects(tree, interval)) && (tss.counts >= minTssCounts))
        {
            // Find the upper and lower intervals
            auto upperIt = tree.upper_bound(interval);
            // Check the it is not empty
            int upperBound, lowerBound;
            if (upperIt == tree.end())
            {
                upperBound = start + tssRegion; // edge case: the upper interval does not exist
            }
            else
            {
                int distanceToUpper = upperIt->lower() - start;                                         // Calcualte the distance to the first bound of the upper interval
                upperBound = distanceToUpper > tssRegion ? start + tssRegion : start + distanceToUpper; // Check if the distance is greater than the tssRegion
            }
            // Check if the lower interval exists
            auto lowerIt = upperIt;
            if (upperIt != tree.begin())
            {
                --lowerIt;                                                                              // Move to the previous interval
                int distanceToLower = start - lowerIt->upper();                                         // Calcualte the distance to the second bound of the upper interval
                lowerBound = distanceToLower > tssRegion ? start - tssRegion : start - distanceToLower; // Check if the distance is greater than the tssRegion
            }
            else
            {
                lowerBound = start - tssRegion; // edge case: the lower interval does not exist
                if (lowerBound < 1)
                {
                    lowerBound = 1; // edge case: the lower interval is 0 or negative
                }
            }

            // Create an interval
            GenomicInterval intervalToAdd = GenomicInterval::closed(lowerBound, upperBound);
            tree.insert(intervalToAdd);
            results.push_back(tss);
        }
    }
}