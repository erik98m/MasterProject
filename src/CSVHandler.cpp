#include "CSVHandler.h"
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

// code inspiration: https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/

std::vector<std::pair<std::string, std::vector<int>>> CSVHandler::read_csv(const std::string& filename)
{
    std::vector<std::pair<std::string, std::vector<int>>> result;
    std::ifstream myFile(filename);

    if (!myFile.is_open()) throw std::runtime_error("Could not open file");

    std::string line, colname;
    double val;

    if (myFile.good())
    {
        std::getline(myFile, line);
        std::stringstream ss(line);

        while (std::getline(ss, colname, ','))
        {
            result.push_back({ colname, std::vector<int>{} });
        }
    }

    while (std::getline(myFile, line))
    {
        std::stringstream ss(line);
        int colIdx = 0;

        while (ss >> val)
        {
            result.at(colIdx).second.push_back(val);
            if (ss.peek() == ',') ss.ignore();
            colIdx++;
        }
    }

    myFile.close();
    return result;
}

void CSVHandler::write_csv(const std::string& filename, const std::vector<std::pair<std::string, std::vector<int>>>& dataset, bool append)
{
    std::ofstream myFile;
    if (append) {
        myFile.open(filename, std::ios::app);
    } else {
        myFile.open(filename);
    }

    if (!append) {
        for (int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).first;
            if (j != dataset.size() - 1) myFile << ",";
        }
        myFile << "\n";
    }

    for (int i = 0; i < dataset.at(0).second.size(); ++i)
    {
        for (int j = 0; j < dataset.size(); ++j)
        {
            myFile << dataset.at(j).second.at(i);
            if (j != dataset.size() - 1) myFile << ",";
        }
        myFile << "\n";
    }

    myFile.close();
}

void CSVHandler::write_matrix_csv(const std::string& filename, const std::vector<std::vector<std::vector<int>>>& data)
{
    std::ofstream myFile;
    myFile.open(filename);

    for (size_t iter = 0; iter < data.size(); ++iter) {
        myFile << "Iteration " << iter << "\n";
        const auto& matrix = data[iter];
        for (const auto& row : matrix) {
            for (size_t j = 0; j < row.size(); ++j) {
                myFile << row[j];
                if (j != row.size() - 1) myFile << ",";
            }
            myFile << "\n";
        }
        myFile << "\n"; // Add a new line between iterations
    }

    myFile.close();
}
