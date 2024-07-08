
#ifndef CSV_HANDLER_H
#define CSV_HANDLER_H

#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream

class CSVHandler {
public:
    static std::vector<std::pair<std::string, std::vector<int>>> read_csv(const std::string& filename);
    static void write_csv(const std::string& filename, const std::vector<std::pair<std::string, std::vector<int>>>& dataset, bool append = false);
    static void write_matrix_csv(const std::string& filename, const std::vector<std::vector<std::vector<int>>>& data);
};
#endif 