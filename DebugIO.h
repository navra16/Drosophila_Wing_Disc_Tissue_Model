#pragma once
#include <string>
struct LambdaField;                 // forward declaration

void dumpLambdaCSV(const LambdaField& field,
                   const std::string& fileName,
                   std::size_t maxToWrite = 0);
