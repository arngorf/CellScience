#ifndef DEBUG_HPP
#define DEBUG_HPP

#include <cstdio>
#include <fstream>
#include <string>

#include "UtilityFunctions.hpp"

class Debug
{
public:
    Debug();

    ~Debug();

    static void Info(const std::string msg);

    static void Warning(const std::string msg);

    static void Error(const std::string msg);

    static void Tick();

    static void Tock(std::string identifier);

private:
    static std::ofstream log;

    static void logMsg(const std::string msg);

    static unsigned long int elapsedTime;
    static unsigned long int previousTimerLine;

    static unsigned int debugLevel;
};

#endif // DEBUG_HPP