#include "Debug.hpp"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>

#include "timer.h"

std::ofstream Debug::log;

unsigned long int Debug::elapsedTime;
unsigned long int Debug::previousTimerLine;

unsigned int Debug::debugLevel = 3;

Debug::Debug()
{

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::ostringstream oss;
    oss << std::put_time(&tm, "%d%m%Y-%H%M%S");
    auto str = oss.str();

    std::cout << "Creating logfile for " << str << std::endl;


    std::string filename = "log_";
    filename.append(str);
    filename.append(".txt");
    for(int i = 0; i < filename.length(); ++i){
        if (filename[i] == '/' || filename[i] == ':')
            filename[i] = '-';
    }

    log.open((filename).c_str());

    if (log.fail()) {
        std::cout << "Creating log file failed: "
                  << (filename).c_str()
                  << std::endl;
        exit(1);
    }

    elapsedTime = 0;

    logMsg("Logging started at time string " + str);
}

Debug::~Debug()
{
    logMsg("Logging ended");
    log.close();
}

void Debug::Info(const std::string msg)
{
    if (debugLevel > 2)
    {
        logMsg("I: " + msg);
    }

}

void Debug::Warning(const std::string msg)
{
    if (debugLevel > 1)
    {
        logMsg("W: " + msg);
    }
}

void Debug::Error(const std::string msg)
{
    if (debugLevel > 0)
    {
    logMsg("E: " + msg);
    }
    log.close();
    exit(1);
}

void Debug::Tick()
{
    elapsedTime = 0;
    DebugTimer::tick();
}

void Debug::Tock(std::string identifier)
{
    unsigned long int newTimeDiff = DebugTimer::tick();

    elapsedTime += newTimeDiff;

    logMsg("I: Debug Timer: Time spend at \"" + identifier + "\" was: " + STR(elapsedTime));
}

void Debug::logMsg(const std::string msg)
{
    std::cout << msg << std::endl;
    log << msg << "\n";
}