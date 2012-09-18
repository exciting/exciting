#include "timer.h"
#include "config.h"

std::map<std::string,double> timer::timers;
std::map<std::string,int> timer::tcount;


timer::timer(std::string tname) : tname(tname)
{
    if (timers.count(tname) == 0) 
    {
        timers[tname] = 0;
        tcount[tname] = 0;
    }
    gettimeofday(&start, NULL);
}
        
timer::~timer()
{
    timeval end;
    gettimeofday(&end, NULL);
    timers[tname] += double(end.tv_sec - start.tv_sec) + double(end.tv_usec - start.tv_usec)/1e6;
    tcount[tname]++;
}
        
void timer::print()
{
    std::map<std::string,double>::iterator it;
    for (it = timers.begin(); it != timers.end(); it++)
        std::cout << it->first << " : " << it->second << " (total) " << it->second/tcount[it->first] << " (average) " << std::endl;
}

void timer::print(std::string tname)
{
    std::cout << timers[tname] << std::endl;
}

std::map<std::string,timer*> ftimers;

extern "C" void FORTRAN(lapw_timer_start)(char *name_, int32_t name_len)
{
    std::string name(name_, name_len);
    ftimers[name] = new timer(name);
}