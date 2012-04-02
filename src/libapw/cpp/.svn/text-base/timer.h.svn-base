#include <sys/time.h>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <stdint.h>

class timer
{
    public:
        
        timer(std::string tname);
        
        ~timer();
        
        static void print();
        
        static void print(std::string tname);
        
    private:
    
        std::string tname;
        timeval start;
        static std::map<std::string,double> timers;
        static std::map<std::string,int> tcount;
};

        
    

