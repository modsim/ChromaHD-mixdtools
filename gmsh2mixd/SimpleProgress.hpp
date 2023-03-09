#ifndef _SIMPLEPROGRESS_HPP_
#define _SIMPLEPROGRESS_HPP_

#include <iostream>
#include <iomanip>

class SimpleProgress
{
private:
    long *_steps;
    long _nextStep;
    long _min, _max;
    
public:
    SimpleProgress(long min, long max, long nsteps)
    {
        _min = min;
        _max = max;
        
        _steps = new long[nsteps];
        
        for(long i=0; i<nsteps; i++)
            _steps[i] = min + (i+1)*(max-min)/nsteps;
        
        _nextStep = 0;
    }
    
    ~SimpleProgress()
    {
        delete[] _steps;
    }
    
    
    inline void printIfHitNext(long current)
    {
        if(current == _steps[_nextStep])
        {
            std::cout << std::fixed << std::setprecision(1) << (double)(current-_min)/(double)(_max-_min)*100.0 << "%   " << std::flush;
            _nextStep++;
        }
    }
    
};


#endif
