#ifndef _SIMPLEPROGRESS_HPP_
#define _SIMPLEPROGRESS_HPP_

#include <iostream>
#include <iomanip>

class SimpleProgress
{
private:
    int *_steps;
    int _nextStep;
    int _min, _max;
    
public:
    SimpleProgress(int min, int max, int nsteps)
    {
        _min = min;
        _max = max;
        
        _steps = new int[nsteps];
        
        for(int i=0; i<nsteps; i++)
            _steps[i] = min + (i+1)*(max-min)/nsteps;
        
        _nextStep = 0;
    }
    
    ~SimpleProgress()
    {
        delete[] _steps;
    }
    
    
    inline void printIfHitNext(int current)
    {
        if(current == _steps[_nextStep])
        {
            std::cout << std::fixed << std::setprecision(1) << (double)(current-_min)/(double)(_max-_min)*100.0 << "%   " << std::flush;
            _nextStep++;
        }
    }
    
};


#endif
