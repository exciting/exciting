#include <memory>
#include "mdarray.h"


int main(int argc, char **argv)
{
    std::vector<mdarray<int,2>*> v;
    for (int i=0; i<4;i++) {
        v.push_back(new mdarray<int,2>());
        v[i]->set_dimensions(20,20);
        v[i]->allocate();
    }
    
    for (int i =0; i<v.size(); i++) delete v[i];
}

