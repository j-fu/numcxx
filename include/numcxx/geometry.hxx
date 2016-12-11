#ifndef NUMCXX_GEOMETRY_H
#define NUMCXX_GEOMETRY_H

#include "tarray1.hxx"
#include "tarray2.hxx"

namespace numcxx
{
    class Geometry
    {
        
    public:
        void set_points       (const std::initializer_list<std::initializer_list<double>> &il) {points=TArray2<double>::create(il);};
        void set_bfaces       (const std::initializer_list<std::initializer_list<int>> &il) {bfaces=TArray2<int>::create(il);};
        void set_bfaceregions (const std::initializer_list<int> &il){bfaceregions=TArray1<int>::create(il);};
        void set_regionpoints (const std::initializer_list<std::initializer_list<double>> &il) {regionpoints=TArray2<double>::create(il);};
        void set_regionnumbers(const std::initializer_list<int> &il){regionnumbers=TArray1<int>::create(il);};
        void set_regionvolumes(const std::initializer_list<double> &il){regionvolumes=TArray1<double>::create(il);};
        
        std::shared_ptr<TArray2<double>> points=nullptr;
        std::shared_ptr<TArray2<int>>    bfaces=nullptr;
        std::shared_ptr<TArray1<int>>    bfaceregions=nullptr;
        std::shared_ptr<TArray2<double>> regionpoints=nullptr;
        std::shared_ptr<TArray1<int>>    regionnumbers=nullptr;
        std::shared_ptr<TArray1<double>> regionvolumes=nullptr;
        static std::shared_ptr<Geometry> create() { return std::make_shared<Geometry>();}
    };
    
}


#endif
