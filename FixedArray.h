#ifndef __FixedArray_h
#define __FixedArray_h

#include <iostream>


/**\class FixedArray
 * A class representing an array of fixed length of specified type.
 */ 
template<typename TValueType, unsigned int VDimension>
class FixedArray
{
    public:      
    
        FixedArray();
        FixedArray(const TValueType v[VDimension]);
        FixedArray(const FixedArray& fa);
        virtual ~FixedArray();
        void Fill(const TValueType&);
        
        FixedArray& operator=(const TValueType v[VDimension]);
        FixedArray& operator=(const FixedArray& fa);

        bool operator==(const FixedArray& fa) const;
        bool operator!=(const FixedArray& fa) const;
     
        TValueType& operator[](short index)                { return array[index]; }
        const TValueType& operator[](short index) const          { return array[index]; }
        TValueType& operator[](unsigned short index)       { return array[index]; }
        const TValueType& operator[](unsigned short index) const { return array[index]; }
        TValueType& operator[](int index)                  { return array[index]; }
        const TValueType& operator[](int index) const            { return array[index]; }
        TValueType& operator[](unsigned int index)         { return array[index]; }
        const TValueType& operator[](unsigned int index) const   { return array[index]; }
        TValueType& operator[](long index)                 { return array[index]; }
        const TValueType& operator[](long index) const           { return array[index]; }
        TValueType& operator[](unsigned long index)        { return array[index]; }
        const TValueType& operator[](unsigned long index) const  { return array[index]; }
        
        FixedArray operator-() const;
        FixedArray operator+(const FixedArray& fa) const;
        FixedArray operator-(const FixedArray& fa) const;
        FixedArray operator*(const TValueType& s) const;
        
        std::ostream& Write(std::ostream& os) const;
        std::istream& Read(std::istream& is);
        
        double Norm() const;
        FixedArray Normalize() const;
        double Distance(const FixedArray& fa) const;
        
    private:
        TValueType array[VDimension];
};

template<typename TValueType, unsigned int VDimension>
std::ostream & operator<<(std::ostream &os, const FixedArray<TValueType, VDimension> &fa);

template<typename TValueType, unsigned int VDimension>
std::istream & operator>>(std::istream &is, FixedArray<TValueType, VDimension> &fa);

template<typename TValueType, unsigned int VDimension>
double dot_product(const FixedArray<TValueType, VDimension> &a, const FixedArray<TValueType, VDimension> &b);

template<typename TValueType>
FixedArray<TValueType, 3> cross_product(const FixedArray<TValueType, 3> &a, const FixedArray<TValueType, 3> &b);

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension> operator*(const TValueType a, const FixedArray<TValueType, VDimension> &b);


#include "FixedArray.txx"

#endif
