#include <cmath>

#include "FixedArray.h"

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>::FixedArray()
{

}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>::FixedArray(const TValueType v[VDimension])
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        this->array[i] = v[i];
    }
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>::FixedArray(const FixedArray& fa)
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        this->array[i] = fa.array[i];
    }
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>::~FixedArray()
{

}

template<typename TValueType, unsigned int VDimension>
bool FixedArray<TValueType, VDimension>::operator==(const FixedArray& fa) const
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        if(this->array[i] != fa[i])
        {
            return(false);
        }
    }
    return(true);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>& FixedArray<TValueType, VDimension>::operator=(const TValueType v[VDimension])
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        this->array[i] = v[i];
    }
    return(*this);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>& FixedArray<TValueType, VDimension>::operator=(const FixedArray& fa)
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        this->array[i] = fa.array[i];
    }
    return(*this);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension> FixedArray<TValueType, VDimension>::operator-() const
{
    return(-1.0* (*this));
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension> FixedArray<TValueType, VDimension>::operator*(const TValueType& s) const
{
    FixedArray fa;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        fa[i] = s*this->array[i];
    }
    
    return(fa);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>  FixedArray<TValueType, VDimension>::operator+(const FixedArray& fa) const
{
    FixedArray added;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        added[i] = this->array[i] + fa[i];
    }
    
    return(added);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension>  FixedArray<TValueType, VDimension>::operator-(const FixedArray& fa) const
{
    FixedArray subtracted;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        subtracted[i] = this->array[i] - fa[i];
    }
    
    return(subtracted);
}

template<typename TValueType, unsigned int VDimension>
bool FixedArray<TValueType, VDimension>::operator!=(const FixedArray& fa) const
{
    return( !(*this == fa) );
}



template<typename TValueType, unsigned int VDimension>
void FixedArray<TValueType, VDimension>::Fill(const TValueType& v)
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        this->array[i] = v;
    }
}

template<typename TValueType, unsigned int VDimension>
std::ostream& FixedArray<TValueType, VDimension>::Write(std::ostream& os) const
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        os << this->array[i] << " ";
    }
    
    return(os);
}

template<typename TValueType, unsigned int VDimension>
std::istream& FixedArray<TValueType, VDimension>::Read(std::istream& is) 
{
    for(unsigned int i = 0; i < VDimension; i++)
    {
        TValueType v;
        is >> v;
        this->array[i] = v;
    }
    return(is);
}
template<typename TValueType, unsigned int VDimension>
double FixedArray<TValueType, VDimension>::Norm() const
{
    double norm = 0.0;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        norm += this->array[i]*this->array[i];
    }
    norm = sqrt(norm);
    return(norm);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension> FixedArray<TValueType, VDimension>::Normalize() const
{
    double norm = this->Norm();
    FixedArray normalized(*this);
    if(norm > 0)
    {
        for(unsigned int i = 0; i < VDimension; i++)
        {
            normalized[i] = this->array[i]/norm;
        }
    }
    return(normalized);
}

template<typename TValueType, unsigned int VDimension>
double FixedArray<TValueType, VDimension>::Distance(const FixedArray& fa) const
{
    double dist = 0.0;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        dist += (this->array[i]-fa[i])*(this->array[i]-fa[i]);
    }
    dist = sqrt(dist);
    return(dist);
}


template<typename TValueType, unsigned int VDimension>
std::ostream & operator<<(std::ostream &os, const FixedArray<TValueType, VDimension> &fa)
{
    return(fa.Write(os));
}

template<typename TValueType, unsigned int VDimension>
std::istream & operator>>(std::istream &is, FixedArray<TValueType, VDimension> &fa)
{
    return(fa.Read(is));
}

template<typename TValueType, unsigned int VDimension>
double dot_product(const FixedArray<TValueType, VDimension> &a, const FixedArray<TValueType, VDimension> &b)
{
    double dot = 0.0;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        dot += a[i]*b[i];
    }
    return(dot);
}

template<typename TValueType>
FixedArray<TValueType, 3> cross_product(const FixedArray<TValueType, 3> &a, const FixedArray<TValueType, 3> &b)
{
    FixedArray<TValueType, 3> cross;
    
    cross[0] = a[1]*b[2] - a[2]*b[1];
    cross[1] = a[2]*b[0] - a[0]*b[2];
    cross[2] = a[0]*b[1] - a[1]*b[0];
    
    return(cross);
}

template<typename TValueType, unsigned int VDimension>
FixedArray<TValueType, VDimension> operator*(const TValueType a, const FixedArray<TValueType, VDimension> &b)
{
    FixedArray<TValueType, VDimension> fa;
    for(unsigned int i = 0; i < VDimension; i++)
    {
        fa[i] = a*b[i];
    }
    
    return(fa);
}


