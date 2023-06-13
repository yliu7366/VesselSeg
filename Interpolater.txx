#include <cmath>

#include "Interpolater.h"

template<typename TValueType, unsigned int VDimension>
int Interpolater<TValueType, VDimension>::GetSpanIndex(int numControlPoints, int order, TValueType t, const TValueListType& knotVector) const
{
    if(t == knotVector[numControlPoints+1])
    {
        return(numControlPoints);
    }
    int low = order;
    int high = numControlPoints+1;
    int mid = (low+high) / 2;
    
    while(t < knotVector[mid] || t >= knotVector[mid+1])
    {
        //std::cout << "Low " << low << " High " << high << std::endl;
        //std::cout << "Knot at mid " << knotVector[mid] << " Knot at mid + 1 " << knotVector[mid+1] << std::endl;
        if(t < knotVector[mid])
        {
            high = mid;
        }
        else
        {
            low = mid;
        }
        mid = (low+high) / 2;
    }
    return(mid);
}

template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::TValueListType 
Interpolater<TValueType, VDimension>::EvaluateBasisFunctions(int spanIndex, int order, TValueType t, const TValueListType& knotVector) const
{
    TValueListType basis(order+1);
    basis[0] = 1.0;
    TValueListType left(order+1);
    TValueListType right(order+1);
    for(unsigned int k = 1; k <= order; k++)
    {
        left[k] = t - knotVector[spanIndex - k + 1];
        right[k] = knotVector[spanIndex + k] - t;
        
        TValueType previous = 0.0;
        for(unsigned int l = 0; l < k; l++)
        {
            TValueType temp = basis[l] / (right[l+1]+left[k-l]);
            basis[l] = previous + right[l+1]*temp;
            previous = left[k-l]*temp; 
        }
        
        basis[k] = previous;
    }
    
    return(basis);
}

template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::TValueMatrixType 
Interpolater<TValueType, VDimension>::EvaluateDerivativeBasisFunctions(int spanIndex, int order, int numDerivatives, TValueType t, const TValueListType& knotVector) const
{
    TValueListType firstDBasis(order+1);
    for(int k = 0; k <= order; k++)
    {
        firstDBasis[k] = 0.0;
    }
    
    TValueMatrixType ndu;
    for(int k = 0; k <= order; k++)
    {
        ndu.push_back( firstDBasis );
    }
    ndu[0][0] = 1.0;

    TValueMatrixType a;
    a.push_back( firstDBasis);
    a.push_back( firstDBasis);

    TValueMatrixType ders;
    for(int k = 0; k <= numDerivatives; k++)
    {
        ders.push_back( firstDBasis );
    }

    TValueListType left(order+1);
    TValueListType right(order+1);
    for(int k = 1; k <= order; k++)
    {
        left[k] = t - knotVector[spanIndex - k + 1];
        right[k] = knotVector[spanIndex + k] - t;
        
        TValueType previous = 0.0;
        for(int l = 0; l < k; l++)
        {
            ndu[k][l] = right[l+1] + left[k-l];
            
            TValueType temp = ndu[l][k-1]/ndu[k][l];
            ndu[l][k] = previous + right[l+1]*temp;
            previous = left[k-l]*temp; 
            
        }
        ndu[k][k] = previous;
    }
   
    for(int k = 0; k <= order; k++)
    {
        ders[0][k] = ndu[k][order];
    }
    
    for(int l = 0; l <= order; l++)
    {
        int index1 = 0;
        int index2 = 1;
        a[0][0] = 1.0;
        for(int k = 1; k <= numDerivatives; k++)
        {
            TValueType d = 0.0;
            int rk = l-k;
            int pk = order-k;
            
            if(l >= k)
            {
                a[index2][0] = a[index1][0] / ndu[pk+1][rk];
                d = a[index2][0] * ndu[rk][pk];
            }
            int j1 = 0;
            if(rk >= -1)
            {
                j1 = 1;
            }
            else
            {
                j1 = -rk;
            }
            
            int j2 = 0;
            if(l-1 <= pk)
            {
                j2 = k-1;
            }
            else
            {
                j2 = order - l;
            }
            for(int j = j1; j <= j2; j++)
            {
                a[index2][j] = (a[index1][j] - a[index1][j-1]) / ndu[pk+1][rk+j];
                d += a[index2][j]*ndu[rk+j][pk];
            }
            if( l <= pk)
            {
                a[index2][k] = -a[index1][k-1]/ndu[pk+1][l];
                d += a[index2][k]*ndu[l][pk];
            }
            ders[k][l] = d;
            int tempIndex = index1;
            index1 = index2;
            index2 = tempIndex;
        }
    }
    
    int fact = order;
    for(int k = 1; k <= numDerivatives; k++)
    {
        for(int j = 0; j <= order; j++)
        {
            ders[k][j] *= fact;
        }
        fact *= order - k;
    }
    
    return(ders);
}

template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::PointType 
Interpolater<TValueType, VDimension>::EvaluateCurvePoint(int numControlPoints, int order, TValueType t, const TValueListType& knotVector, const PointListType& controlPoints) const
{
    PointType thePoint;
    thePoint.Fill(0.0);
    int spanIndex = this->GetSpanIndex(numControlPoints, order, t, knotVector);
    //std::cout << "Span Index " << spanIndex << " Knot at Index " << knotVector[spanIndex] << std::endl;
    TValueListType basisFunctions = this->EvaluateBasisFunctions(spanIndex, order, t, knotVector);
    for(unsigned int i = 0; i <= order; i++)
    {
        thePoint = thePoint + basisFunctions[i]*controlPoints[spanIndex -order+i];
    }
    
    return(thePoint);
}

template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::PointType 
Interpolater<TValueType, VDimension>::EvaluateDerivatePoint(int numControlPoints, int order, TValueType t, const TValueListType& knotVector, const PointListType& controlPoints) const
{
    PointType thePoint;
    thePoint.Fill(0.0);
    int spanIndex = this->GetSpanIndex(numControlPoints, order, t, knotVector);
    //std::cout << "Span Index " << spanIndex << " Knot at Index " << knotVector[spanIndex] << std::endl;
    int derivativeOrder = 1;
    TValueMatrixType basisFunctions = this->EvaluateDerivativeBasisFunctions(spanIndex, order, derivativeOrder, t, knotVector);
   
    for(unsigned int i = 0; i <= order; i++)
    {
        thePoint = thePoint + basisFunctions[derivativeOrder][i]*controlPoints[spanIndex -order+i];
    }
    
    return(thePoint);
}

template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::PointListType 
Interpolater<TValueType, VDimension>::EvaluateCircle(int numAngles, TValueType radius, const PointType& point, VectorType pointTangent, VectorType pointNormal, VectorType& referenceNormal) const
{
    PointListType circlePoints;
   
    
    //If normal is zero do Gram schmidt to find another normal
    if(pointNormal.Norm() <= 0)
    {
        VectorType temp;
        temp.Fill(0.0);
        if(fabs(pointTangent[2]) <= fabs(pointTangent[1]) && fabs(pointTangent[2]) <= fabs(pointTangent[0]))
        {
            temp[2] = 1.0;
        }
        else if(fabs(pointTangent[1]) <= fabs(pointTangent[2]) && fabs(pointTangent[1]) <= fabs(pointTangent[0]))
        {
            temp[1] = 1.0;
        }
        else 
        {
            temp[0] = 1.0;
        }

        TValueType proj = dot_product(temp, pointTangent);
        VectorType normal = temp - proj*pointTangent;
        normal = normal.Normalize();
        pointNormal = normal;
    }

    PointType pointBinormal = cross_product(pointTangent, pointNormal);
    
    //project currentNormal vector onto reference Normal vector
    TValueType dotX = dot_product(pointNormal, referenceNormal);
    TValueType dotY = dot_product(pointBinormal, referenceNormal);
    
    
    TValueType angleOffset = atan2(dotY,dotX);
    PointType newNormal = pointNormal*cos(angleOffset) + pointBinormal*sin(angleOffset);
    PointType newBiNormal = cross_product(pointTangent, newNormal);
    
    for(int i = 0; i < numAngles; i++)
    {
        TValueType currentAngle = TValueType(i)/TValueType(numAngles) * 2.0 * M_PI;
        PointType currentVector = radius * (newNormal*cos(currentAngle) + newBiNormal*sin(currentAngle));
        PointType currentPoint = point + currentVector;
        
        circlePoints.push_back(currentPoint);
    }
    
    referenceNormal = newNormal;
    return(circlePoints);
}


template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::TValueListType 
Interpolater<TValueType, VDimension>::GetEqualKnots(int order, int numKnots) const
{
    TValueListType knotVector(numKnots);
    for(int i = 0; i <= order; i++)
    {
        knotVector[i] = 0.0;
    }
    int numMiddle = numKnots-2*order-1;
    for(int i = 1; i < numMiddle; i++)
    {
        knotVector[i+order] = TValueType(i)/TValueType(numMiddle);
    }
    for(int i = 0; i <= order; i++)
    {
        knotVector[i+order+numMiddle] = 1.0;
    } 
    
    return(knotVector);
}

template<typename TValueType, unsigned int VDimension>
typename Interpolater<TValueType, VDimension>::TValueListType 
Interpolater<TValueType, VDimension>::GetEndWeightKnots(int order, int numKnots) const
{
    TValueListType knotVector(numKnots);

    for(int i = 0; i <= order; i++)
    {
        knotVector[i] = 0;
    }
    int numMiddle = numKnots-2*order-1;
    for(int i = 1; i < numMiddle/2; i++)
    {
        knotVector[i+order] = 0.0;
    }
    for(int i = numMiddle/2; i < numMiddle; i++)
    {
        knotVector[i+order] = 1.0;
    }
    for(int i = 0; i <= order; i++)
    {
        knotVector[i+order+numMiddle] = 1.0;
    } 
    
    return(knotVector);
}


