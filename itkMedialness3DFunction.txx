#include <cmath>
#include <limits>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <numeric>

#include "itkMedialness3DFunction.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkCrossHelper.h"

namespace itk
{


template < typename TInputImage, typename TCoordRep >
Medialness3DFunction<TInputImage, TCoordRep>::Medialness3DFunction()
{
    typename VectorLinearInterpolateImageFunction<TInputImage>::Pointer defaultIterpolater = VectorLinearInterpolateImageFunction<TInputImage>::New();
    defaultIterpolater->Register();
    this->gradientInterpolator = defaultIterpolater;
    
}

template < typename TInputImage, typename TCoordRep >
Medialness3DFunction<TInputImage, TCoordRep>::~Medialness3DFunction()
{
    this->gradientInterpolator = 0;
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>::
SetInputImage(const TInputImage * ptr)
{
    Superclass::SetInputImage(ptr);
    this->gradientInterpolator->SetInputImage(this->m_Image);

    typename TInputImage::SpacingType pixelSpacing = this->m_Image->GetSpacing();
    typedef std::vector<TCoordRep> SpacingListType;
    SpacingListType spacingList;
    spacingList.push_back(pixelSpacing[0]);
    spacingList.push_back(pixelSpacing[1]);
    spacingList.push_back(pixelSpacing[2]);
    
    this->minSpacing = *std::min_element(spacingList.begin(), spacingList.end(), std::less_equal<TCoordRep>() );
    this->maxSpacing = *std::max_element(spacingList.begin(), spacingList.end(), std::less_equal<TCoordRep>() );
    
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>::
SetInterpolater(VectorInterpolaterType* i)
{
    this->gradientInterpolator = i;
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>::
SetRadii(TCoordRep start, TCoordRep stop, unsigned int numRadii, 
unsigned int maxNumAngles)
{
    this->radii.clear();
    for(int i = 0; i < numRadii; i++)
    {
        TCoordRep currentRadius = start + (TCoordRep(i)/TCoordRep(numRadii))*(stop-start);
        this->radii.push_back(currentRadius);
    } 
    
    this->maxNumAngles = maxNumAngles;
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>::
SetRadii(ListType radii,  unsigned int maxNumAngles)
{
    this->radii = radii;
    std::sort(this->radii.begin(), this->radii.end());
    this->maxNumAngles = maxNumAngles;
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>::
SetNormal(const VectorType& normal)
{
    this->normal = normal;
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>::
ComputeInPlaneVector(VectorType& inPlane1, VectorType& inPlane2) const
{
    TCoordRep n1 = this->normal[0];
    TCoordRep n2 = this->normal[1];
    TCoordRep n3 = this->normal[2];
    
    //Get ready to find orthogonal vector to normal
    //Do this in a numerically stable way
    VectorType uTilde;
    uTilde.Fill(0.0);
 
    if(fabs(n3) <= fabs(n2) && fabs(n3) <= fabs(n1))
    {
        uTilde[2] = 1.0;
    }
    else if(fabs(n2) <= fabs(n3) && fabs(n2) <= fabs(n1))
    {
        uTilde[1] = 1.0;
    }
    else if(fabs(n1) <= fabs(n3) && fabs(n1) <= fabs(n2))
    {
        uTilde[0] = 1.0;
    }
    else
    {
        //throw("Problem with Normal Vector");
    }
    
    //Do Gram Schmidt orthogonalization
    TCoordRep proj = uTilde * this->normal;
    inPlane1 = uTilde - proj*this->normal;
    inPlane1.Normalize();

    typedef CrossHelper<VectorType> CrossType;
    CrossType cross;
    inPlane2 = cross(this->normal, inPlane1);
}

template < typename TInputImage, typename TCoordRep >
typename Medialness3DFunction<TInputImage, TCoordRep>::TOutput Medialness3DFunction<TInputImage, TCoordRep>
::EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const
{
    TCoordRep medialness = 0.0;
    TCoordRep medialnessNorm = 0.0;
    TCoordRep radius = 0.0;

    for(int i = 0; i < this->radii.size(); i++)
    {
        TCoordRep currentRadius = this->radii[i];
        int numAngles = int(ceil( 2.0*M_PI * (currentRadius/this->minSpacing) + 2));
        
        if(numAngles > this->maxNumAngles)
        {
            numAngles = this->maxNumAngles;
        }
        
        TCoordRep currentMedialness = 0.0;
        TCoordRep currentMedialnessNorm = 0.0;
       
        
        for(int k = 0; k < numAngles && numAngles >= 4; k++)
        {
            TCoordRep currentAngle = 2.0*M_PI * (double(k) / double(numAngles));
            VectorType currentPointVector = currentRadius*(inPlane1 *cos(currentAngle) + inPlane2 * sin(currentAngle));

            VectorType currentNormalVector = -1.0*inPlane1 *cos(currentAngle) - 1.0*inPlane2 * sin(currentAngle);
            currentNormalVector.Normalize();

            PointType currentPoint = point + currentPointVector;
            
            if(this->gradientInterpolator->IsInsideBuffer(currentPoint) == true )
            {
                VectorType currentGradient = this->gradientInterpolator->Evaluate(currentPoint);

                if(currentGradient.GetNorm() > 1.0e-6)
                {
                    TCoordRep gradNorm = currentGradient.GetNorm();
                    TCoordRep currentDot = currentNormalVector*currentGradient;                   
                    TCoordRep currentDotNorm = currentDot/gradNorm;

                    currentMedialness += currentDot;
                    currentMedialnessNorm += currentDotNorm;
                   
                }
            }
        }
        currentMedialness = currentMedialness / double(numAngles);
        currentMedialnessNorm = currentMedialnessNorm / double(numAngles);
        
        if(currentMedialness > medialness)
        {
            medialness = currentMedialness;
            medialnessNorm = currentMedialnessNorm;
            radius = currentRadius - this->maxSpacing;
            if(radius < this->minSpacing/2.0 )
            {
                radius = this->minSpacing/2.0;
            }
        }
    }

    
    TOutput medRadPair;
    medRadPair.first = medialnessNorm;
    medRadPair.second = radius;
    
    return(medRadPair);
}

template < typename TInputImage, typename TCoordRep >
typename Medialness3DFunction<TInputImage, TCoordRep>::TOutput Medialness3DFunction<TInputImage, TCoordRep>
::EvaluateWithNormal (const PointType &point, const VectorType &normal) 
{
    this->SetNormal(normal);
    VectorType inPlane1;
    VectorType inPlane2;
    this->ComputeInPlaneVector(inPlane1, inPlane2);
    return(this->EvaluateWithPlane(point, inPlane1, inPlane2));
}

template < typename TInputImage, typename TCoordRep >
typename Medialness3DFunction<TInputImage, TCoordRep>::TOutput Medialness3DFunction<TInputImage, TCoordRep>
::Evaluate(const PointType &point) const
{
    VectorType inPlane1;
    VectorType inPlane2;
    this->ComputeInPlaneVector(inPlane1, inPlane2);
    return(this->EvaluateWithPlane(point, inPlane1, inPlane2));
}

template < typename TInputImage, typename TCoordRep >
typename Medialness3DFunction<TInputImage, TCoordRep>::TOutput Medialness3DFunction<TInputImage, TCoordRep>::
EvaluateAtContinuousIndex (const ContinuousIndexType &index) const
{
    PointType point;
    this->m_Image->TransformContinuousIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}

template < typename TInputImage, typename TCoordRep >
typename Medialness3DFunction<TInputImage, TCoordRep>::TOutput Medialness3DFunction<TInputImage, TCoordRep>::
EvaluateAtIndex (const IndexType &index) const
{
    PointType point;
    this->m_Image->TransformIndexToPhysicalPoint( index, point);
    return(this->Evaluate(point));
}

template < typename TInputImage, typename TCoordRep >
void Medialness3DFunction<TInputImage, TCoordRep>
::PrintSelf(std::ostream& os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);
    os << "Radii ";
    std::ostream_iterator<TCoordRep> osIterNeigh(os, " ");
    std::copy(this->radii.begin(), this->radii.end(), osIterNeigh);
    os << std::endl;
    os << "Max Number of Angles " << this->maxNumAngles << std::endl;
    os << "Normal " << this->normal << std::endl;
}

template < typename TInputImage, typename TCoordRep >
MedialnessMin3DFunction<TInputImage, TCoordRep>::MedialnessMin3DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
MedialnessMin3DFunction<TInputImage, TCoordRep>::~MedialnessMin3DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
typename MedialnessMin3DFunction<TInputImage, TCoordRep>::TOutput MedialnessMin3DFunction<TInputImage, TCoordRep>
::EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const
{
    TCoordRep medialness = 0.0;
    TCoordRep medialnessNorm = 0.0;
    TCoordRep radius = 0.0;

    for(int i = 0; i < this->radii.size(); i++)
    {
        TCoordRep currentRadius = this->radii[i];
        int numAngles = int(ceil( M_PI * (currentRadius/this->minSpacing) + 2));
        
        if(numAngles > this->maxNumAngles)
        {
            numAngles = this->maxNumAngles;
        }
        
        TCoordRep currentMedialness = 0.0;
        TCoordRep currentMedialnessNorm = 0.0;
        
        //std::cout << "Medialness radius " << currentRadius << std::endl;
        //std::cout << "Medialness angles " << numAngles << std::endl;;
        
        for(int k = 0; k < numAngles; k++)
        {
            TCoordRep currentAngle = M_PI * (double(k) / double(numAngles));
            VectorType currentPointVector = currentRadius*(inPlane1 *cos(currentAngle) + inPlane2 * sin(currentAngle));

            VectorType currentNormalVector = -1.0*inPlane1 *cos(currentAngle) - 1.0*inPlane2 * sin(currentAngle);
            currentNormalVector.Normalize();

            PointType currentPointPos = point + currentPointVector;
            PointType currentPointNeg = point - currentPointVector;
            
            if(this->gradientInterpolator->IsInsideBuffer(currentPointPos) == true && this->gradientInterpolator->IsInsideBuffer(currentPointNeg) == true)
            {
                VectorType currentGradientPos = this->gradientInterpolator->Evaluate(currentPointPos);
                VectorType currentGradientNeg = this->gradientInterpolator->Evaluate(currentPointNeg);

                if(currentGradientPos.GetNorm() > 1.0e-6 && currentGradientNeg.GetNorm() > 1.0e-6 )
                {
                    TCoordRep gradPosNorm = currentGradientPos.GetNorm();
                    TCoordRep gradNegNorm = currentGradientNeg.GetNorm();
                    
                    TCoordRep maxGrad = std::max(gradPosNorm, gradNegNorm);
                    TCoordRep minGrad = std::min(gradPosNorm, gradNegNorm);
                    TCoordRep gradRatio = maxGrad/minGrad;
                    
                    TCoordRep currentDotPos = currentNormalVector*currentGradientPos;
                    TCoordRep currentDotNeg = currentNormalVector*currentGradientNeg;
                    TCoordRep currentDotGrad = currentGradientNeg*currentGradientPos;
                    
                    TCoordRep currentDotPosNorm = currentDotPos/gradPosNorm;
                    TCoordRep currentDotNegNorm = currentDotNeg/gradNegNorm;
                    TCoordRep currentDotGradNorm = currentDotGrad/(gradNegNorm*gradPosNorm);
                    
                    //std::cout << "Radius " << currentRadius << " Angle " << currentAngle << std::endl;
                    //std::cout << "Dot Pos " << currentDotPosNorm << std::endl;
                    //std::cout << "Dot Neg " << currentDotNegNorm << std::endl;
                    //std::cout << "Dot Grad " << currentDotGradNorm << std::endl;
                    //std::cout << "Grad Ratio " << gradRatio << std::endl;
                    
                    //if(currentDotPosNorm > 0.5 && currentDotNegNorm < -0.5 && currentDotGradNorm < -0.5 &&  gradRatio < 2)
                    //{
                        //currentMedialness += std::min(currentDotPos, -currentDotNeg) ;
                        //currentMedialnessNorm += (std::min(currentDotPosNorm, -currentDotNegNorm)  - currentDotGradNorm + 1.0/gradRatio)/3.0;
                        //currentMedialness += std::min(currentDotPos, -currentDotNeg)  - currentDotGrad;
                        
                        currentMedialnessNorm += std::min(currentDotPosNorm, -currentDotNegNorm);
                        currentMedialness += std::min(currentDotPos, -currentDotNeg) ;
                    //}
                }
            }
        }
        currentMedialness = currentMedialness / double(numAngles);
        currentMedialnessNorm = currentMedialnessNorm / double(numAngles);
        
        if(currentMedialness > medialness)
        {
            medialness = currentMedialness;
            medialnessNorm = currentMedialnessNorm;
            radius = currentRadius - this->maxSpacing;
            if(radius < this->minSpacing/2.0 )
            {
                radius = this->minSpacing/2.0;
            }
        }
        
        //std::cout << "Point " << point << " Medialness " << currentMedialness << " Medialness Norm " << currentMedialnessNorm <<  " Radius " << currentRadius << std::endl;
    }

    
    TOutput medRadPair;
    medRadPair.first = medialnessNorm;
    //medRadPair.first = medialness;
    medRadPair.second = radius;
    
    return(medRadPair);
}

template < typename TInputImage, typename TCoordRep >
MedialnessVarRadius3DFunction<TInputImage, TCoordRep>::MedialnessVarRadius3DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
MedialnessVarRadius3DFunction<TInputImage, TCoordRep>::~MedialnessVarRadius3DFunction()
{
}


template < typename TInputImage, typename TCoordRep >
typename MedialnessVarRadius3DFunction<TInputImage, TCoordRep>::TOutput MedialnessVarRadius3DFunction<TInputImage, TCoordRep>
::EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const
{
    TCoordRep medialness = 0.0;
    TCoordRep medialnessNorm = 0.0;
    TCoordRep radius = 0.0;
    
    //compute max radius
    //decide on number of angles
    //loop through angles and then radius
    //find radius with maximum response
    //collect radius and response
    //see if all maximum response radii are within ratio
    //if so report average radius
    
    
    int numAngles = this->maxNumAngles;
    
    typedef std::vector<TCoordRep> ListType;
    ListType maxRadiusList(numAngles);
    
    for(int k = 0; k < numAngles && numAngles >= 4; k++)
    {
        TCoordRep currentAngle = 2.0*M_PI * (double(k) / double(numAngles));
        VectorType currentNormalVector = -1.0*inPlane1 *cos(currentAngle) - 1.0*inPlane2 * sin(currentAngle);
        currentNormalVector.Normalize();
            
        TCoordRep largestDot = 0.0;
        TCoordRep largestDotNorm = 0.0;
        TCoordRep largestRadius = 0.0;
        for(int i = 0; i < this->radii.size(); i++)
        {
            TCoordRep currentRadius = this->radii[i];
            VectorType currentPointVector = currentRadius*(inPlane1 *cos(currentAngle) + inPlane2 * sin(currentAngle));

            PointType currentPoint = point + currentPointVector;

            if(this->gradientInterpolator->IsInsideBuffer(currentPoint) == true )
            {
                VectorType currentGradient = this->gradientInterpolator->Evaluate(currentPoint);

                if(currentGradient.GetNorm() > 1.0e-6)
                {
                    TCoordRep currentDot = currentNormalVector*currentGradient;
                    TCoordRep gradNorm = currentGradient.GetNorm();                   
                    TCoordRep currentDotNorm = currentDot/gradNorm;
                    
                    if(currentDot>largestDot)
                    {
                        largestDot = currentDot;
                        largestDotNorm = currentDotNorm;
                        largestRadius = currentRadius;
                    }

                   
                }
            }
        }
        medialness += largestDot;
        medialnessNorm += largestDotNorm;
        maxRadiusList[k] = largestRadius;
    }
    
        
    TCoordRep maxRadius = *std::max_element(maxRadiusList.begin(), maxRadiusList.end(), std::less_equal<TCoordRep>() );
    TCoordRep minRadius = *std::min_element(maxRadiusList.begin(), maxRadiusList.end(), std::less_equal<TCoordRep>() );
    TCoordRep averageRadius = std::accumulate(maxRadiusList.begin(), maxRadiusList.end(), 0 );
    averageRadius = averageRadius/TCoordRep(numAngles);
    
    TCoordRep radiusRatio = 0.0;
    if(minRadius == maxRadius && maxRadius <= 0.0)
    {
        radiusRatio = 1.0;
    }
    else if(minRadius != maxRadius && maxRadius <= 0.0)
    {
        radiusRatio = 0.0;
    }
    else
    {
        radiusRatio = minRadius/maxRadius;
    }

    medialnessNorm = radiusRatio*medialnessNorm/TCoordRep(numAngles);
    //medialnessNorm = medialnessNorm/TCoordRep(numAngles);
    
    TOutput medRadPair;
    medRadPair.first = medialnessNorm;
    medRadPair.second = averageRadius;
    
    return(medRadPair);
}

template < typename TInputImage, typename TCoordRep >
MedialnessNormal3DFunction<TInputImage, TCoordRep>::MedialnessNormal3DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
MedialnessNormal3DFunction<TInputImage, TCoordRep>::~MedialnessNormal3DFunction()
{
}

template < typename TInputImage, typename TCoordRep >
typename MedialnessNormal3DFunction<TInputImage, TCoordRep>::TOutput MedialnessNormal3DFunction<TInputImage, TCoordRep>
::EvaluateWithPlane(const PointType &point, const VectorType &inPlane1, const VectorType &inPlane2) const
{
    TCoordRep medialness = 0.0;
    TCoordRep medialnessNorm = 0.0;
    TCoordRep radius = 0.0;

    for(int i = 0; i < this->radii.size(); i++)
    {
        TCoordRep currentRadius = this->radii[i];
        int numAngles = int(ceil( 2.0*M_PI * (currentRadius/this->minSpacing) + 2));
        
        if(numAngles > this->maxNumAngles)
        {
            numAngles = this->maxNumAngles;
        }
        
        TCoordRep currentMedialness = 0.0;
        TCoordRep currentMedialnessNorm = 0.0;
               
        for(int k = 0; k < numAngles && numAngles >= 4; k++)
        {
            TCoordRep currentAngle = 2.0*M_PI * (double(k) / double(numAngles));
            VectorType currentPointVector = currentRadius*(inPlane1 *cos(currentAngle) + inPlane2 * sin(currentAngle));

            PointType currentPoint = point + currentPointVector;
            
            if(this->gradientInterpolator->IsInsideBuffer(currentPoint) == true )
            {
                VectorType currentGradient = this->gradientInterpolator->Evaluate(currentPoint);

                if(currentGradient.GetNorm() > 1.0e-6)
                {
                    TCoordRep gradNorm = currentGradient.GetNorm();
                    TCoordRep currentDot = this->normal*currentGradient;                   
                    TCoordRep currentDotNorm = currentDot/gradNorm;

                    currentMedialness += 1.0 - fabs(currentDot);
                    currentMedialnessNorm += 1.0 - fabs(currentDotNorm);
                   
                }
            }
        }
        currentMedialness = currentMedialness / double(numAngles);
        currentMedialnessNorm = currentMedialnessNorm / double(numAngles);
        
        if(currentMedialnessNorm > medialness)
        {
            medialness = currentMedialness;
            medialnessNorm = currentMedialnessNorm;
            radius = currentRadius - this->maxSpacing;
            if(radius < this->minSpacing/2.0 )
            {
                radius = this->minSpacing/2.0;
            }
        }
    }

    
    TOutput medRadPair;
    medRadPair.first = medialnessNorm;
    medRadPair.second = radius;
    
    return(medRadPair);
}



}//end namespace itk


