#ifndef __Interpolater_h
#define __Interpolater_h

#include <iostream>
#include <vector>

#include "FixedArray.h"


/**\class Interpolater
 *  This class which implements spline interpolation
 * see "The NURBS Book" by Les Piegl and Wayne Tiller
 */ 
template<typename TValueType, unsigned int VDimension>
class Interpolater
{
    public:      
        typedef std::vector<TValueType> TValueListType;
        typedef std::vector<TValueListType> TValueMatrixType;
        typedef FixedArray<TValueType, VDimension> PointType;
        typedef FixedArray<TValueType, VDimension> VectorType;
        typedef std::vector<PointType> PointListType;
    
        /**
        * Determine the knot span index (see Algorithm A2.1)
        * @param numControlPoints number of control points
        * @param order spline order
        * @param t where on curve to evaluate in [0,1]
        * @param knotVector spline knots
        * @return span index
        */
        int GetSpanIndex(int numControlPoints, int order, TValueType t, const TValueListType& knotVector) const;
        
        /**
        * Evaluate basis functions (see Algorithm A2.2)
        * @param spanIndex span index
        * @param order spline order
        * @param t where on curve to evaluate in [0,1]
        * @param knotVector spline knots
        * @return a point
        */
        TValueListType EvaluateBasisFunctions(int spanIndex, int order, TValueType t, const TValueListType& knotVector) const;
        
        /**
        * Evaluate derivative of basis functions (see Algorithm A2.3)
        * @param spanIndex span index
        * @param order spline order
        * @param t where on curve to evaluate in [0,1]
        * @param knotVector spline knots
        * @return a point
        */
        TValueMatrixType EvaluateDerivativeBasisFunctions(int spanIndex, int order, int numDerivatives, TValueType t, const TValueListType& knotVector) const;
        
        /**
        * Evaluate point on the curve defined by the following (see Algorithm A4.1)
        * @param numControlPoints number of control points
        * @param order spline order
        * @param t where on curve to evaluate in [0,1]
        * @param knotVector spline knots
        * @param controlPoints control points
        * @return point at t
        */
        PointType EvaluateCurvePoint(int numControlPoints, int order, TValueType t, const TValueListType& knotVector, const PointListType& controlPoints) const;
        
        /**
        * Evaluate point on the derivative of the curve defined by the following 
        * @param numControlPoints number of control points
        * @param order spline order
        * @param t where on curve to evaluate in [0,1]
        * @param knotVector spline knots
        * @param controlPoints control points
        * @return point at t
        */
        PointType EvaluateDerivatePoint(int numControlPoints, int order, TValueType t, const TValueListType& knotVector, const PointListType& controlPoints) const;
        
        /**
        * Compute points in 3d along a circle at specified point which is in the
        * plane of the normal and binormal (normal cross tangent) vector.  We
        * attempt to prevent twists from centerpoint to centerpoint
        * @param numAngles number of points to compute
        * @param radius the radius from origin to compute circle
        * @param point center of circle
        * @param pointTangent tangent vector at point
        * @param pointNormal normal vector at point
        * @param referenceNormal normal vector at preceeding point
        * @return a list of points in a circle with specified radius and
        *  centered at specified point in the specified plane
        */
        PointListType EvaluateCircle(int numAngles, TValueType radius, const PointType& point, VectorType pointTangent, VectorType pointNormal, PointType& referenceNormal) const;
        
        /**
        * Get spline knots with equal weighting
        * @param order spline order
        * @param numKnots number of knots
        * @return the knots
        */
        TValueListType GetEqualKnots(int order, int numKnots) const;
        
        /**
        * Get spline knots with 0 or 1 weighting
        * @param order spline order
        * @param numKnots number of knots
        * @return the knots
        */
        TValueListType GetEndWeightKnots(int order, int numKnots) const;
        
    private:

};



#include "Interpolater.txx"

#endif
