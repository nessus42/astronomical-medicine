/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkLogImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2006/03/19 04:36:56 $
  Version:   $Revision: 1.20 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMakeNanMaskFilter_h
#define __itkMakeNanMaskFilter_h

#include "itkUnaryFunctorImageFilter.h"
#include "vnl/vnl_math.h"

namespace itk
{
  
/** \class MakeNaNMaskFilter
 * \brief Computes the vcl_log(x) pixel-wise
 *
 * \ingroup IntensityImageFilters  Multithreaded
 */
namespace Function {  
  
template< class TInput, class TOutput>
class IsNaN
{
public:
  IsNaN() {}

  ~IsNaN() {}

  bool operator!=( const IsNaN & ) const
  {
    return false;
  }

  bool operator==( const IsNaN & other ) const
  {
    return !(*this != other);
  }

  inline TOutput operator()( const TInput & A )
  {
    return (TOutput)std::isnan(A);
  }
}; // end class IsNan

} // end namespace Function

template <class TInputImage, class TOutputImage>
class ITK_EXPORT MakeNaNMaskFilter :
    public
UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                        Function::IsNaN< typename TInputImage::PixelType, 
                                       typename TOutputImage::PixelType>   >
{
public:
  /** Standard class typedefs. */
  typedef MakeNaNMaskFilter  Self;
  typedef UnaryFunctorImageFilter<TInputImage,TOutputImage, 
                                  Function::IsNaN< typename TInputImage::PixelType, 
                                                 typename TOutputImage::PixelType> > Superclass;
  typedef SmartPointer<Self>   Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputConvertibleToDoubleCheck,
    (Concept::Convertible<typename TInputImage::PixelType, double>));
  itkConceptMacro(DoubleConvertibleToOutputCheck,
    (Concept::Convertible<double, typename TOutputImage::PixelType>));
  /** End concept checking */
#endif

protected:
  MakeNaNMaskFilter() {}
  virtual ~MakeNaNMaskFilter() {}

private:
  MakeNaNMaskFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namespace itk


#endif
