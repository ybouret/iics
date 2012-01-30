#ifndef ICP_ARRAYS_INCLUDED
#define ICP_ARRAYS_INCLUDED 1
#include "types.h"

#ifdef __cplusplus
extern "C" {
#endif
    


size_t icp_bytes(); /*!< current allocated bytes */

/*******************************************************************************
 ** 1D
 ******************************************************************************/

/*! create a new 1D array */
/**
 \param ilo lower index
 \param ihi upper index
 \return a valid array upon success, NULL upon failure.
 
 a = icp_create_array1D(ilo,ihi);
 a[ilo]=...; a[ihi]=...
 */
real_t *icp_create_array1D( indx_t  ilo, indx_t ihi );

/*! delete a 1D array */
void    icp_delete_array1D( real_t *a, indx_t ilo, indx_t ihi );

/*******************************************************************************
 ** 2D
 ******************************************************************************/
/*! create a new 2D array */
/**
 \param xmin lower x index
 \param xmax upper x index
 \param ymin lower y index
 \param ymax upper y index
 \return a valid array upon success, NULL upon failure
 \warning be careful about the index order !
 
 a = icp_create_array2D(xmin,xmax,ymin,ymax);
 for( j=ymin;j<=ymax;++j)
 {
 for( i=xmin;i<=xmax;++i)
 {
 a[j][i] = ...;
 }
 }
 
 */
real_t **icp_create_array2D( indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax );

/*! delete a 2D array */
void     icp_delete_array2D( real_t **a,  indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax );

/*******************************************************************************
 ** 3D
 ******************************************************************************/
real_t ***icp_create_array3D( indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax, indx_t zmin, indx_t zmax );
void      icp_delete_array3D( real_t ***a, indx_t xmin, indx_t xmax, indx_t ymin, indx_t ymax, indx_t zmin, indx_t zmax );
#ifdef __cplusplus
}
#endif
#endif

