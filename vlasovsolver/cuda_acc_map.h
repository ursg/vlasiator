bool map3DCuda(Realf **blockDatas,
               vmesh::GlobalID **blockIDs,
               vmesh::LocalID *nBlocks,
               Real *intersections,
               const uint nCells,
               const Realf blockSize[3],
               const  vmesh::LocalID gridLength[3],
               const Realf gridMinLimits[3],
               const float sparsityThreshold);





/*! A namespace for storing indices into an array which contains 
 * the intersections for accelerations.*/
namespace AccelerationIntersections {
   enum {
      X,    /*!< Map in x-coordinate. Intersection x-coordinate at i,j,k=0 */
      X_DI, /*!<Map in x-coordinate. Change in x-coordinate for a change in i index of 1*/
      X_DJ, /*!<Map in x-coordinate. Change in x-coordinate for a change in j index of 1*/
      X_DK, /*!<Map in x-coordinate. Change in x-coordinate for a change in k index of 1*/
      Y,    /*!< Map in y-coordinate. Intersection y-coordinate at i,j,k=0 */
      Y_DI, /*!<Map in y-coordinate. Change in y-coordinate for a change in i index of 1*/
      Y_DJ, /*!<Map in y-coordinate. Change in y-coordinate for a change in j index of 1*/
      Y_DK, /*!<Map in y-coordinate. Change in y-coordinate for a change in k index of 1*/
      Z,    /*!< Map in z-coordinate. Intersection z-coordinate at i,j,k=0 */
      Z_DI, /*!<Map in z-coordinate. Change in z-coordinate for a change in i index of 1*/
      Z_DJ, /*!<Map in z-coordinate. Change in z-coordinate for a change in j index of 1*/
      Z_DK, /*!<Map in z-coordinate. Change in z-coordinate for a change in k index of 1*/
      N_INTERSECTIONS
   };
}
