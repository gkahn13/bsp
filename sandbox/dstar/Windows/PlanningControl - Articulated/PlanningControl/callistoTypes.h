#pragma once

/*! The result of #CAL_GetObjectType if the object is a box. \ingroup Retrieval */
#define CAL_BOX                        1
/*! The result of #CAL_GetObjectType if the object is a sphere. \ingroup Retrieval */
#define CAL_SPHERE                     2
/*! The result of #CAL_GetObjectType if the object is a cylinder. \ingroup Retrieval */
#define CAL_CYLINDER                   3
/*! The result of #CAL_GetObjectType if the object is a cone. \ingroup Retrieval */
#define CAL_CONE                       4
/*! The result of #CAL_GetObjectType if the object is a triangle group. \ingroup Retrieval */
#define CAL_TRIANGLES                  5
/*! The result of #CAL_GetObjectType if the object is a polygon. \ingroup Retrieval */
#define CAL_POLYGON                    6
/*! The result of #CAL_GetObjectType if the object is a polygon group. \ingroup Retrieval */
#define CAL_POLYGONGROUP               7
/*! The result of #CAL_GetObjectType if the object is a polyline. \ingroup Retrieval */
#define CAL_POLYLINE                   8
/*! The result of #CAL_GetObjectType if the object is a tetrahedron. \ingroup Retrieval */
#define CAL_TETRAHEDRON                9
/*! The result of #CAL_GetObjectType if the object is an elevation grid. \ingroup Retrieval */
#define CAL_ELEVATIONGRID             10
/*! The result of #CAL_GetObjectType if the object is a user drawn object. \ingroup Retrieval */
#define CAL_USERDRAWN                 11

/*! Retrieve information about the last collision check/closest pair/penetration depth. 
\brief Structure to receive collision info.
\ingroup AdvancedGlobal
*/

struct SCALResult
{
  /*! The first resulting object. */
  int objID0;
  /*! The second resulting object. */
  int objID1;
  /*! The group ID of the first resulting object. */
  int groupID0;
  /*! The group ID of the second resulting object. */
  int groupID1;
  /*! The first point involved in the check. Only applicable in case of closest pair and penetration depth. */
  CAL_scalar vector0[3];
  /*! The second point involved in the check. Only applicable in case of closest pair and penetration depth. */
  CAL_scalar vector1[3];
  /*! The penetration distance or distance between the groups, depending on the type of check. */
  CAL_scalar distance;
};

/*! The structure used with CAL_GetGroup to retrieve group information. 
\brief Structure to retrieve a group.
\ingroup Retrieval 
*/
struct SCALGroup
{
  /*! The number of direct subgroups this group has. */
  int nrChildren;
  /*! The number of objects this group contains (not including the potential subgroups). */
  int nrObjects;
  /*! The visibility settings for all of the 4 views. */
  bool visible[4];
  /*! Is the group collision capable (is it in the collision checker)? */
  bool colCapable;
  /*! The current world matrix of the group. */
  CAL_matrix matrix;
  /*! The label of the group. */
  char* label;
};

/*! A structure used with CAL_GetObject to retrieve information about a box. 
\brief Structure to retrieve a box.
\ingroup Retrieval 
*/
struct SCALBox
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The width of the box. */
  CAL_scalar width;
  /*! The height of the box. */
  CAL_scalar height;
  /*! The depth of the box. */
  CAL_scalar depth;
};

/*! A structure used with CAL_GetObject to retrieve information about a sphere.
\brief Structure to retrieve a sphere.
\ingroup Retrieval 
*/
struct SCALSphere
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The radius of the sphere. */
  CAL_scalar radius;
};

/*! A structure used with CAL_GetObject to retrieve information about a cylinder.
\brief Structure to retrieve a cylinder.
\ingroup Retrieval 
*/
struct SCALCylinder
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The radius of the cylinder. */
  CAL_scalar radius;
  /*! The height of the cylinder. */
  CAL_scalar height;
};

/*! A structure used with CAL_GetObject to retrieve information about a cone.
\brief Structure to retrieve a cone.
\ingroup Retrieval 
*/
struct SCALCone
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The radius of the cone. */
  CAL_scalar radius;
  /*! The height of the cone. */
  CAL_scalar height;
};

/*! A structure used with CAL_GetObject to retrieve information about a group of triangles.
\brief Structure to retrieve a triangle mesh.
\ingroup Retrieval 
*/
struct SCALTriangles
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The number of triangles in this object. */
  int nrTriangles;
  /*! The points of the triangles, 9 values for every triangle. */
  CAL_scalar *points;
};

/*! A structure used with CAL_GetObject to retrieve information about a polygon.
\brief Structure to retrieve a polygon.
\ingroup Retrieval 
*/
struct SCALPolygon
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The number of points the polygon contains. */
  int nrPoints;
  /*! The coordinates of the polygon, 3 values for every point. */
  CAL_scalar *points;
};

/*! A structure used with CAL_GetObject to retrieve information about a polygon group.
\brief Structure to retrieve a polygon group.
\ingroup Retrieval 
*/
struct SCALPolygonGroup
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The number of polygons. */
  int nrPolygons;
  /*! The number of points per polygon. The length of this list is determined by nrPolygons. */
  int *nrPoints;
  /*! The coordinates of the polygons, 3 values for every point. */
  CAL_scalar *points;
};

/*! A structure used with CAL_GetObject to retrieve information about a poly line.
\brief Structure to retrieve a polyline.
\ingroup Retrieval 
*/
struct SCALPolyline
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The number of lines. */
  int nrLines;
  /*! The number of points per line. The length of this list is determined by nrLines. */
  int *nrPoints;
  /*! The coordinates of the lines, 3 values per point. */
  CAL_scalar *points;
};

/*! A structure used with CAL_GetObject to retrieve information about a tetrahedron.
\brief Structure to retrieve a tetrahedron.
\ingroup Retrieval 
*/
struct SCALTetrahedron
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The coordinates of the tetrahedron. The length of this list is 3*4=12 values. */
  CAL_scalar *points;
};

/*! A structure used with CAL_GetObject to retrieve information about a user drawn object.
\brief Structure to retrieve a user drawn object.
\ingroup Retrieval 
*/
struct SCALUserDrawn
{
  /*! The ID of the group the object is in. */
  int groupID;
  /*! The current world matrix of the object. */
  CAL_matrix matrix;
  /*! The pointer to the callback function */
  CAL_UserDrawnCallback callback;
};

/*! \file callistoTypes.h
\brief Structures to retrieve data from Callisto.
*/