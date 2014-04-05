#pragma once

// the resulting manual of this file can be found at http://www.cs.uu.nl/people/dennis/callisto/docs/html/index.html

/*! \mainpage Callisto 3.04b
 *
 * \section intro Introduction
 *
 * The Callisto functions are divided in two parts: \ref Basic "basic" and \ref Advanced "advanced".
   The basic functions include functions to create groups and objects, to set their
   transformations and do collision checks. The advanced functions add multiple views,
   extra types of collission checks, clearance, animation and group/object retrieval.
   You are advised to first read the <A href="http://www.cs.uu.nl/people/dennis/callisto/introduction.html">introduction</A> 
   and take a look at the <A href="http://www.cs.uu.nl/people/dennis/callisto/using.html">short example</A>.
 */

/*!
  \defgroup Errors Error return codes in Callisto
*/

/*!
  \defgroup Basic Basic functions
*/
/*!
  \defgroup Advanced Advanced functions
*/

/*!
  \defgroup Global Global functions
  \ingroup Basic
*/
/*!
  \defgroup Objectdef Object definition
  \ingroup Basic
*/
/*!
  \defgroup Group Group manipulation
  \ingroup Basic
*/
/*!
  \defgroup Objectman Object manipulation
  \ingroup Basic
*/
/*!
  \defgroup Collision Collision detection
  \ingroup Basic
*/

/*!
  \defgroup AdvancedGlobal Global functions
  \ingroup Advanced
*/
/*!
  \defgroup AdvancedObject Object functions
  \ingroup Advanced
*/
/*!
  \defgroup AdvancedGroup Group functions
  \ingroup Advanced
*/
/*!
  \defgroup Retrieval Group/Object retrieval
  \ingroup Advanced
*/

/*!
  \defgroup Animation Animation functions
*/

/*!
  \defgroup Statistics Statistical information
*/

/*! CAL_scalar value used. */
typedef float CAL_scalar;

/*! CAL_matrix value used. */
typedef CAL_scalar CAL_matrix[4][4];

/*! Callback function type for key press. */
typedef void (* CAL_KeypressCallback) (char key, bool pressed);

/*! Callback function type for object selection. */
typedef void (* CAL_ObjectSelectCallback) (int objID, CAL_scalar hitPoint[3], void *data);

/*! Callback function type for a user drawn object \ingroup AdvancedObject */
typedef void (*CAL_UserDrawnCallback)(void *userdef, float time);

/*! CAL_NULL value. */
#define CAL_NULL 0

/*! The function was succesful. \ingroup Errors */
#define CAL_SUCCESS                        0
/*! Callistio is not yet initialized, use CAL_Initialisation. \ingroup Errors */
#define CAL_NOT_INIT                      -1
/*! The group does not exist. \ingroup Errors */
#define CAL_NO_SUCH_GROUP                 -2
/*! The parent group does not exist. \ingroup Errors */
#define CAL_NO_SUCH_PARENT_GROUP          -3
/*! The root group is a system group and cannot be manipulated. \ingroup Errors */
#define CAL_CANNOT_MANIPULATE_ROOT_GROUP  -4
/*! The object does not exist. \ingroup Errors */
#define CAL_NO_SUCH_OBJECT                -5
/*! One of the parameters has an illegal value. \ingroup Errors */
#define CAL_ILLEGAL_VALUE                 -6
/*! Cannot collision check because both groups are part of the same tree-branch. \ingroup Errors */
#define CAL_GROUPS_IN_SAME_SUBTREE        -7
/*! This label does not exist. \ingroup Errors */
#define CAL_LABEL_NOT_FOUND               -8
/*! This label already exists. \ingroup Errors */
#define CAL_LABEL_ALREADY_EXISTS          -9
/*! The file you try to read is not valid or could not be found. \ingroup Errors */
#define CAL_FILE_ERROR                   -10
/*! You try to do something with the visualiser while it is not running or it is suspended. \ingroup Errors */
#define CAL_VIS_NOT_RUNNING              -11
/*! The clone you are trying to make cannot be put in the source's subtree. \ingroup Errors */
#define CAL_CLONE_IN_SUBGROUP            -12
/*! You try to do something with a view, while it is not visible. \ingroup Errors */
#define CAL_VIEW_NOT_VISIBLE             -13
/*! Cannot calculate penetration depth, because objects do not overlap. \ingroup Errors */
#define CAL_GROUPS_DO_NOT_OVERLAP        -14
/*! Group not capable of collision checks. \ingroup Errors */
#define CAL_GROUP_NOT_COL_CAPABLE        -15
/*! Cannot load texture, check filename/path and fileformat. \ingroup Errors */
#define CAL_TEXTURE_ERROR                -16
/*! Key state already defined. \ingroup Errors */
#define CAL_ILLEGAL_KEY_STATE            -17
/*! You try to set properties of a group/object that has a motion defined. \ingroup Errors */
#define CAL_IS_DYNAMIC                   -18
/*! You try to use a file with the incorrect extension. \ingroup Errors */
#define CAL_INVALID_EXTENSION            -19
/*! Either provide visibility for all keystates of a group, or none. \ingroup Errors */
#define CAL_CANNOT_SET_VISIBILITY        -20
/*! Statistical information is not enabled. \ingroup Errors */
#define CAL_STATISTICSNOTENABLED         -21

/*! Show the navigation compass. \ingroup AdvancedGlobal */
#define CAL_SHOWCOMPASS          1
/*! Hide the navigation compass. \ingroup AdvancedGlobal */
#define CAL_HIDECOMPASS          2
/*! Show the grid. \ingroup AdvancedGlobal */
#define CAL_SHOWGRID             4
/*! Hide the grid. \ingroup AdvancedGlobal */
#define CAL_HIDEGRID             8
/*! Use perspective projection (default) in visualisation window. \ingroup AdvancedGlobal */
#define CAL_PERSPPROJ           16
/*! Use orthogonal projection in visualisation window.. \ingroup AdvancedGlobal */
#define CAL_ORTHOPROJ           32
/*! Show the position of the camera in the bottom right corner of the screen. \ingroup AdvancedGlobal */
#define CAL_SHOWSTATUSTEXT      64
/*! Hide the position of the camera in the bottom right corner of the screen. \ingroup AdvancedGlobal */
#define CAL_HIDESTATUSTEXT     128
/*! Turn on backface culling (hide polygons that have their normal pointed away from the eye). \ingroup AdvancedGlobal */
#define CAL_BACKFACECULLINGON  256
/*! Turn off backface culling. \ingroup AdvancedGlobal */
#define CAL_BACKFACECULLINGOFF 512
/*! Prevent the user from navigating using the mouse and/or keyboard. \ingroup AdvancedGlobal */
#define CAL_LOCKNAVIGATION    1024
/*! Allow the user to navigate. \ingroup AdvancedGlobal */
#define CAL_UNLOCKNAVIGATION  2048
/*! Disables keyboard response of Callisto navigation. \ingroup AdvancedGlobal */
#define CAL_LOCKKEYBOARDNAVIGATION 4096
/*! Enables keyboard response of Callisto navigation. \ingroup AdvancedGlobal */
#define CAL_UNLOCKKEYBOARDNAVIGATION 8192
/*! Enable the animation of the viewpoint. \ingroup AdvancedGlobal */
#define CAL_ENABLEANIMATION   16384
/*! Disable the animation of the viewpoint (viewKeystates are ignored). \ingroup AdvancedGlobal */
#define CAL_DISABLEANIMATION  32768

/*! Enable statistical information of collision checks. \ingroup AdvancedGlobal */
#define CAL_ENABLESTATISTICS  1
/*! Disable statistical information of collision checks. \ingroup AdvancedGlobal */
#define CAL_DISABLESTATISTICS 2

/*! Use with CAL_SetGroupMotionOptions or CAL_SetObjectMotionOptions. Makes the group motion cyclic. \ingroup AdvancedGroup */
#define CAL_CYCLIC           1
/*! Use with CAL_SetGroupMotionOptions or CAL_SetObjectMotionOptions. Makes the group motion non-cyclic. \ingroup AdvancedGroup */
#define CAL_NONCYCLIC        2
/*! Use with CAL_SetGroupMotionOptions or CAL_SetObjectMotionOptions. Interpolates between keystates (using linear interpolation for translation and SLERP for quaternion interpolation. \ingroup AdvancedGroup */
#define CAL_INTERPOLATION    4
/*! Use with CAL_SetGroupMotionOptions or CAL_SetObjectMotionOptions. Does not interpolate keystates (group/object jumps). \ingroup AdvancedGroup */
#define CAL_NOINTERPOLATION  8

#define CAL_FALSE 0
#define CAL_TRUE 1
#define CAL_USEPARENT 2

/*! This function initializes Callisto, starts the output window and GUI.
	\param visualisation Set to false if you dont want a visualisation window, default is TRUE.
  \param guiVisible Set to false if you dont want a GUI, only applicable if visualisation=true.
  \param extendedInterface Set to true if you want the extended interface.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_Initialisation (bool visualisation=true, bool guiVisible=true, bool extendedInterface=false);

/*! This function ends Callisto and cleans up memory.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_End ();

/*! This function stops the visualisation until CAL_ResumeVisualisation is called.
  This function can be used to temporarily stop visualisation, for example when you need all processor power to do some calculations.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SuspendVisualisation ();

/*! This function resumes the visualisation after CAL_SuspendVisualisation.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_ResumeVisualisation ();

/*! Shows a view. Every view has its own unique ID.
    There are at most 4 views (0..3). 0 being the main view. This view cannot be switched on/off.
  \param viewID The ID of the view.
  \param title The title of the window of the view, default is no title.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_ShowView (int viewID, char* title="");

/*! Hides a view. Every view has its own unique ID.
    There are at most 4 views (0..3). 0 being the main view. This view cannot be switched on/off.
  \param viewID The ID of the view.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_HideView (int viewID);

/*! Change the view parameters.
  \param viewID The ID of the view.
  \param options Legal values are #CAL_SHOWCOMPASS, #CAL_SHOWGRID etc. Multiple parameters can be changed by using the | operator.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_SetViewOptions (int viewID, int options);

/*! This function sets the position and orientation of the camera in the visualisation window.
    It is positioned at eye and looks to the point defined with view. By default the up vector is (0,1,0)
    You can change this if necessary.
  \param viewID The ID of the view to set the view parameters for, set to 0 if you don't know what this is.
  \param eyeX The x position of the camera.
  \param eyeY The y position of the camera.
  \param eyeZ The z position of the camera.
  \param viewX The x position of the look at point.
  \param viewY The y position of the look at point.
  \param viewZ The z position of the look at point.
  \param upX The x direction of the up vector, default is 0.
  \param upY The y direction of the up vector, default is 1.
  \param upZ The z direction of the up vector, default is 0.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SetViewParams (int viewID, CAL_scalar eyeX, CAL_scalar eyeY, CAL_scalar eyeZ, CAL_scalar viewX, CAL_scalar viewY, CAL_scalar viewZ, CAL_scalar upX=0, CAL_scalar upY=1, CAL_scalar upZ=0);

/*! This function gets the current position and orientation of the camera in a visualisation window.
    It is positioned at eye and looks to the point defined with view. 
  \param viewID The ID of the view get the view parameters from.
  \param *eyeX The x position of the camera.
  \param *eyeY The y position of the camera.
  \param *eyeZ The z position of the camera.
  \param *viewX The x position of the look at point.
  \param *viewY The y position of the look at point.
  \param *viewZ The z position of the look at point.
  \param *upX The x direction of the up vector, this is optional.
  \param *upY The y direction of the up vector, this is optional.
  \param *upZ The z direction of the up vector, this is optional.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_GetViewParams (int viewID, float *eyeX, float *eyeY, float *eyeZ, float *viewX, float *viewY, float *viewZ, float *upX=CAL_NULL, float *upY=CAL_NULL, float *upZ=CAL_NULL);

/*! Adds a keystate for the view animation.
  \param viewID The ID of the view.
  \param time The time stamp for which the configuration accounts.
  \param eye Three scalar values that determine the position of the eye point.
  \param view Three scalar values that determine the position of the look-at point.
  \param up Three scalar values that determine the up vector.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_AddViewKeyState (int viewID, float time, CAL_scalar *eye=CAL_NULL, CAL_scalar *view=CAL_NULL, CAL_scalar *up=CAL_NULL);

/*! Determines whether a view animation is cyclic or not.
  \param viewID The ID of the view.
  \param cyclic Set to true for a cyclic motion, or false for a non-cyclic motion.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetViewCyclic (int viewID, bool cyclic);

/*! Deletes all existing view key states.
  \param viewID The ID of the view.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_ClearViewKeyStates (int viewID);

/*! Sets the background color of the visualisation window.
  \param viewID The ID of the view get the view parameters from, set to 0 if you don't know what this is.
  \param r The red component of the color (0...1).
  \param g The green component of the color (0...1).
  \param b The blue component of the color (0...1).
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SetBackgroundColor (int viewID, float r, float g, float b);

/*! Loads a texture to memory. The texture must be in the .ppm or the .png format. Not supported by CAL_Polygon, CAL_PolygonGroup and CAL_Triangles.
  \param textureID The ID you want to give the texture, there is room for 500 textures numbered (0..499).
  \param fileName The filename+path of the texture. The texture has to be in .ppm or .png format.
  \param alpha The alpha value of the texture.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_LoadTexture (int textureID, char* fileName, float alpha);

/*! Add a texture from memory.
  \param textureID The ID you want to give the texture, there is room for 500 textures numbered (0..499).
  \param width The width in pixels of the texture. Must be a power of 2.
  \param height The height in pixels of the texture. Must be a power of 2.
  \param tex The textures in RGB format.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_SetTextureFromMem (int textureID, int width, int height, unsigned char *tex);

/*! Saves the content of a view in .bmp format to disk.
  \param viewID The ID of the view to save, set to 0 if you don't know what this is.
  \param fileName The file name of the .bmp.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_ScreenCapture (int viewID, char* fileName);

/*! Load a scene in XML or VRML format from disk.
  \param *fileName The name of the file to load.
  \param parentID The parent group to put the loaded scene in, use 0 for no parent.
  \param **error String with possible error string. Can be omitted.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_LoadScene (char* fileName, int parentID, char **error=CAL_NULL);

/*! Save a (part of a) scene in XML or VRML format from disk.
  \param *fileName The name of the file to load.
  \param viewID The view id which needs to be saved (only visible objects are saved)
  \param groupID The group id of the group that needs to be saved
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Global
*/
int CAL_SaveScene (char* fileName, int viewID=0, int groupID=0);

/*! Change the time, all dynamic groups and objects will adapt to this time.
  \param time The time;
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetTime (float time);

/*! Creates an object group.
  \param groupID This is set to the group ID the group gets.
  \param parentID The parent group of this group, use 0 if this group does not need a parent.
  \param colCheck True if this group is used for collision checks, false if not.
  \param label The label of the group as shown in the GUI, default is no label.
  \param collapsed The group will appear collapsed in the interface.
  \param visible Set to either CAL_TRUE, CAL_FALSE or CAL_USEPARENT.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_CreateGroup (int* groupID, int parentID, bool colCheck, char* label="", bool collapsed=false, int visible=CAL_USEPARENT);

/*! Delets a group, its child groups and all its objects.
  \param groupID The ID of the group to delete.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_DestroyGroup (int groupID);

/*! Deletes a groups child groups and all its objects.
  \param groupID The ID of the group to empty.
  \param subGroups Set to TRUE if subgroups of group should also be removed, default is FALSE.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_EmptyGroup (int groupID, bool subGroups=false);

/*! Move a group to a new parent.
  \param groupID The ID of the group to move.
  \param parentID The ID of the new parent.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_MoveGroup (int groupID, int parentID);

/*! Place a group and all its objects at a new position.
  \param groupID The ID of the group to translate.
  \param x The x component of the new position.
  \param y The y component of the new position.
  \param z The z component of the new position.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_SetGroupPosition (int groupID, CAL_scalar x, CAL_scalar y, CAL_scalar z);

/*! Place a group and all its objects in a new orientation using Euler angles (in radians).
  \param groupID The ID of the group to translate.
  \param xRot The rotation orientation with respect to the x-axis;
  \param yRot The rotation orientation with respect to the y-axis;
  \param zRot The rotation orientation with respect to the z-axis;
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_SetGroupOrientation (int groupID, CAL_scalar xRot, CAL_scalar yRot, CAL_scalar zRot);

/*! Scale a group and all its objects.
  \param groupID The ID of the group to scale.
  \param xScale The scaling factor in the x direction.
  \param yScale The scaling factor in the y direction.
  \param zScale The scaling in the z direction.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_SetGroupScaling (int groupID, CAL_scalar xScale, CAL_scalar yScale, CAL_scalar zScale);

/*! Rotate a group and all its objects using a quaternion.
  \param groupID The ID of the group to rotate.
  \param x The x component of the quaternion.
  \param y The y component of the quaternion.
  \param z The z component of the quaternion.
  \param w The w component of the quaternion.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_SetGroupQuaternion (int groupID, CAL_scalar x, CAL_scalar y, CAL_scalar z, CAL_scalar w);

/*! Spherically expand a group. Collision will occur with the set of points whose distance to the objects is at most the clearance.
  \param groupID The ID of the group to translate.
  \param c The clearance.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Group
*/
int CAL_SetGroupClearance (int groupID, CAL_scalar c);

/*! Sets the color of a group.
  \param groupID The ID of the group.
  \param r The red component of the color (0...1).
  \param g The green component of the color (0...1).
  \param b The blue component of the color (0...1).
  \param a The alpha value of the color (0...1), 0 is fully transparant, 1 is fully opaque. Default is 1.
  \param subGroups Set to TRUE if subgroups of group should also get new color, default is FALSE.
  \param sID The ID of the surface which gets this color, default is surface 0.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_SetGroupColor (int groupID, float r, float g, float b, float a=1, bool subGroups=false, int sID=0);

/*! Set the texture for a group.
  \param groupID The ID of the group.
  \param textureID The ID of the texture.
  \param xtile The number of times to repeat the texture in the x-direction.
  \param ytile The number of times to repeat the texture in the y-direction.
  \param subGroups Set to TRUE if subgroups of group should also get new color. Default is FALSE.
  \param sID The ID of the surface which gets this texture, default is surface 0.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_SetGroupTexture (int groupID, int textureID, float xtile, float ytile, bool subGroups=false, int sID=0);

/*! Sets the active surface of the group.
  \param groupID The ID of the group
  \param subGroups Set to true if subgroups also have to change the active surface.
  \param sID The ID of the new active surface.
  \ingroup AdvancedGroup
*/
int CAL_SetGroupActiveSurface (int groupID, bool subGroups, int sID);

/*! Sets the collision check capability of a group.
  \param groupID The ID of the group.
  \param colCapable True if the group can be collision checked.
  \param subGroups Set to true if subgroups also have to change their collision check capability recursively.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_SetGroupCollisionCheckCapability (int groupID, bool colCapable, bool subGroups=false);

/*! Sets the visibility of a group in a particular view. This does not effect collision checks.
  \param groupID The ID of the group.
  \param viewID The ID of the view.
  \param visible True if the group should be visible. False for invisibility.
  \param subGroups Set to true if subgroups also have to change their visibility recursively.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_SetGroupVisibility (int groupID, int viewID, bool visible, bool subGroups=false);

/*! Change the label of a group.
  \param groupID The ID of the group.
  \param label The new label
  \ingroup AdvancedGroup
*/
int CAL_SetGroupLabel (int groupID, char *label);

/*! Adds a keystate for the group animation.
  \param groupID The ID of the group.
  \param time The time stamp for which the configuration accounts.
  \param position Three scalar values that determine the position of the group at this keystate.
  \param orientation The orientation of the group, this can consist of Euler angles or a quaternion, depending on the last parameter.
  \param scaling The scaling of the group, consisting of three values.
  \param isQuat Set this to false if the second paramter consists of Euler angles.
  \param visible If set the main view group visibility adapts to this parameter, if not set the group value (by CAL_SetGroupVisibility) is used for the whole animation.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_AddGroupKeyState (int groupID, float time, CAL_scalar *position, CAL_scalar *orientation=CAL_NULL, CAL_scalar *scaling=CAL_NULL, bool isQuat=true, int visible=2);

/*! Change the group motion parameters.
  \param groupID The ID of the group.
  \param options Legal values are: CAL_CYCLIC, CAL_NONCYCLIC, CAL_INTERPOLATION, CAL_NOINTERPOLATION.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetGroupMotionOptions (int groupID, int options);

/*! Deletes all existing group key states.
  \param groupID The ID of the group.
  \param subGroups Set to TRUE if subgroups of group should also be cleared.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_ClearGroupKeyStates (int groupID, bool subGroups=false);

/*! Clones a group including all objects.
  \param groupIDNew Set to the group ID of the clone.
  \param groupID The ID of the group to clone.
  \param parentID The parent group of the clone (0 for no parent).
  \param colCheck Set to false if clone group does not need collision checks.
  \param label The label of the group, default is no label.
  \param cloneObjs Optional flag that indicates whether to clone the objects in the group as well.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_CloneGroup (int* groupIDNew, int groupID, int parentID, bool colCheck, char* label="", bool cloneObjs=true);

/*! Clones a group including all objects and child groups.
  \param newgroupIDs List of new group IDs, can be NULL.
  \param groupID The ID of the group to clone.
  \param parentID The parent group of the clone (0 for no parent).
  \param nr The size of the ids and labels lists, can be 0.
  \param ids List of id's corresponding to the labels list, can be NULL.
  \param labels Optional list of new labels for the clone, defaults to NULL.
  \param cloneObjs Optional flag that indicates whether to clone the objects in the group as well, defaults to true.
  \param keepColCap Optional flag that indicates whether the collision capabilities should be preserved. Default clones do not have collision check capabilities.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGroup
*/
int CAL_CloneGroupRecursive (int* newgroupIDs, int groupID, int parentID, int nr, int* ids, char** labels=CAL_NULL, bool cloneObjs=true, bool keepColCap=false);

/*! Move an object to another group.
  \param objID The ID of the object to destroy.
  \param groupID The ID of the new group.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedObject
*/
int CAL_MoveObject (int objID, int groupID);

/*! Destroys an object.
  \param objID The ID of the object to destroy.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_DestroyObject (int objID);

/*! Place an object at a new position.
  \param objID The ID of the object to translate.
  \param x The x component of the new position.
  \param y The y component of the new position.
  \param z The z component of the new position.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectPosition (int objID, CAL_scalar x, CAL_scalar y, CAL_scalar z);

/*! Place an object in a new orientation using Euler angles (in radians).
  \param objID The ID of the object to translate.
  \param xRot The orientation with respect to the x-axis;
  \param yRot The orientation with respect to the y-axis;
  \param zRot The orientation with respect to the z-axis;
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectOrientation (int objID, CAL_scalar xRot, CAL_scalar yRot, CAL_scalar zRot);

/*! Scale an object.
  \param objID The ID of the object to scale.
  \param xScale The scaling factor in the x direction.
  \param yScale The scaling factor in the y direction.
  \param zScale The scaling factor in the z direction.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectScaling (int objID, CAL_scalar xScale, CAL_scalar yScale, CAL_scalar zScale);

/*! Rotate an object using a quaternion.
  \param objID The ID of the object to rotate.
  \param x The x component of the quaternion.
  \param y The y component of the quaternion.
  \param z The z component of the quaternion.
  \param w The w component of the quaternion.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectQuaternion (int objID, CAL_scalar x, CAL_scalar y, CAL_scalar z, CAL_scalar w);

/*! Set the WORLD matrix of an object.
  \param objID The ID of the object.
  \param matrix Should be of type CAL_matrix;
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectWorldMatrix (int objID, CAL_matrix *matrix);

/*! Spherically expand an object. Collision will occur with the set of points whose distance to the object is at most the clearance.
  \param objID The ID of the object to translate.
  \param c The clearance.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectClearance (int objID, CAL_scalar c);

/*! Sets the color of an object.
  \param objID The ID of the object.
  \param r The red component of the color (0...1).
  \param g The green component of the color (0...1).
  \param b The blue component of the color (0...1).
  \param a The alpha value of the color (0...1), 0 is fully transparant, 1 is fully opaque. Default is 1.
  \param sID The ID of the surface which gets this color. Default value is surface 0.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectman
*/
int CAL_SetObjectColor (int objID, float r, float g, float b, float a=1, int sID=0);

/*! Set the texture for an object.
  \param objID The ID of the object.
  \param textureID The ID of the texture.
  \param xtile The number of times to repeat the texture in the x-direction.
  \param ytile The number of times to repeat the texture in the y-direction.
  \param sID The ID of the surface which gets this texture, default is surface 0.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedObject
*/
int CAL_SetObjectTexture (int objID, int textureID, float xtile, float ytile, int sID=0);

/*! Sets the active surface of the group.
  \param objID The ID of the object.
  \param sID The ID of the new active surface.
  \ingroup AdvancedObject
*/
int CAL_SetObjectActiveSurface (int objID, int sID);

/*! Adds a keystate for the object animation.
  \param objID The ID of the object.
  \param time The time stamp for which the configuration accounts.
  \param position Three scalar values that determine the position of the object at this keystate.
  \param orientation The orientation of the object, this can consist of Euler angles or a quaternion, depending on the last parameter.
  \param scaling The scaling of the group, consisting of three values.
  \param isQuat Set this to false if the second paramter consists of Euler angles.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_AddObjectKeyState (int objID, float time, CAL_scalar *position, CAL_scalar *orientation=CAL_NULL, CAL_scalar *scaling=CAL_NULL, bool isQuat=true);

/*! Change the objectmotion parameters.
  \param objID The ID of the object.
  \param options Legal values are: CAL_CYCLIC, CAL_NONCYCLIC, CAL_INTERPOLATION, CAL_NOINTERPOLATION.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_SetObjectMotionOptions (int objID, int options);

/*! Deletes all existing object key states.
  \param objID The ID of the object.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Animation
*/
int CAL_ClearObjectKeyStates (int objID);

/*! Returnes the ID of a group or object with a certain label.
  \param ID This is set to the ID.
  \param label The label of the group/object.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetID (int* ID, char* label);

/*! Get the properties of a group.
  \param groupID The ID of the group.
  \param CALGroup Pointer to an SCALGroup-structure.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetGroup (int groupID, void *CALGroup);

/*! Returnes the ID of the nr'th childgroup.
  \param groupID The ID of the group.
  \param nr The nr of the childgroup.
  \param childGroupID This value is set to the nr'th childgroup.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetGroupChildID (int groupID, int nr, int *childGroupID);

/*! Returnes the ID of the nr'th object.
  \param groupID The ID of the group.
  \param nr The nr of the object.
  \param objectID This value is set to the nr'th object.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetGroupObjectID (int groupID, int nr, int *objectID);


/*! Get the type of an object (#CAL_BOX, #CAL_CYLINDER etc.).
  \param objID The ID of the object.
  \param objType This will be set to the object type according to the values in callistoTypes.h.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetObjectType (int objID, int* objType);

/*! Get the WORLD matrix of an object.
  \param objID The ID of the object.
  \param matrix Should be of type CAL_matrix;
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetObjectWorldMatrix (int objID, CAL_matrix *matrix);

/*! Get the properties of an object.
  \param objID The ID of the object.
  \param SCALObj Pointer to an SCAL-object (SCALBox, SCALSphere etc.). This has to be of the right type (the type can be retrieved by using CAL_GetObjectType).
         Note that to retrieve a CAL_ELEVATIONGRID object, SCALPolygonGroup needs to be used.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Retrieval
*/
int CAL_GetObject (int objID, void* SCALObj);

/*! Set a callback function, this is called when the user pressed a key.
  \param cb The adress of the callback function.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedObject
*/
int CAL_SetKeypressCallback (CAL_KeypressCallback cb);

/*! Set a callback function, this is called when the user selects an object by using the mouse and pressing SPACE.
  \param cb The adress of the callback function.
  \param data A user definable data pointer that is supplied as the last argument of the callback function
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedObject
*/
int CAL_SetObjectSelectCallback (CAL_ObjectSelectCallback cb, void* userdef = NULL);

/*! Create a box.
  \param groupID The group ID to put the object in.
  \param xw The width of the object.
  \param yw The height of the object.
  \param zw The depth of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateBox (int groupID, CAL_scalar xw, CAL_scalar yw, CAL_scalar zw, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* label="");

/*! Create a sphere.
  \param groupID The group ID to put the object in.
  \param r The radius of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateSphere (int groupID, CAL_scalar r, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* label="");

/*! Create a cylinder.
  \param groupID The group ID to put the object in.
  \param r The radius of the object.
  \param h The height of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateCylinder (int groupID, CAL_scalar r, CAL_scalar h, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* label="");

/*! Create a sphere.
  \param groupID The group ID to put the object in.
  \param r The radius of the bottom of the object.
  \param h The height of the object.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateCone (int groupID, CAL_scalar r, CAL_scalar h, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* label="");

/*! Create a triangle mesh.
  \param groupID The group ID to put the object in.
  \param nt The number of triangles the object consists of.
  \param *p List of coordinates of the object. The size of the list must be 9*nt.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateTriangles (int groupID, int nt, CAL_scalar* p, int* objID=CAL_NULL, char* label="");

/*! Create a polygon. Must be planar and convex.
  \param groupID The group ID to put the object in.
  \param np The number of points the object consists of.
  \param *p List of coordinates of the object. The size of the list must be 3*np.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreatePolygon (int groupID, int np, CAL_scalar* p, int* objID=CAL_NULL, char* label="");

/*! Create a polygon group. All must be planar and convex.
  \param groupID The group ID to put the object in.
  \param npol The number of polygons the object consists of.
  \param npoints The number of points each polygon consists of.
  \param *points List of coordinates of all polygons, its length should be npol * npoints * 3.
  \param objID Set to the object ID.
  \param label Label of the object, default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreatePolygonGroup (int groupID, int npol, int *npoints, CAL_scalar *points, int* objID=CAL_NULL, char* label="");

/*! Create a polyline. This is for the visualisation only.
  \param groupID The group ID to put the object in.
  \param nl The number of lines.
  \param *np The number of points each line consists of.
  \param *p List of coordinates of the lines. Its size should be is 3 * np * nl.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreatePolyline (int groupID, int nl, int *np, CAL_scalar *p, int* objID=CAL_NULL, char* label="");

/*! Create a tetrahedron.
  \param groupID The group ID to put the object in.
  \param *p List of coordinates of the object. The size of the list must be 3*4 since a tetrahydron consists of 4 points.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateTetrahedron (int groupID, CAL_scalar* p, int* objID=CAL_NULL, char* label="");

/*! Create an elevation grid on the XZ plane.
  \param groupID The group ID to put the object in.
  \param xDim The number of x coordinates.
  \param zDim The number of z coordinates.
  \param xStep The stepsize along the x-axis.
  \param zStep The stepsize along the z-axis.
  \param p The height parameters. The length of p should be (xDim+1)*(zDim+1).
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Objectdef
*/
int CAL_CreateElevationGrid (int groupID, int xDim, int zDim, float xStep, float zStep, CAL_scalar* p, int *objID=CAL_NULL, char* label="");

/*! Create a user drawn object. This object is only useable in the visualisation, and not for collision detection.
  \param groupID The group ID to put the object in.
  \param callback The user drawn callback. This is called from inside the render loop. The user can execute arbitrary openGL code.
  \param userdef A user defined pointer that is passed to the callback.
  \param x The x position of the object.
  \param y The y position of the object.
  \param z The z position of the object.
  \param objID Set to the object ID.
  \param label Label of the object. Default is no label.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedObject
*/
int CAL_CreateUserDrawn (int groupID, CAL_UserDrawnCallback callback, void *userdef, CAL_scalar x, CAL_scalar y, CAL_scalar z, int* objID=CAL_NULL, char* label="");

/*! Check whether a point collides with a group.
  \param groupID The group ID to check with.
  \param x The x coordinate of the point to check.
  \param y The y coordinate of the point to check.
  \param z The z coordinate of the point to check.
  \param multiple Flag whether to find all results, or just the first (faster).
  \param *nrCols The number of collisions.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Collision
*/
int CAL_CheckPointCollision (int groupID, CAL_scalar x, CAL_scalar y, CAL_scalar z, bool multiple, int *nrCols);

/*! Check whether a line collides with a group.
  \param groupID The group ID to check with.
  \param x0 The x coordinate of the first point to check.
  \param y0 The y coordinate of the first point to check.
  \param z0 The z coordinate of the first point to check.
  \param x1 The x coordinate of the second point to check.
  \param y1 The y coordinate of the second point to check.
  \param z1 The z coordinate of the second point to check.
  \param multiple Flag whether to find all results, or just the first (faster).
  \param *nrCols The number of collisions.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Collision
*/
int CAL_CheckLineCollision (int groupID, CAL_scalar x0, CAL_scalar y0, CAL_scalar z0, CAL_scalar x1, CAL_scalar y1, CAL_scalar z1, bool multiple, int *nrCols);

/*! Check whether two groups collide. Groups cannot be each others subgroups.
  \param group0 The ID of the first group.
  \param group1 The ID of the second group.
  \param multiple Flag whether to find all results, or just the first (faster).
  \param *nrCols The number of collisions.
	\returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup Collision
*/
int CAL_CheckGroupCollision (int group0, int group1, bool multiple, int *nrCols);

/*! Get the two positions where the distance between two groups is smallest.
  \param groupID0 The ID of the first group.
  \param groupID1 The ID of the second group.
  \param nrPairs The number of results.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_GetClosestPairs (int groupID0, int groupID1, int *nrPairs);

/*! Get the two positions where the penetration of two groups is largest.
  \param groupID0 The ID of the first group.
  \param groupID1 The ID of the second group.
  \param nrPairs The number of results.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_GetPenetrationDepths (int groupID0, int groupID1, int *nrPairs);

/*! Get the results of the objects involved in the last collision check/penetration depth/closest pair.
  \param userResults A list of type SCALResult which contains the results. The client is reponsible of creating the list with the right size (size == count). The results for closest pairs and penetration depths are sorted on distance.
  \returns The function returns #CAL_SUCCESS on success, and an \ref Errors "errorcode" on failure.
  \ingroup AdvancedGlobal
*/
int CAL_GetResults (void* userResults);

/*! Enable statistics. Statistical information will be gathered for collision checks, closest pairs and penetration depth. Both the number of calls and the total time spent in these functions will be administered. Beware that to gather the statistics themselved a small amount of time is spent.
  \param enable This value has to be either equal to either CAL_ENABLESTATISTICS or CAL_DISABLESTATISTICS.
  \ingroup Statistics
*/
int CAL_GatherStatistics (int enable);

/*! Save the statistical information to a comma delimited file.
  \param groupID The ID of the group.
  \param fileName The path and name of the file to write the information to.
  \ingroup Statistics
*/
int CAL_SaveGroupStatistics (int groupID, char *fileName);

/*! \file callisto.h
\brief All functions and constants of Callisto.
*/