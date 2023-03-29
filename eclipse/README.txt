The GDSC-SMLM code was developed using the Eclipse IDE:
https://eclipse.org/

The code can be built using Maven. See the README.md for details. However using
Eclipse is preferred for code development to provide a debugging environment.

You will need the Maven and Git Eclipse plugins. The standard Eclipse IDE for
Java developers has these.

Set-up the project
------------------

Import the project into Eclipse (File>Import)
Select: Maven > Existing Maven projects

This will import the project but may not link it to the source Git repository.
Right-click on the project name and select 'Team > Share'. If you share it back to
the same location it will attach to the source Git repository.

Code formatting
---------------

The Eclipse code format rules (in this directory) can be loaded using:

Eclipse Preferences : Java > Code Style > Formatter
Eclipse Preferences : Java > Code Style > Clean Up

Click 'Import...' to load the provided rules.

Running the code
----------------

Build the project.

---
Optional:

To debug the ImageJ 3D Viewer from Eclipse relies on being able to find the native
runtime libraries for Open GL. These should now be included by Maven. If this fails
then do a manual install:

1. Download JOGL
   http://jogamp.org/wiki/index.php/Downloading_and_installing_JOGL
2. Extract all the files
3. Copy jogl-all.jar and gluegen-rt.jar to a folder
4. Copy the 'natives' jars for the current platform to the same location,
   e.g. jogl-all-natives-linux-amd64.jar and gluegen-rt-natives-linux-amd64.jar
5. Add the jogl-all.jar and gluegen-rt.jar as external jars to the build path in the
   Eclipse project
---

To run the code execute the 'main' function in the Smlm_Plugin class within the gdsc-smlm-ij_
project. This will:

- Configure the ImageJ plugins path to the build path
- Run ImageJ
- Register all the GDSC SMLM plugins in ImageJ and show a GDSC SMLM plugins window


Create a new Run configuration in the gdsc-smlm-ij_ project.

Select Smlm_Plugin as the main class.

Run the code.

The same target can be used to debug the code.
